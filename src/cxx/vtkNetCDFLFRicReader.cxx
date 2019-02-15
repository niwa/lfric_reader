#include "vtkNetCDFLFRicReader.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>

#include <unordered_map>

#define CALL_NETCDF(call) \
  { \
    int errorcode = call; \
    if (errorcode != NC_NOERR) \
    { \
      vtkErrorMacro(<< "NetCDF Error: " << nc_strerror(errorcode)); \
    } \
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkNetCDFLFRicReader)

//----------------------------------------------------------------------------
vtkNetCDFLFRicReader::vtkNetCDFLFRicReader()
{
  if (DEBUG)
  {
    this->DebugOn();
  }

  vtkDebugMacro("Entering vtkNetCDFLFRicReader constructor..." << endl);

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->FileName = nullptr;
  this->UseCartCoords = 0;
  this->Fields.clear();
  this->TimeSteps.clear();
  this->NumberOfLevels = 0;
  this->NumberOfFaces2D = 0;
  this->NumberOfEdges2D = 0;

  vtkDebugMacro("Finished vtkNetCDFLFRicReader constructor" << endl);
}

//----------------------------------------------------------------------------
vtkNetCDFLFRicReader::~vtkNetCDFLFRicReader()
{
  vtkDebugMacro("Entering vtkNetCDFLFRicReader destructor..." << endl);

  this->SetFileName(nullptr);
  this->Fields.clear();
  this->TimeSteps.clear();

  vtkDebugMacro("Finished vtkNetCDFLFRicReader destructor" << endl);
}

//----------------------------------------------------------------------------
void vtkNetCDFLFRicReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(nullptr)") << endl;

  os << indent << "NumberOfTimeSteps: "
     << this->TimeSteps.size() << endl;

  os << indent << "NumberOfFields: "
     << this->Fields.size() << endl;

  os << indent << "NumberOfCellArrays: "
     << this->GetNumberOfCellArrays() << endl;
}

//----------------------------------------------------------------------------

// Retrieve time steps and field names from netCDF file
// Convention is that this funtion returns 1 on success, 0 otherwise
int vtkNetCDFLFRicReader::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  vtkDebugMacro("Entering RequestInformation..." << endl);

  if(this->FileName == nullptr)
  {
    vtkErrorMacro("FileName not set.");
    return 0;
  }
  vtkDebugMacro("FileName=" << this->FileName << endl);

  // VTK object for handing over time step information
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int ncid;
  CALL_NETCDF(nc_open(this->FileName, NC_NOWRITE, &ncid));

  // time_counter can contain 0
  size_t NumberOfTimeSteps = getNCDim(ncid, "time_counter");
  vtkDebugMacro("Number of time steps in file=" << NumberOfTimeSteps << endl);

  // Read time steps array, if available, and tell the pipeline
  this->TimeSteps.clear();
  if (NumberOfTimeSteps > 0)
  {
    this->TimeSteps = getNCVarDouble(ncid, "time_instant",
                                     {0}, {NumberOfTimeSteps});

    // Tell the pipeline which steps are available and their range
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
		 this->TimeSteps.data(), static_cast<int>(this->TimeSteps.size()));

    double timeRange[2] = {this->TimeSteps.front(), this->TimeSteps.back()};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

    vtkDebugMacro("timeRange=" << timeRange[0] << " " << timeRange[1] << endl);
  }
  else
  {
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

    vtkDebugMacro("Only single time step available" << endl);
  }

  // Get variable names and populate "Fields" map, ignoring UGRID mesh definitions
  for (std::string name : getNCVarNames(ncid))
  {
    vtkDebugMacro("Considering variable name=" << name << endl);

    if (name.compare(0, 6, "Mesh2d") != 0 &&
        name.compare(0, 11, "half_levels") != 0 &&
        name.compare(0, 11, "full_levels") != 0 &&
        name.compare(0, 12, "time_instant") != 0)
    {
      // If field is not in list, insert and default to "don't load"
      std::map<std::string,bool>::iterator it = this->Fields.find(name);
      if (it == this->Fields.end())
      {
        this->Fields.insert(it, std::pair<std::string,bool>(name,false));
      }
    }
  }
  vtkDebugMacro("Number of data arrays found=" << this->Fields.size() << endl);

  CALL_NETCDF(nc_close(ncid));

  vtkDebugMacro("Finished RequestInformation" << endl);
  
  return 1;

}

//----------------------------------------------------------------------------
// Create VTK grid and load data for requested time step (defaults to first one)
// Convention is that this function returns 1 on success, or 0 otherwise
int vtkNetCDFLFRicReader::RequestData(vtkInformation *vtkNotUsed(request),
				      vtkInformationVector **vtkNotUsed(inputVector),
				      vtkInformationVector *outputVector)
{
  vtkDebugMacro("Entering RequestData..." << endl);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkUnstructuredGrid *outputGrid = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if(this->FileName == nullptr)
  {
    vtkErrorMacro("FileName not set.");
    return 0;
  }
  vtkDebugMacro("FileName=" << FileName << endl);

  // Default to first time step if no specific request has been made, or if
  // there is only a single time step in the file
  size_t timestep = 0;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
  {
    double requested_time = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    for (timestep = 0; timestep < this->TimeSteps.size(); timestep++)
    {
      if (this->TimeSteps[timestep] >= requested_time) break;
    }
    vtkDebugMacro("Requested time=" << requested_time << endl);
    vtkDebugMacro("Returning timestep=" << timestep << endl);
  }
  else
  {
    vtkDebugMacro("Defaulting to timestep=" << timestep << endl);
  }

  int ncid;
  CALL_NETCDF(nc_open(this->FileName, NC_NOWRITE, &ncid));

  // Read UGRID description from file and create unstructured VTK grid
  if (!this->CreateVTKGrid(ncid, outputGrid))
  {
    vtkErrorMacro("Could not create VTK grid.");
    return 0;
  }

  // Load requested field data for requested time step
  if (!this->LoadFields(ncid, outputGrid, timestep))
  {
    vtkErrorMacro("Could not load field data.");
    return 0;
  }

  CALL_NETCDF(nc_close(ncid));

  vtkDebugMacro("Finished RequestData" << endl);

  return 1;
}

size_t vtkNetCDFLFRicReader::getNCDim(const int ncid, const char * dimname)
{
  vtkDebugMacro("getNCDim: dimname=" << dimname << endl);
  int dimid;
  size_t dim;
  CALL_NETCDF(nc_inq_dimid(ncid, dimname, &dimid));
  CALL_NETCDF(nc_inq_dimlen(ncid, dimid, &dim));
  return dim;
}

std::vector<std::string> vtkNetCDFLFRicReader::getNCVarNames(const int ncid)
{
  std::vector<std::string> varNames;
  varNames.clear();

  int numVars;
  CALL_NETCDF(nc_inq_nvars(ncid, &numVars));
  vtkDebugMacro("getNCVarNames: Number of variables=" << numVars << endl);

  varNames.resize(numVars);

  for (int iVar = 0; iVar < numVars; iVar++)
  {
    char name[NC_MAX_NAME+1];
    CALL_NETCDF(nc_inq_varname(ncid, iVar, name));
    varNames[iVar] = name;
  }

  return varNames;
}

std::vector<double> vtkNetCDFLFRicReader::getNCVarDouble(const int ncid,
                    const char * varname,
                    const std::initializer_list<size_t> start,
                    const std::initializer_list<size_t> count)
{
  vtkDebugMacro("getNCVarDouble: varname=" << varname << endl);

  // Find variable by name
  int varId;
  CALL_NETCDF(nc_inq_varid(ncid, varname, &varId));

  // Compute total number of elements to read
  size_t size = 1;
  for (size_t n : count)
  {
    size *= n;
  }

  std::vector<double> varData;
  varData.resize(size);

  vtkDebugMacro("getNCVarDouble: reading size=" << size << " elements" << endl);

  // NetCDF will automatically convert non-double numeric data into double
  // This function will also check index ranges automatically
  // initialiser_list.begin() should give us access to a contiguous size_t array
  CALL_NETCDF(nc_get_vara_double(ncid, varId, start.begin(), count.begin(), varData.data()));

  return varData;
}

std::vector<unsigned long long> vtkNetCDFLFRicReader::getNCVarULongLong(const int ncid,
                                const char * varname,
                                const std::initializer_list<size_t> start,
                                const std::initializer_list<size_t> count)
{
  vtkDebugMacro("getNCVarULongLong: varname=" << varname << endl);

  int varId;
  CALL_NETCDF(nc_inq_varid(ncid, varname, &varId));

  size_t size = 1;
  for (size_t n : count)
  {
    size *= n;
  }

  std::vector<unsigned long long> varData;
  varData.resize(size);

  vtkDebugMacro("getNCVarULongLong: reading size=" << size << " elements" << endl);

  // netCDF will automatically convert non-double numeric data into unsigned long long
  CALL_NETCDF(nc_get_vara_ulonglong(ncid, varId, start.begin(), count.begin(), varData.data()));

  return varData;
}

void vtkNetCDFLFRicReader::mirror_points(vtkSmartPointer<vtkUnstructuredGrid> grid) {

  vtkDebugMacro("Entering mirror_points..." << endl);

  // Compute xy grid dimensions
  double gridBounds[6];
  grid->GetBounds(gridBounds);
  double gridDx = gridBounds[1] - gridBounds[0];
  double gridDy = gridBounds[3] - gridBounds[2];

  // Need to keep track of points that have already been duplicated,
  // to avoid degenerate points
  std::unordered_map<vtkIdType, vtkIdType> mirrorPointsX;
  std::unordered_map<vtkIdType, vtkIdType> mirrorPointsY;
  std::unordered_map<vtkIdType, vtkIdType> mirrorPointsXY;
  std::unordered_map<vtkIdType, vtkIdType>::const_iterator mirrorPointsIt;

  vtkPoints * gridPoints = grid->GetPoints();
  vtkSmartPointer<vtkIdList> oldCellPoints = vtkSmartPointer<vtkIdList>::New();

  // Search entire grid
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); cellId++)
  {

    // Compute xy cell dimensions
    double cellBounds[6];
    grid->GetCellBounds(cellId, cellBounds);
    double cellDx = cellBounds[1] - cellBounds[0];
    double cellDy = cellBounds[3] - cellBounds[2];

    // Find cells that span across the grid
    bool spanX = cellDx > 0.5*gridDx;
    bool spanY = cellDy > 0.5*gridDy;

    if (spanX or spanY)
    {

      grid->GetCellPoints(cellId, oldCellPoints);

      vtkIdType newCellPoints[8];

      // Check each cell vertex and mirror if needed
      for (vtkIdType pointIdIndex = 0; pointIdIndex < 8; pointIdIndex++)
      {

	vtkIdType thisPointId = oldCellPoints->GetId(pointIdIndex);
	double thisPointCoords[3];
	grid->GetPoint(thisPointId, thisPointCoords);

	// Mirror corner point
	if (spanX and spanY and thisPointCoords[0] < 0 and thisPointCoords[1] < 0)
        {
	  // Keep track of mirrored points to avoid degeneracy; insert a new point if
	  // no mirror point has been created yet
	  mirrorPointsIt = mirrorPointsXY.find(thisPointId);
	  if (mirrorPointsIt == mirrorPointsXY.end())
          {
            vtkIdType newPointId = gridPoints->InsertNextPoint(-thisPointCoords[0],
                                                               -thisPointCoords[1],
                                                                thisPointCoords[2]);
            mirrorPointsXY.insert({thisPointId, newPointId});
            newCellPoints[pointIdIndex] = newPointId;
	  }
          else
          {
            newCellPoints[pointIdIndex] = mirrorPointsIt->second;
          }
	}
	// Mirror point on left domain boundary
	else if (spanX && thisPointCoords[0] < 0)
        {
          mirrorPointsIt = mirrorPointsX.find(thisPointId);
	  if (mirrorPointsIt == mirrorPointsX.end())
          {
            vtkIdType newPointId = gridPoints->InsertNextPoint(-thisPointCoords[0],
                                                                thisPointCoords[1],
                                                                thisPointCoords[2]);
            mirrorPointsX.insert({thisPointId, newPointId});
            newCellPoints[pointIdIndex] = newPointId;
	  }
          else {
            newCellPoints[pointIdIndex] = mirrorPointsIt->second;
          }
	}
	// Mirror points on bottom domain boundary
	else if (spanY && thisPointCoords[1] < 0)
        {
          mirrorPointsIt = mirrorPointsY.find(thisPointId);
	  if (mirrorPointsIt == mirrorPointsY.end())
          {
            vtkIdType newPointId = gridPoints->InsertNextPoint(thisPointCoords[0],
                                                              -thisPointCoords[1],
                                                               thisPointCoords[2]);
            mirrorPointsY.insert({thisPointId, newPointId});
            newCellPoints[pointIdIndex] = newPointId;
	  }
          else
          {
            newCellPoints[pointIdIndex] = mirrorPointsIt->second;
          }
	}
	// No mirror point needed
	else
        {
	  newCellPoints[pointIdIndex] = thisPointId;
	}

      }
      grid->ReplaceCell(cellId, 8, newCellPoints);
    }
  }

  vtkDebugMacro("mirrorPointsX: " << mirrorPointsX.size() << endl);
  vtkDebugMacro("mirrorPointsY: " << mirrorPointsY.size() << endl);
  vtkDebugMacro("mirrorPointsXY: " << mirrorPointsXY.size() << endl);

  vtkDebugMacro("Finished mirror_points" << endl);
}

// Read UGRID description from netCDF file and build VTK grid
// The VTK grid replicates the "full level face grid" in the
// LFRic output file, data that is stored on the other grids
// is mapped onto this VTK grid
int vtkNetCDFLFRicReader::CreateVTKGrid(const int ncid, vtkUnstructuredGrid *grid)
{
  vtkDebugMacro("Entering CreateVTKGrid..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure is not available.");
    return 0;
  }

  // Read various dimensions, keep some of them in our object
  const size_t nnodes = getNCDim(ncid, "nMesh2d_full_levels_node");
  vtkDebugMacro("nnodes=" << nnodes << endl);

  this->NumberOfLevels = getNCDim(ncid, "full_levels");
  vtkDebugMacro("NumberOfLevels=" << this->NumberOfLevels << endl);

  this->NumberOfFaces2D = getNCDim(ncid, "nMesh2d_full_levels_face");
  vtkDebugMacro("NumberOfFaces2D=" << this->NumberOfFaces2D << endl);

  this->NumberOfEdges2D = getNCDim(ncid, "nMesh2d_full_levels_edge");
  vtkDebugMacro("NumberOfEdges2D=" << this->NumberOfEdges2D << endl);

  const size_t nverts_per_face = getNCDim(ncid, "nMesh2d_full_levels_vertex");
  vtkDebugMacro("nverts_per_face=" << nverts_per_face << endl);

  // x (lon) and y (lat) coordinates of grid nodes
  const std::vector<double> node_coords_x = getNCVarDouble(ncid, "Mesh2d_full_levels_node_x",
                                                          {0}, {nnodes});
  const std::vector<double> node_coords_y = getNCVarDouble(ncid, "Mesh2d_full_levels_node_y",
                                                          {0}, {nnodes});

  // Vertical level heights
  const std::vector<double> levels = getNCVarDouble(ncid, "full_levels",
                                                    {0}, {this->NumberOfLevels});

  // Node connectivity
  const std::vector<unsigned long long> face_nodes = getNCVarULongLong(ncid,
                                        "Mesh2d_full_levels_face_nodes",
                                        {0,0}, {this->NumberOfFaces2D,nverts_per_face});

  vtkDebugMacro("Setting VTK points..." << endl);

  // Construct VTK grid points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  if (this->UseCartCoords)
  {
    const double deg2rad = vtkMath::Pi()/180.0;
    for (size_t ilevel = 0; ilevel < this->NumberOfLevels; ilevel++)
    {
      for (size_t inode = 0; inode < nnodes; inode++)
      {
        // Avoid zero level height
        const double level_shift = 0.5;

        const double coords[3] = {node_coords_x[inode],
                                  node_coords_y[inode],
                                  levels[ilevel] + level_shift};

        // Convert from lon-lat-rad to Cartesian coordinates
        const double xyz[3] = {coords[2]*cos(coords[0]*deg2rad)*cos(coords[1]*deg2rad),
                               coords[2]*sin(coords[0]*deg2rad)*cos(coords[1]*deg2rad),
                               coords[2]*sin(coords[1]*deg2rad)};

        points->InsertNextPoint(xyz);
      }
    }
  }
  else
  {
    for (size_t ilevel = 0; ilevel < this->NumberOfLevels; ilevel++)
    {
      for (size_t inode = 0; inode < nnodes; inode++)
      {
        // Shift by 180 degrees to enable detection of grid periodicity
        const double coords[3] = {node_coords_x[inode]-180.0,
                                  node_coords_y[inode],
                                  levels[ilevel]};
        points->InsertNextPoint(coords);
      }
    }
  }
  grid->SetPoints(points);

  vtkDebugMacro("Setting VTK cells..." << endl);

  // Build up grid cells vertical layer-wise
  grid->Allocate(static_cast<vtkIdType>(this->NumberOfFaces2D*(this->NumberOfLevels-1)));

  // Number of cells in the vertical = number of full levels - 1
  std::vector<vtkIdType> cell_verts;
  cell_verts.resize(2*nverts_per_face);
  for (size_t ilevel = 0; ilevel < this->NumberOfLevels-1; ilevel++)
  {
    for (size_t iface = 0; iface < this->NumberOfFaces2D; iface++)
    {
      for (size_t ivertex = 0; ivertex < nverts_per_face; ivertex++)
      {
        cell_verts[ivertex] = static_cast<vtkIdType>(face_nodes[iface*nverts_per_face + ivertex] + ilevel*nnodes);
        cell_verts[ivertex+nverts_per_face] = static_cast<vtkIdType>(cell_verts[ivertex] + nnodes);
      }
      grid->InsertNextCell(VTK_HEXAHEDRON, static_cast<vtkIdType>(2*nverts_per_face),
                           cell_verts.data());
    }
  }

  if (!this->UseCartCoords)
  {
    mirror_points(grid);
  }

  vtkDebugMacro("Finished CreateVTKGrid" << endl);

  return 1;
}

// Read field data from netCDF file and add to the VTK grid
int vtkNetCDFLFRicReader::LoadFields(const int ncid, vtkUnstructuredGrid *grid, const size_t timestep)
{
  vtkDebugMacro("Entering LoadFields..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure not available.");
    return 0;
  }

  // Need to identify time dimension
  int time_counter_dimid;
  CALL_NETCDF(nc_inq_dimid(ncid, "time_counter", &time_counter_dimid));

  // Get edge-face connectivity for handling W2 horizontal fields
  // We have to assume here that the edges in the half-level edge and
  // half-level face grids coincide
  const std::vector<unsigned long long> edge_faces = getNCVarULongLong(ncid,
                                        "Mesh2d_half_levels_edge_face_links",
                                        {0,0}, {this->NumberOfEdges2D,2});

  SetProgressText("Loading Field Data");

  for (std::map<std::string,bool>::iterator it = this->Fields.begin();
       it != this->Fields.end(); it++)
  {
    // Load variable?
    if (it->second)
    {
      std::string name = it->first;

      vtkDebugMacro("Reading variable " << name << endl);

      // Get variable ID
      int ivar;
      CALL_NETCDF(nc_inq_varid(ncid, name.c_str(), &ivar));

      // Gather variable info
      int ndims;
      int dimids[NC_MAX_VAR_DIMS];
      CALL_NETCDF(nc_inq_var(ncid, ivar, NULL, NULL, &ndims, dimids, NULL));

      // Check if field is time-dependent and store the index of time dimension
      size_t time_dim = -1;
      for (size_t idim = 0; idim < static_cast<size_t>(ndims); idim++ )
      {
        if ( dimids[idim] == time_counter_dimid )
        {
          time_dim = idim;
          vtkDebugMacro("Found time dimension in time_dim=" << time_dim << endl);
          break;
        }
      }

      // Work out field dimensions to create a read buffer of sufficient size
      std::vector<size_t> start;
      start.resize(ndims);
      std::vector<size_t> count;
      count.resize(ndims);
      size_t total_size = 1;
      for (size_t idim = 0; idim < static_cast<size_t>(ndims); idim++)
      {
        if ( idim == time_dim )
        {
          start[idim] = timestep;
          count[idim] = 1;
        }
        else
        {
          size_t dimlen;
          CALL_NETCDF(nc_inq_dimlen(ncid, dimids[idim], &dimlen));
          start[idim] = 0;
          count[idim] = dimlen;
          total_size *= dimlen;
        }
      }
      vtkDebugMacro("Computed total size of field total_size=" << total_size << endl);

      std::vector<double> read_buffer;
      read_buffer.resize(total_size);

      CALL_NETCDF(nc_get_vara_double(ncid, ivar, start.data(), count.data(), read_buffer.data()));

      mesh_types mesh_type;

      // Find out which mesh type is used for this variable
      size_t mesh_name_len;
      CALL_NETCDF(nc_inq_attlen(ncid, ivar, "mesh", &mesh_name_len));
      char *mesh_name = new char[mesh_name_len+1];
      CALL_NETCDF(nc_get_att_text(ncid, ivar, "mesh", mesh_name));

      if ( strncmp(mesh_name, "Mesh2d_half_levels", 18) == 0 )
      {
        mesh_type = half_level_face;
        vtkDebugMacro("Field is defined on half level face mesh" << endl);
        if (total_size != (this->NumberOfLevels-1)*this->NumberOfFaces2D)
        {
          vtkErrorMacro("Unexpected field size.");
        }
      }
      else if ( strncmp(mesh_name, "Mesh2d_edge_half_levels", 23) == 0 )
      {
        mesh_type = half_level_edge;
        vtkDebugMacro("Field is defined on half level edge mesh" << endl);
        if (total_size != (this->NumberOfLevels-1)*this->NumberOfEdges2D)
        {
          vtkErrorMacro("Unexpected field size.");
        }
      }
      else if ( strncmp(mesh_name, "Mesh2d_full_levels", 18) == 0 )
      {
        mesh_type = full_level_face;
        vtkDebugMacro("Field is defined on full level face mesh" << endl);
        if (total_size != this->NumberOfLevels*this->NumberOfFaces2D)
        {
          vtkErrorMacro("Unexpected field size.");
        }
      }
      // Support for nodal grid (or other grids) will be added as needed
      else
      {
        vtkErrorMacro("Unknown mesh.");
        return 0;
      }

      delete [] mesh_name;

      vtkDebugMacro("Setting vtkDoubleArray for this field..." << endl);

      // Create vtkDoubleArray for field data, vector components are stored separately
      vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
      field->SetNumberOfComponents(1);
      field->SetNumberOfTuples((this->NumberOfLevels-1)*this->NumberOfFaces2D);
      field->SetName(name.c_str());

      switch(mesh_type)
      {
        case half_level_face :
          vtkDebugMacro("half level face mesh: no interpolation needed" << endl);
          for (size_t i = 0; i < (this->NumberOfLevels-1)*this->NumberOfFaces2D; i++)
          {
	    field->SetComponent(static_cast<vtkIdType>(i), 0, read_buffer[i]);
          }
          break;
        case full_level_face :
          vtkDebugMacro("full level face mesh: averaging top and bottom faces" << endl);
          for (size_t i = 0; i < (this->NumberOfLevels-1)*this->NumberOfFaces2D; i++)
          {
            const double bottomval = read_buffer[i];
            const double topval = read_buffer[i+this->NumberOfFaces2D];
            field->SetComponent(static_cast<vtkIdType>(i), 0, 0.5*(bottomval+topval));
          }
          break;
        case half_level_edge:
          vtkDebugMacro("half level edge mesh: averaging four edges" << endl);
          field->Fill(0.0);
          for (size_t edge = 0; edge < this->NumberOfEdges2D; edge++)
	  {
            // Each edge is shared by 2 faces
            for (size_t side = 0; side < 2; side++)
	    {
              // Look up face ID, then evaluate entire vertical column
              const size_t face = static_cast<size_t>(edge_faces[edge*2+side]);
              for (size_t level = 0; level < (this->NumberOfLevels-1); level++)
	      {
                const vtkIdType cellId = static_cast<vtkIdType>(face+level*this->NumberOfFaces2D);
                double fieldval = field->GetComponent(cellId, 0);
                fieldval += 0.25*read_buffer[edge+level*this->NumberOfEdges2D];
                field->SetComponent(cellId, 0, fieldval);
              }
            }
          }
          break;
        case nodal:
          vtkErrorMacro("nodal mesh: not currently supported" << endl);
          field->Fill(0.0);
      }

      grid->GetCellData()->AddArray(field);

      this->UpdateProgress(static_cast<float>(ivar)/
                           static_cast<float>(this->Fields.size()));
    }
  }

  vtkDebugMacro("Finished LoadFields" << endl);

  return 1;

}

void vtkNetCDFLFRicReader::SetUseCartCoords(const int status)
{
  this->UseCartCoords = status;
  // Notify pipeline
  this->Modified();
}

int vtkNetCDFLFRicReader::GetNumberOfCellArrays()
{
  return this->Fields.size();
}

const char* vtkNetCDFLFRicReader::GetCellArrayName(const int index)
{
  // Using a map to store fields is convenient for most purposes, but
  // ParaView also wants us to provide array names by index. Input files
  // shouldn't contain too many arrays, so we'll simply count map keys.
  int counter = 0;
  for (std::map<std::string,bool>::iterator it = this->Fields.begin();
                                            it != this->Fields.end();
                                            it++, counter++)
  {
    if (counter == index)
    {
      return it->first.c_str();
    }
  }
  // Return nullptr by default - happens if index is out of range
  vtkDebugMacro("GetCellArrayName: out of range index=" << index << endl);
  return nullptr;
}

int vtkNetCDFLFRicReader::GetCellArrayStatus(const char* name)
{
  int status;
  std::map<std::string,bool>::iterator it = this->Fields.find(name);
  if (it != this->Fields.end())
  {
    status = static_cast<int>(it->second);
  }
  else
  {
    status = 0;
    vtkDebugMacro("GetCellArrayStatus: no array with name=" << name << endl);
  }
  return status;
}

void vtkNetCDFLFRicReader::SetCellArrayStatus(const char* name, const int status)
{
  std::map<std::string,bool>::iterator it = this->Fields.find(name);
  if (it != this->Fields.end())
  {
    it->second = static_cast<bool>(status);
    this->Modified();
  }
  else
  {
    // Ignore unknown names
    vtkDebugMacro("SetCellArrayStatus: no array with name=" << name << endl);
  }
}
