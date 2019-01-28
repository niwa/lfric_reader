#include "vtkNetCDFLFRicReader.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>

#include <vtk_netcdf.h>
#include <math.h>

#define CALL_NETCDF(call) \
  { \
    int errorcode = call; \
    if (errorcode != NC_NOERR) \
    { \
      vtkErrorMacro(<< "netCDF Error: " << nc_strerror(errorcode)); \
      return 0; \
    } \
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkNetCDFLFRicReader);

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
  this->TimeSteps = nullptr;
  this->NumberOfTimeSteps = 0;
  this->NumberOfLevels = 0;
  this->NumberOfFaces2D = 0;

  vtkDebugMacro("Finished vtkNetCDFLFRicReader constructor" << endl);
}

//----------------------------------------------------------------------------
vtkNetCDFLFRicReader::~vtkNetCDFLFRicReader()
{
  vtkDebugMacro("Entering vtkNetCDFLFRicReader destructor..." << endl);

  this->SetFileName(nullptr);
  delete []this->TimeSteps;
  this->TimeSteps = nullptr;

  vtkDebugMacro("Finished vtkNetCDFLFRicReader destructor" << endl);
}

//----------------------------------------------------------------------------
void vtkNetCDFLFRicReader::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(nullptr)") << endl;

  os << indent << "NumberOfTimeSteps: "
     << this->NumberOfTimeSteps << endl;

}

//----------------------------------------------------------------------------

// Retrieve time steps from netCDF file:
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

  int mode = NC_NOWRITE;
  int ncid;
  CALL_NETCDF(nc_open(this->FileName, mode, &ncid));

  // Get number of time steps
  int dimid;
  CALL_NETCDF(nc_inq_dimid(ncid, "time_counter", &dimid));
  size_t size;
  CALL_NETCDF(nc_inq_dimlen(ncid, dimid, &size));
  this->NumberOfTimeSteps = size;
  vtkDebugMacro("Number of time steps in file=" << this->NumberOfTimeSteps << endl);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Read time steps array, if available
  if (this->NumberOfTimeSteps > 0)
  {
    if(this->TimeSteps != nullptr)
    {
      delete []this->TimeSteps;
    }
    this->TimeSteps = new double[this->NumberOfTimeSteps];

    int varid;
    CALL_NETCDF(nc_inq_varid(ncid, "time_instant", &varid));
    size_t start[] = {0};
    size_t count[] = {this->NumberOfTimeSteps};
    CALL_NETCDF(nc_get_vara_double(ncid, varid, start, count, this->TimeSteps));

    // Tell the pipeline which steps are available and their range
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
		 this->TimeSteps, static_cast<int>(this->NumberOfTimeSteps));
    double timeRange[2] = {this->TimeSteps[0],
			   this->TimeSteps[this->NumberOfTimeSteps - 1]};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
    vtkDebugMacro("timeRange=" << timeRange[0] << " " << timeRange[1] << endl);
  }
  else
  {
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
    vtkDebugMacro("Only single time step available" << endl);
  }

  CALL_NETCDF(nc_close(ncid));

  vtkDebugMacro("Finished RequestInformation" << endl);
  
  return 1;

}

//----------------------------------------------------------------------------
// Create VTK grid and load data for a specific time step (defaults to first
// time step if no request has been made)
// Convention is that this function returns 1 on success, or 0 otherwise
int vtkNetCDFLFRicReader::RequestData(vtkInformation *vtkNotUsed(request),
				      vtkInformationVector **vtkNotUsed(inputVector),
				      vtkInformationVector *outputVector)
{
  vtkDebugMacro("Entering RequestData..." << endl);

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkUnstructuredGrid *outputGrid = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Default to first time step if no specific request has been made
  double time = 0.0;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
  {
    time = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  }

  vtkDebugMacro("Requested time=" << time << endl);

  if(this->FileName == nullptr)
  {
    vtkErrorMacro("FileName not set.");
    return 0;
  }
  vtkDebugMacro("FileName=" << FileName << endl);

  int mode = NC_NOWRITE;
  int ncid;
  CALL_NETCDF(nc_open(this->FileName, mode, &ncid));

  // Read UGRID description from file and create unstructured grid
  if (!this->CreateVTKGrid(ncid, outputGrid))
  {
    vtkErrorMacro("Could not create VTK grid.");
  }

  // Load requested field data for requested time step
  // FIXME: add time step here
  if (!this->LoadFields(ncid, outputGrid))
  {
    vtkErrorMacro("Could not load field data.");
  }

  CALL_NETCDF(nc_close(ncid));

  vtkDebugMacro("Finished RequestData" << endl);

  return 1;

}

size_t vtkNetCDFLFRicReader::getNCDim(const int ncid, const char * dimname) {
  vtkDebugMacro("getNCDim: dimname=" << dimname << endl);
  int dimid;
  size_t dim;
  CALL_NETCDF(nc_inq_dimid(ncid, dimname, &dimid));
  CALL_NETCDF(nc_inq_dimlen(ncid, dimid, &dim));
  // FIXME: could return dimid, dim tuple here
  return dim;
}

int vtkNetCDFLFRicReader::getNCVar(const int ncid, const char * varname, const nc_type expect_vartype,
                                   const int expect_ndims, const size_t start [], const size_t count [], void * buffer)
{
  vtkDebugMacro("getNCVar: varname=" << varname << endl);
  int varid;
  CALL_NETCDF(nc_inq_varid(ncid, varname, &varid));
  int actual_ndims, dimids[NC_MAX_VAR_DIMS];
  nc_type actual_vartype;
  CALL_NETCDF(nc_inq_var(ncid, varid, NULL, &actual_vartype, &actual_ndims, dimids, NULL));
  // FIXME: check dimids here
  if( (actual_vartype != expect_vartype) || (actual_ndims != expect_ndims) ) // || (dimids[0] != nnodes_dimid) )
  {
    vtkErrorMacro("Variable " << varname << ": unexpected type, shape, or dimension." << endl);
    return 0;
  }
  CALL_NETCDF(nc_get_vara(ncid, varid, start, count, buffer));
  return 1;
}

// Read UGRID description from netCDF file and build VTK grid
int vtkNetCDFLFRicReader::CreateVTKGrid(const int ncid, vtkUnstructuredGrid *grid)
{
  vtkDebugMacro("Entering CreateVTKGrid..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure is not available.");
    return 0;
  }

  // Get number of nodes
  //  int nnodes_dimid;
  size_t nnodes = getNCDim(ncid, "nMesh2d_full_levels_node");
  vtkDebugMacro("nnodes=" << nnodes << endl);

  // Get number of levels
  //  int nlevels_dimid;
  size_t nlevels = getNCDim(ncid, "full_levels_faces");
  this->NumberOfLevels = nlevels;
  vtkDebugMacro("nlevels=" << nlevels << endl);

  // Get number of faces
  //  int nfaces_dimid;
  size_t nfaces = getNCDim(ncid, "nMesh2d_full_levels_face");
  this->NumberOfFaces2D = nfaces;
  vtkDebugMacro("nfaces=" << nfaces << endl);

  // Get number of vertices per face
  //  int nverts_per_face_dimid;
  size_t nverts_per_face = getNCDim(ncid, "nMesh2d_full_levels_vertex");
  vtkDebugMacro("nverts_per_face=" << nverts_per_face << endl);

  // Create memory buffers for grid description
  float *node_coords_x = new float[nnodes];
  float *node_coords_y = new float[nnodes];
  float *levels = new float[nlevels];
  int *face_nodes = new int[nfaces*nverts_per_face];

  // Get x coordinates
  {
    size_t start[] = {0};
    size_t count[] = {nnodes};
    getNCVar(ncid, "Mesh2d_full_levels_node_x", NC_FLOAT, 1, start, count, node_coords_x);
  }

  // Get y coordinates
  {
    size_t start[] = {0};
    size_t count[] = {nnodes};
    getNCVar(ncid, "Mesh2d_full_levels_node_y", NC_FLOAT, 1, start, count, node_coords_y);
  }

  // Get vertical levels
  {
    size_t start[] = {0};
    size_t count[] = {nlevels};
    getNCVar(ncid, "full_levels_faces", NC_FLOAT, 1, start, count, levels);
  }

  // Node connectivity
  {
    size_t start[] = {0,0};
    size_t count[] = {nfaces,nverts_per_face};
    getNCVar(ncid, "Mesh2d_full_levels_face_nodes", NC_INT, 2, start, count, face_nodes);
  }

  vtkDebugMacro("Setting VTK points..." << endl);

  // Set grid points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (size_t ilevel = 0; ilevel < nlevels; ilevel++)
  {
    for (size_t inode = 0; inode < nnodes; inode++)
    {
      // Avoid zero level height
      const float level_shift = 0.5;

      const float coords[3] = {node_coords_x[inode],
                               node_coords_y[inode],
			       levels[ilevel] + level_shift};

      // Convert from lon-lat-rad to xyz coordinates
      const float xyz[3] = {coords[2]*cosf(coords[0]/180.0*M_PI)*cosf(coords[1]/180.0*M_PI),
			    coords[2]*sinf(coords[0]/180.0*M_PI)*cosf(coords[1]/180.0*M_PI),
                            coords[2]*sinf(coords[1]/180.0*M_PI)};

      points->InsertNextPoint(xyz);
    }
  }
  grid->SetPoints(points);

  vtkDebugMacro("Setting VTK cells..." << endl);

  // Build up grid cells vertical layer-wise
  grid->Allocate(static_cast<vtkIdType>(nfaces*(nlevels-1)));
  // Number of cells in the vertical = number of half levels in file
  for (size_t ilevel = 0; ilevel < nlevels-1; ilevel++)
  {
    for (size_t iface = 0; iface < nfaces; iface++)
    {
      vtkIdType cell_verts[2*nverts_per_face];
      for (size_t ivertex = 0; ivertex < nverts_per_face; ivertex++)
      {
        cell_verts[ivertex] = face_nodes[iface*nverts_per_face + ivertex] + ilevel*nnodes;
	cell_verts[ivertex+nverts_per_face] = cell_verts[ivertex] + nnodes;
      }
      grid->InsertNextCell(VTK_HEXAHEDRON, static_cast<vtkIdType>(2*nverts_per_face), cell_verts);
    }
  }

  delete []node_coords_x;
  delete []node_coords_y;
  delete []levels;
  delete []face_nodes;

  vtkDebugMacro("Finished CreateVTKGrid" << endl);

  return 1;

}

// Read field data from netCDF file and add to the VTK grid
int vtkNetCDFLFRicReader::LoadFields(const int ncid, vtkUnstructuredGrid *grid)
{
  vtkDebugMacro("Entering LoadFields..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure not available.");
    return 0;
  }

  // Work out how many variables the netCDF file contains
  int numVars;
  CALL_NETCDF(nc_inq_nvars(ncid, &numVars));

  vtkDebugMacro("numVars=" << numVars << endl);

  int time_counter_dimid;
  CALL_NETCDF(nc_inq_dimid(ncid, "time_counter", &time_counter_dimid));

  if (numVars > 0)
  {

    for (int ivar = 0; ivar < numVars; ivar++)
    {
      char name[NC_MAX_NAME+1];
      CALL_NETCDF(nc_inq_varname(ncid, ivar, name));

      // Don't read UGRID mesh description and time axis
      if ( (strstr(name, "Mesh2d_") == nullptr) &&
	   (strstr(name, "half_levels_") == nullptr) &&
	   (strstr(name, "full_levels_") == nullptr) &&
	   (strstr(name, "time_instant") == nullptr) )
      {

        vtkDebugMacro("Reading variable " << name << endl);

        // Gather variable info
	int ndims;
        int dimids[NC_MAX_VAR_DIMS];
        nc_type vartype;
        CALL_NETCDF(nc_inq_var(ncid, ivar, NULL, &vartype, &ndims, dimids, NULL));
        if( (vartype != NC_DOUBLE) || (ndims > 3) )
        {
          vtkErrorMacro("Expected all field variables to use double precision and up to 3 dimensions.");
          return 0;
        }

	// Check if field is time-dependent and store the index of time dimension
	size_t time_dim = -1;
	for (size_t idim = 0; idim < ndims; idim++ )
	{
	  if ( dimids[idim] == time_counter_dimid )
	  {
	    time_dim = idim;
            vtkDebugMacro("Found time dimension in time_dim=" << time_dim << endl);
	    break;
	  }
	}

	// Work out field dimensions to create a read buffer of sufficient size
        size_t start[ndims];
        size_t count[ndims];
        size_t total_size = 1;
	for (size_t idim = 0; idim < ndims; idim++)
	{
	  size_t dimlen;
	  CALL_NETCDF(nc_inq_dimlen(ncid, dimids[idim], &dimlen));
	  if ( idim == time_dim )
	  {
	    // FIXME:
	    // Work out correct time index, something along the lines of:
	    // if (timeValues->GetValue(static_cast<vtkIdType>(start[0])) >= time) break;
	    // Or maybe it is possible to get a time index out of VTK's Information class?
	    start[idim] = 0;
	    count[idim] = 1;
	  }
	  else
	  {
            start[idim] = 0;
            count[idim] = dimlen;
	    total_size *= dimlen;
	  }
	}
        vtkDebugMacro("Computed total size of field total_size=" << total_size << endl);

	double *read_buffer = new double[total_size];

        CALL_NETCDF(nc_get_vara_double(ncid, ivar, start, count, read_buffer));

	// Check if variable is defined on full levels or half levels
	bool half_level_mesh;
	size_t mesh_name_len;
	CALL_NETCDF(nc_inq_attlen(ncid, ivar, "mesh", &mesh_name_len));
	char mesh_name[mesh_name_len];
	CALL_NETCDF(nc_get_att_text(ncid, ivar, "mesh", mesh_name));
	if ( strcmp(mesh_name, "Mesh2d_half_levels") == 0 )
	{
          half_level_mesh = true;
          vtkDebugMacro("Field is defined on half level mesh" << endl);
          if (total_size != (this->NumberOfLevels-1)*this->NumberOfFaces2D)
          {
            vtkErrorMacro("Unexpected field size.");
          }
	}
	else if ( strcmp(mesh_name, "Mesh2d_full_levels") == 0 )
	{
          half_level_mesh = false;
          vtkDebugMacro("Field is defined on full level mesh" << endl);
          if (total_size != this->NumberOfLevels*this->NumberOfFaces2D)
          {
            vtkErrorMacro("Unexpected field size.");
          }
	}
	else
	{
	  vtkErrorMacro("Unknown mesh.");
          return 0;
	}

        vtkDebugMacro("Setting vtkDoubleArray for this field..." << endl);

	// Create vtkDoubleArray for field data, vector components are stored separately
        vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
        field->SetNumberOfComponents(1);
        field->SetNumberOfTuples((this->NumberOfLevels-1)*this->NumberOfFaces2D);
        field->SetName(name);

	// Average each cell over top and bottom face in case of full level mesh
        if (half_level_mesh)
	{
          vtkDebugMacro("half level mesh: no interpolation needed" << endl);
          for (vtkIdType i = 0; i < (this->NumberOfLevels-1)*this->NumberOfFaces2D; i++) {
	    field->SetComponent(i, 0, read_buffer[i]);
          }
        }
        else
	{
          vtkDebugMacro("full level mesh: averaging top and bottom faces" << endl);
          for (vtkIdType i = 0; i < (this->NumberOfLevels-1)*this->NumberOfFaces2D; i++) {
            double bottomval = read_buffer[i];
            double topval = read_buffer[i+this->NumberOfFaces2D];
	    field->SetComponent(i, 0, 0.5*(bottomval+topval));
          }
        }

	grid->GetCellData()->AddArray(field);

	delete []read_buffer;
      }
    }
  }

  vtkDebugMacro("Finished LoadFields" << endl);

  return 1;

}
