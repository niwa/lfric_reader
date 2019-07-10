#include "vtkNetCDFLFRicReader.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkUnsignedCharArray.h>

#include <unordered_map>
#include <algorithm>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkNetCDFLFRicReader)

//----------------------------------------------------------------------------
vtkNetCDFLFRicReader::vtkNetCDFLFRicReader()
{

  // Build with -DCMAKE_BUILD_TYPE=Debug to activate this
#ifdef DEBUG
  this->DebugOn();
#endif

  vtkDebugMacro("Entering vtkNetCDFLFRicReader constructor..." << endl);

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->FileName = nullptr;
  this->UseCartCoords = 0;
  this->VerticalScale = 1.0;
  this->VerticalBias = 1.0;
  this->Fields.clear();
  this->TimeSteps.clear();
  this->NumberOfLevelsGlobal = 0;
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
int vtkNetCDFLFRicReader::CanReadFile(const char* fileName)
{
  vtkDebugMacro("Entering CanReadFile..." << endl);

  if(fileName == nullptr)
  {
    vtkErrorMacro("CanReadFile: FileName not set." << endl);
    return 0;
  }
  vtkDebugMacro("fileName=" << fileName << endl);

  netCDFLFRicFile inputFile(fileName);
  if (not inputFile.IsFileOpen())
  {
    vtkErrorMacro("CanReadFile: Failed to open file " << fileName);
    return 0;
  }

  // Full-level mesh is needed to construct VTK grid
  bool valid = false;
  for (std::string const &name : inputFile.GetVarNames())
  {
    valid |= (name.compare("Mesh2d_full_levels") == 0);
  }

  if (valid)
  {
    vtkDebugMacro("Finished CanReadFile (file is valid)" << endl);
    return 1;
  }
  else
  {
    vtkDebugMacro("Finished CanReadFile (file is NOT valid)" << endl);
    return 0;
  }
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
    vtkErrorMacro("FileName not set." << endl);
    return 0;
  }
  vtkDebugMacro("FileName=" << this->FileName << endl);

  // Open input file
  netCDFLFRicFile inputFile(this->FileName);
  if (not inputFile.IsFileOpen())
  {
    vtkErrorMacro("Failed to open file " << this->FileName << endl);
    return 0;
  }

  // Get variable names and populate "Fields" map, ignoring UGRID mesh definitions
  // Also try and find variable with time step times
  std::string timeVarName;
  timeVarName.clear();
  for (std::string const &varName : inputFile.GetVarNames())
  {
    vtkDebugMacro("Considering variable name=" << varName << endl);

    if (inputFile.VarHasAtt(varName, "long_name"))
    {
      // Identify time variable
      if (inputFile.GetAttText(varName, "long_name") == "Time axis")
      {
        timeVarName = varName;
      }
      // Fields should have "units" and "mesh" attributes apart from "long_name"
      else if (inputFile.VarHasAtt(varName, "units") and
               inputFile.VarHasAtt(varName, "mesh"))
      {
        // If field is not in list, insert and default to "don't load"
        std::map<std::string,bool>::const_iterator it = this->Fields.find(varName);
        if (it == this->Fields.end())
        {
          this->Fields.insert(it, std::pair<std::string,bool>(varName,false));
        }
      }
    }
  }
  vtkDebugMacro("Number of data arrays found=" << this->Fields.size() << endl);
  vtkDebugMacro("Name of time variable (if any)=" << timeVarName << endl);

  // Unlike time variable, "time_counter" dimension is always present
  // but may be zero
  size_t NumberOfTimeSteps = inputFile.GetDimLen("time_counter");
  vtkDebugMacro("Number of time steps in file=" << NumberOfTimeSteps << endl);

  // If there are time steps (possibly just one), we always require a time variable
  if (NumberOfTimeSteps > 0 and timeVarName.empty())
  {
    vtkErrorMacro("RequestInformation: Time steps but no time variable found in file." << endl);
    return 0;
  }

  // Get VTK object for handing over time step information
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Read time variable if applicable, and tell the pipeline
  this->TimeSteps.clear();
  if (NumberOfTimeSteps > 0)
  {
    this->TimeSteps = inputFile.GetVarDouble(timeVarName, {0},
                                             {NumberOfTimeSteps});

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

  // Reader supports pieces (partitioning) for parallel operation
  // Data is partitioned by vertical layers, so we can produce
  // only as many pieces as there are layers
  outInfo->Set(vtkAlgorithm::CAN_HANDLE_PIECE_REQUEST(), 1);
  this->NumberOfLevelsGlobal = inputFile.GetDimLen("full_levels")-1;
  vtkDebugMacro("NumberOfLevelsGlobal (max number of pieces)=" <<
                this->NumberOfLevelsGlobal << endl);

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
    vtkErrorMacro("FileName not set." << endl);
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

  // Find out required partition and ghost levels
  int piece = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  int numGhosts = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  vtkDebugMacro("piece=" << piece << " numPieces=" << numPieces <<
                " numGhosts=" << numGhosts << endl);

  if (numPieces > this->NumberOfLevelsGlobal)
  {
    vtkErrorMacro("Pipeline requested " << numPieces <<
                  " pieces but reader can only provide " <<
                  this->NumberOfLevelsGlobal << endl);
    return 0;
  }

  // Distribute vertical levels evenly across pieces. If there
  // is a non-zero remainder, add one layer to first few pieces
  // Use int here to handle potentially negative results when
  // ghost layers are added.
  int numLevels = this->NumberOfLevelsGlobal/numPieces;
  int remainder = this->NumberOfLevelsGlobal%numPieces;
  int startLevel = numLevels*piece + remainder;
  if (piece < remainder)
  {
    numLevels++;
    startLevel = numLevels*piece;
  }
  vtkDebugMacro("startLevel=" << startLevel << " numLevels=" << numLevels << endl);

  // Add ghost levels but limit them to available level range
  // We need to mark ghost cells in VTK grid, so keep track of
  // ghost levels
  int numGhostsBelow = std::min(numGhosts, startLevel);
  int numGhostsAbove = std::min(numGhosts, static_cast<int>(this->NumberOfLevelsGlobal)-
                                (startLevel+numLevels));
  startLevel -= numGhostsBelow;
  numLevels += numGhostsBelow + numGhostsAbove;
  vtkDebugMacro("numGhostsBelow=" << numGhostsBelow <<
                " numGhostsAbove=" << numGhostsAbove <<
                " startLevel=" << startLevel << " numLevels=" <<
                numLevels << endl);

  // Sanity checks
  size_t stopLevel = static_cast<size_t>(startLevel+numLevels);
  if ((startLevel < 0) || (stopLevel > this->NumberOfLevelsGlobal))
  {
    vtkErrorMacro("Erroneous level range encountered: " << startLevel <<
                  "..." << (stopLevel-1) << endl);
    return 0;
  }
  if ((numGhostsBelow < 0) || (numGhostsAbove < 0))
  {
    vtkErrorMacro("Erroneous ghost levels encountered: " << numGhostsBelow <<
                  " " << numGhostsAbove << endl);
    return 0;
  }

  netCDFLFRicFile inputFile(this->FileName);
  if (not inputFile.IsFileOpen())
  {
    vtkErrorMacro("Failed to open file " << this->FileName << endl);
    return 0;
  }

  // Read UGRID description from file and create unstructured VTK grid
  if (!this->CreateVTKGrid(inputFile, outputGrid,
                           static_cast<size_t>(startLevel),
                           static_cast<size_t>(numLevels),
                           static_cast<size_t>(numGhostsAbove),
                           static_cast<size_t>(numGhostsBelow)))
  {
    vtkErrorMacro("Could not create VTK grid." << endl);
    return 0;
  }

  // Load requested field data for requested time step
  if (!this->LoadFields(inputFile, outputGrid, timestep,
                        static_cast<size_t>(startLevel),
                        static_cast<size_t>(numLevels)))
  {
    vtkErrorMacro("Could not load field data." << endl);
    return 0;
  }

  vtkDebugMacro("Finished RequestData" << endl);

  return 1;
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
int vtkNetCDFLFRicReader::CreateVTKGrid(netCDFLFRicFile& inputFile, vtkUnstructuredGrid *grid,
                                        const size_t startLevel, const size_t numLevels,
                                        const size_t numGhostsAbove, const size_t numGhostsBelow)
{
  vtkDebugMacro("Entering CreateVTKGrid..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure is not available." << endl);
    return 0;
  }

  this->SetProgressText("Creating VTK Grid");
  this->UpdateProgress(0.0);

  // Get names of various coordinate arrays
  std::vector<std::string> nodeCoordVarNames =
    inputFile.GetAttTextSplit("Mesh2d_full_levels", "node_coordinates");

  std::vector<std::string> edgeCoordVarNames =
    inputFile.GetAttTextSplit("Mesh2d_full_levels", "edge_coordinates");

  std::string faceNodeConnVarName = inputFile.GetAttText("Mesh2d_full_levels",
                                                         "face_node_connectivity");

  // Read various dimensions, keep some of them in our object
  std::string nodeDimName = inputFile.GetVarDimName(nodeCoordVarNames[0], 0);
  const size_t nnodes = inputFile.GetDimLen(nodeDimName);
  vtkDebugMacro("nodeDimName=" << nodeDimName << " nnodes=" << nnodes << endl);

  std::string faceDimName = inputFile.GetVarDimName(faceNodeConnVarName, 0);
  this->NumberOfFaces2D = inputFile.GetDimLen(faceDimName);
  vtkDebugMacro("faceDimName=" << faceDimName << " NumberOfFaces2D=" <<
                this->NumberOfFaces2D << endl);

  std::string edgeDimName = inputFile.GetVarDimName(edgeCoordVarNames[0], 0);
  this->NumberOfEdges2D = inputFile.GetDimLen(edgeDimName);
  vtkDebugMacro("edgeDimName=" << edgeDimName << " NumberOfEdges2D=" <<
                this->NumberOfEdges2D << endl);

  std::string vertexDimName = inputFile.GetVarDimName(faceNodeConnVarName, 1);
  const size_t nverts_per_face = inputFile.GetDimLen(vertexDimName);
  vtkDebugMacro("vertexDimName=" << vertexDimName << " nverts_per_face=" <<
                nverts_per_face << endl);

  // x (lon) and y (lat) coordinates of grid nodes
  const std::vector<double> node_coords_x = inputFile.GetVarDouble(nodeCoordVarNames[0],
                                                          {0}, {nnodes});
  const std::vector<double> node_coords_y = inputFile.GetVarDouble(nodeCoordVarNames[1],
                                                          {0}, {nnodes});

  // Vertical vertex heights for this piece
  const std::vector<double> levels = inputFile.GetVarDouble("full_levels",
                                                            {startLevel}, {numLevels+1});

  // Node connectivity
  const std::vector<long long> face_nodes = inputFile.GetVarLongLong(
                               faceNodeConnVarName,
                               {0,0}, {this->NumberOfFaces2D,nverts_per_face});

  // Work out horizontal grid ranges to distinguish cubed-sphere and biperiodic grids
  // This could be based on netCDF variable attributes if available
  double xmin = *std::min_element(node_coords_x.begin(), node_coords_x.end());
  double xmax = *std::max_element(node_coords_x.begin(), node_coords_x.end());
  double ymin = *std::min_element(node_coords_y.begin(), node_coords_y.end());
  double ymax = *std::max_element(node_coords_y.begin(), node_coords_y.end());

  double xshift;
  double yshift;

  // We expect a longitude range of 0..360 degrees for a cubed-sphere grid
  if (xmin >= 0.0 and xmax <= 360.0)
  {
    xshift = 180.0;
    yshift = 0.0;
    vtkDebugMacro("Detected cubed-sphere grid, setting xshift=" <<
                  xshift << " yshift=" << yshift << endl);
  }
  else
  {
    xshift = 0.5*(xmin + xmax);
    yshift = 0.5*(ymin + ymax);
    vtkDebugMacro("Detected grid other than cubed-sphere, setting xshift=" <<
                  xshift << " yshift=" << yshift << endl);
  }

  vtkDebugMacro("Setting VTK points..." << endl);

  // Construct VTK grid points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  if (this->UseCartCoords)
  {
    const double deg2rad = vtkMath::Pi()/180.0;
    for (size_t ilevel = 0; ilevel < (numLevels+1); ilevel++)
    {
      for (size_t inode = 0; inode < nnodes; inode++)
      {
        const double coords[3] = {node_coords_x[inode],
                                  node_coords_y[inode],
                                  this->VerticalScale*(levels[ilevel] + this->VerticalBias)};

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
    for (size_t ilevel = 0; ilevel < (numLevels+1); ilevel++)
    {
      for (size_t inode = 0; inode < nnodes; inode++)
      {
        // Shift horizontal coordinates to enable resolution of grid periodicity
        const double coords[3] = {node_coords_x[inode]-xshift,
                                  node_coords_y[inode]-yshift,
                                  this->VerticalScale*(levels[ilevel] + this->VerticalBias)};
        points->InsertNextPoint(coords);
      }
    }
  }
  grid->SetPoints(points);

  this->UpdateProgress(0.5);

  vtkDebugMacro("Setting VTK cells..." << endl);

  // Build up grid cells for this piece vertical layer-wise
  grid->Allocate(static_cast<vtkIdType>(this->NumberOfFaces2D*numLevels));

  std::vector<vtkIdType> cell_verts;
  cell_verts.resize(2*nverts_per_face);
  for (size_t ilevel = 0; ilevel < numLevels; ilevel++)
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

  // Mark ghost cells as "duplicate cell"
  if (numGhostsAbove > 0 || numGhostsBelow > 0)
  {
    grid->AllocateCellGhostArray();
    vtkUnsignedCharArray * ghosts = grid->GetCellGhostArray();
    vtkIdType cellId = 0;
    for (size_t ilevel = 0; ilevel < numLevels; ilevel++)
    {
      if ((ilevel < numGhostsBelow) || (ilevel > (numLevels-numGhostsAbove-1)))
      {
        for (size_t iface = 0; iface < this->NumberOfFaces2D; iface++)
        {
          ghosts->SetValue(cellId, vtkDataSetAttributes::DUPLICATECELL);
          cellId++;
        }
      }
      else
      {
        cellId += this->NumberOfFaces2D;
      }
    }
  }

  // Need to replicate (mirror) points in periodic grid
  if (!this->UseCartCoords)
  {
    mirror_points(grid);
  }

  this->UpdateProgress(1.0);

  vtkDebugMacro("Finished CreateVTKGrid" << endl);

  return 1;
}

// Read field data from netCDF file and add to the VTK grid
int vtkNetCDFLFRicReader::LoadFields(netCDFLFRicFile& inputFile, vtkUnstructuredGrid* grid,
                                     const size_t timestep, const size_t startLevel,
                                     const size_t numLevels)
{
  vtkDebugMacro("Entering LoadFields..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure not available." << endl);
    return 0;
  }

  // This code is currently disabled until edge-centered field can be handled correctly
  // Get edge-face connectivity for handling W2 horizontal fields
  // We have to assume here that the edges in the half-level edge and
  // half-level face grids coincide
  // const std::vector<long long> hl_edge_faces = inputFile.GetVarLongLong(
  //                              "Mesh2d_half_levels_edge_face_links",
  //                              {0,0}, {this->NumberOfEdges2D,2});

  this->SetProgressText("Loading Field Data");

  int fieldCount = 0;
  for (std::pair<std::string,bool> const &field : this->Fields)
  {
    // Load variable?
    if (field.second)
    {
      std::string varName = field.first;

      vtkDebugMacro("Reading variable " << varName << endl);

      // Find out if there is a time dimension
      bool hasTimeDim = inputFile.VarHasDim(varName, "time_counter");
      if (hasTimeDim)
      {
        vtkDebugMacro("Found time dimension" << endl);
      }
      else
      {
        vtkDebugMacro("Found NO time dimension" << endl);
      }

      std::vector<double> read_buffer;

      // Find out which mesh type is used for this variable, we already
      // know that this variable has the "mesh" attribute
      mesh_types mesh_type;
      std::string mesh_name = inputFile.GetAttText(varName, "mesh");

      if ( mesh_name == "Mesh2d_half_levels" )
      {
        mesh_type = half_level_face;
        vtkDebugMacro("Field is defined on half level face mesh" << endl);
        if (hasTimeDim)
        {
          read_buffer = inputFile.GetVarDouble(varName, {timestep,startLevel,0}, {1,numLevels,this->NumberOfFaces2D});
        }
        else
        {
          read_buffer = inputFile.GetVarDouble(varName, {startLevel,0}, {numLevels,this->NumberOfFaces2D});
        }
      }
      else if ( mesh_name == "Mesh2d_edge_half_levels" )
      {
        mesh_type = half_level_edge;
        vtkDebugMacro("Field is defined on half level edge mesh" << endl);
        if (hasTimeDim)
        {
          read_buffer = inputFile.GetVarDouble(varName, {timestep,startLevel,0}, {1,numLevels,this->NumberOfEdges2D});
        }
        else
        {
          read_buffer = inputFile.GetVarDouble(varName, {startLevel,0}, {numLevels,this->NumberOfEdges2D});
        }
      }
      else if ( mesh_name == "Mesh2d_full_levels" )
      {
        mesh_type = full_level_face;
        vtkDebugMacro("Field is defined on full level face mesh" << endl);
        if (hasTimeDim)
        {
          read_buffer = inputFile.GetVarDouble(varName, {timestep,startLevel,0}, {1,numLevels+1,this->NumberOfFaces2D});
        }
        else
        {
          read_buffer = inputFile.GetVarDouble(varName, {startLevel,0}, {numLevels+1,this->NumberOfFaces2D});
        }
      }
      // Support for nodal grid (or other grids) will be added as needed
      else
      {
        vtkErrorMacro("Unknown mesh." << endl);
        return 0;
      }

      vtkDebugMacro("Setting vtkDoubleArray for this field..." << endl);

      // Create vtkDoubleArray for field data, vector components are stored separately
      vtkSmartPointer<vtkDoubleArray> dataField = vtkSmartPointer<vtkDoubleArray>::New();
      dataField->SetNumberOfComponents(1);
      dataField->SetNumberOfTuples(numLevels*this->NumberOfFaces2D);
      dataField->SetName(varName.c_str());

      switch(mesh_type)
      {
        case half_level_face :
          vtkDebugMacro("half level face mesh: no interpolation needed" << endl);
          for (size_t i = 0; i < numLevels*this->NumberOfFaces2D; i++)
          {
            dataField->SetComponent(static_cast<vtkIdType>(i), 0, read_buffer[i]);
          }
          break;
        case full_level_face :
          vtkDebugMacro("full level face mesh: averaging top and bottom faces" << endl);
          for (size_t i = 0; i < numLevels*this->NumberOfFaces2D; i++)
          {
            const double bottomval = read_buffer[i];
            const double topval = read_buffer[i+this->NumberOfFaces2D];
            dataField->SetComponent(static_cast<vtkIdType>(i), 0, 0.5*(bottomval+topval));
          }
          break;
        case half_level_edge:
          vtkWarningMacro("WARNING: edge-centered fields cannot be handled correctly at the moment" << endl);
          vtkDebugMacro("half level edge mesh: averaging four edges" << endl);
          dataField->Fill(0.0);

          // This code is currently disabled until edge-centered field can be handled correctly
          // It also needs to be modified to handle piece requests (startLevel/stopLevel)
          // for (size_t edge = 0; edge < this->NumberOfEdges2D; edge++)
          // {
          //   // Each edge is shared by 2 faces
          //   for (size_t side = 0; side < 2; side++)
          //   {
          //     // Look up face ID, then evaluate entire vertical column
          //     const size_t face = static_cast<size_t>(hl_edge_faces[edge*2+side]);
          //     for (size_t level = 0; level < (this->NumberOfLevels-1); level++)
          //     {
          //       const vtkIdType cellId = static_cast<vtkIdType>(face+level*this->NumberOfFaces2D);
          //       double fieldval = dataField->GetComponent(cellId, 0);
          //       // Pick up data in edge grid
          //       fieldval += 0.25*read_buffer[edge+level*this->NumberOfEdges2D];
          //       dataField->SetComponent(cellId, 0, fieldval);
          //     }
          //   }
          // }

          break;
        case nodal:
          vtkErrorMacro("nodal mesh: not currently supported" << endl);
          dataField->Fill(0.0);
      }

      grid->GetCellData()->AddArray(dataField);

      fieldCount++;
      this->UpdateProgress(static_cast<float>(fieldCount)/
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

void vtkNetCDFLFRicReader::SetVerticalScale(const double value)
{
  this->VerticalScale = value;
  this->Modified();
}

void vtkNetCDFLFRicReader::SetVerticalBias(const double value)
{
  this->VerticalBias = value;
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
  std::map<std::string,bool>::const_iterator it;
  for (it = this->Fields.begin(); it != this->Fields.end();
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
  std::map<std::string,bool>::const_iterator it = this->Fields.find(name);
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
