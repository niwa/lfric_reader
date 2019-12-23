#include "vtkNetCDFLFRicReader.h"
#include "netCDFLFRicReaderUtils.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDataSetAttributes.h>

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
  this->UseIndexAsVertCoord = 0;
  this->VerticalScale = 1.0;
  this->VerticalBias = 1.0;
  this->CellFields.clear();
  this->PointFields.clear();
  this->TimeSteps.clear();
  this->OutputMode = 0;

  vtkDebugMacro("Finished vtkNetCDFLFRicReader constructor" << endl);
}

//----------------------------------------------------------------------------
vtkNetCDFLFRicReader::~vtkNetCDFLFRicReader()
{
  vtkDebugMacro("Entering vtkNetCDFLFRicReader destructor..." << endl);

  this->SetFileName(nullptr);
  this->CellFields.clear();
  this->PointFields.clear();
  this->TimeSteps.clear();

  vtkDebugMacro("Finished vtkNetCDFLFRicReader destructor" << endl);
}

//----------------------------------------------------------------------------
void vtkNetCDFLFRicReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(nullptr)") << endl;

  os << indent << "OutputMode: "
     << this->GetOutputMode() << endl;

  os << indent << "NumberOfTimeSteps: "
     << this->TimeSteps.size() << endl;

  os << indent << "NumberOfCellFields: "
     << this->CellFields.size() << endl;

  os << indent << "NumberOfPointFields: "
     << this->PointFields.size() << endl;

  os << indent << "NumberOfCellArrays: "
     << this->GetNumberOfCellArrays() << endl;

  os << indent << "NumberOfPointArrays: "
     << this->GetNumberOfPointArrays() << endl;

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

  // Need at least a valid UGRID mesh description
  const UGRIDMeshDescription mesh = inputFile.GetMesh2DDescription();
  if (mesh.numTopologies > 0)
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

  // Look for UGRID horizontal unstructured mesh
  this->mesh2D = inputFile.GetMesh2DDescription();
  if (this->mesh2D.numTopologies == 0)
  {
    vtkErrorMacro("Failed to determine 2D UGRID mesh description." << endl);
    return 0;
  }

  // Look for vertical axis - may be absent
  this->zAxis = inputFile.GetZAxisDescription(this->mesh2D.isLFRicXIOSFile,
                                              this->mesh2D.meshType);
  vtkDebugMacro("Vertical axis:" << zAxis.axisVar << endl);

  // Look for time axis - may be absent
  this->tAxis = inputFile.GetTAxisDescription();
  vtkDebugMacro("Time axis:" << tAxis.axisVar << endl);

  // If this is an LFRic output file, add fields that are defined on all supported meshes,
  // otherwise accept only fields on the single mesh that is expected to be of half-level type
  if (this->mesh2D.isLFRicXIOSFile)
  {
    inputFile.UpdateFieldMap(this->CellFields, "face", "nMesh2d_half_levels_face", halfLevelFaceMesh,
                             "half_levels", tAxis.axisDim);
    inputFile.UpdateFieldMap(this->CellFields, "face", "nMesh2d_full_levels_face", fullLevelFaceMesh,
                             "full_levels", tAxis.axisDim);
    inputFile.UpdateFieldMap(this->PointFields, "edge", "nMesh2d_edge_half_levels_edge", halfLevelEdgeMesh,
                             "half_levels", tAxis.axisDim);
  }
  else
  {
    inputFile.UpdateFieldMap(this->PointFields, "face", this->mesh2D.faceDim, halfLevelFaceMesh,
                             zAxis.axisDim, tAxis.axisDim);
  }
  vtkDebugMacro("Number of cell fields found: " << this->CellFields.size() << endl);
  vtkDebugMacro("Number of point fields found: " << this->PointFields.size() << endl);

  // Get VTK object for handing over time step information
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Read time variable if applicable, and tell the pipeline
  if (this->tAxis.numTimeSteps > 0)
  {
    this->TimeSteps = inputFile.GetVarDouble(this->tAxis.axisVar, {0},
                                             {this->tAxis.numTimeSteps});

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

  if (numPieces > this->zAxis.numLevels)
  {
    vtkErrorMacro("Pipeline requested " << numPieces <<
                  " pieces but reader can only provide " <<
                  this->zAxis.numLevels << endl);
    return 0;
  }

  // Distribute vertical levels evenly across pieces. If there
  // is a non-zero remainder, add one layer to first few pieces
  // Use int here to handle potentially negative results when
  // ghost layers are added.
  int numLevels = this->zAxis.numLevels/numPieces;
  int remainder = this->zAxis.numLevels%numPieces;
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
  int numGhostsAbove = std::min(numGhosts, static_cast<int>(this->zAxis.numLevels)-
                                (startLevel+numLevels));
  startLevel -= numGhostsBelow;
  numLevels += numGhostsBelow + numGhostsAbove;
  vtkDebugMacro("numGhostsBelow=" << numGhostsBelow <<
                " numGhostsAbove=" << numGhostsAbove <<
                " startLevel=" << startLevel << " numLevels=" <<
                numLevels << endl);

  // Sanity checks
  size_t stopLevel = static_cast<size_t>(startLevel+numLevels);
  if ((startLevel < 0) || (stopLevel > this->zAxis.numLevels))
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

  // Full VTK grid
  if (this->OutputMode == 0)
  {
    // Read UGRID description from file and create full unstructured VTK grid
    if (!this->CreateVTKGrid(inputFile, outputGrid,
                             static_cast<size_t>(startLevel),
                             static_cast<size_t>(numLevels),
                             static_cast<size_t>(numGhostsAbove),
                             static_cast<size_t>(numGhostsBelow)))
    {
      vtkErrorMacro("Could not create VTK grid." << endl);
      return 0;
    }

    // Load requested field data for requested time step - currently only cell fields
    if (!this->LoadFields(inputFile, outputGrid, this->CellFields,
                          timestep,
                          static_cast<size_t>(startLevel),
                          static_cast<size_t>(numLevels)))
    {
      vtkErrorMacro("Could not load cell field data." << endl);
      return 0;
    }
  }
  // VTK points only for W2 visualisation
  else if (this->OutputMode == 1)
  {
    if (!this->CreateVTKPoints(inputFile, outputGrid,
                               static_cast<size_t>(startLevel),
                               static_cast<size_t>(numLevels),
                               static_cast<size_t>(numGhostsAbove),
                               static_cast<size_t>(numGhostsBelow)))
    {
      vtkErrorMacro("Could not create VTK points." << endl);
      return 0;
    }

    // Load W2 fields
    if (!this->LoadFields(inputFile, outputGrid, this->PointFields,
                          timestep,
                          static_cast<size_t>(startLevel),
                          static_cast<size_t>(numLevels)))
    {
      vtkErrorMacro("Could not load point field data." << endl);
      return 0;
    }
  }

  vtkDebugMacro("Finished RequestData" << endl);

  return 1;
}

//----------------------------------------------------------------------------
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

  //
  // Load node coordinates and face-node connectivities
  //

  std::vector<double> node_coords_x = inputFile.GetVarDouble(this->mesh2D.nodeCoordXVar,
                                                             {0}, {this->mesh2D.numNodes});

  std::vector<double> node_coords_y = inputFile.GetVarDouble(this->mesh2D.nodeCoordYVar,
                                                             {0}, {this->mesh2D.numNodes});

  std::vector<long long> face_nodes = inputFile.GetVarLongLong(
                         this->mesh2D.faceNodeConnVar,
                         {0,0}, {this->mesh2D.numFaces, this->mesh2D.numVertsPerFace});

  // Correct node IDs if non-zero start index
  if (this->mesh2D.faceNodeStartIdx != 0)
  {
    for (size_t idx = 0; idx < face_nodes.size(); idx++)
    {
      face_nodes[idx] -= this->mesh2D.faceNodeStartIdx;
    }
    vtkDebugMacro("Corrected face-node connectivity for start index=" <<
                   this->mesh2D.faceNodeStartIdx << endl);
  }

  //
  // Determine vertical vertex heights
  //

  std::vector<double> levels;
  // Treat 2D case as a single 3D layer
  if (this->UseIndexAsVertCoord or this->zAxis.numLevels == 1)
  {
    // Use level indices
    levels.resize(numLevels+1);
    for (size_t idx = 0; idx < (numLevels+1); idx++)
    {
      levels[idx] = static_cast<double>(startLevel + idx);
    }
  }
  else if (this->mesh2D.meshType == fullLevelFaceMesh)
  {
    // Load vertex heights from file
    levels = inputFile.GetVarDouble(this->zAxis.axisVar,
                                    {startLevel}, {numLevels+1});
  }
  else if (this->mesh2D.meshType == halfLevelFaceMesh)
  {
    std::vector<double> halfLevels = inputFile.GetVarDouble(this->zAxis.axisVar,
                                     {0}, {this->zAxis.numLevels});

    // Extrapolate half-level heights at both ends
    const double firstLevel = 2.0*halfLevels[0]-halfLevels[1];
    halfLevels.insert(halfLevels.begin(),firstLevel);

    const double lastLevel = 2.0*halfLevels[this->zAxis.numLevels-1] -
                                 halfLevels[this->zAxis.numLevels-2];
    halfLevels.push_back(lastLevel);

    // Vertices are in the middle between half-level heights
    levels.resize(numLevels+1);
    for (size_t idx = 0; idx < (numLevels+1); idx++)
    {
      levels[idx] = 0.5*(halfLevels[startLevel+idx]+halfLevels[startLevel+idx+1]);
    }
  }
  else
  {
    vtkErrorMacro("Cannot determine vertical vertex heights." << endl);
    return 0;
  }

  // Apply scale and bias
  for (size_t iLevel = 0; iLevel < levels.size(); iLevel++)
  {
    levels[iLevel] = this->VerticalScale*(levels[iLevel] + this->VerticalBias);
  }

  this->UpdateProgress(0.25);

  //
  // Resolve grid periodicity by mirroring nodes
  //

  // Keep track of node count, in case that nodes are added to resolve
  // grid periodicity
  size_t numNodesCurrent = this->mesh2D.numNodes;

  if (!this->UseCartCoords)
  {
    vtkDebugMacro("Resolving grid periodicity..." << endl);

    resolvePeriodicGrid(node_coords_x, node_coords_y, face_nodes,
                        this->mesh2D.numFaces, this->mesh2D.numVertsPerFace);
    // Update node count to allow for added nodes
    numNodesCurrent = node_coords_x.size();
  }

  this->UpdateProgress(0.5);

  //
  // Construct VTK grid points
  //

  vtkDebugMacro("Setting VTK points..." << endl);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(numNodesCurrent*(numLevels+1));
  vtkDataArray * pointLocs = points->GetData();

  // Use vtkArrayDispatch mechanism to fill vtkDataArray directly, this is
  // a lot faster than using Insert or Set methods
  SetPointLocationWorker pointsWorker(node_coords_x, node_coords_y, levels,
                                      this->UseCartCoords);
  typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals> pointsDispatcher;
  if (!pointsDispatcher::Execute(pointLocs, pointsWorker))
  {
    pointsWorker(pointLocs);
  }

  grid->SetPoints(points);

  this->UpdateProgress(0.75);

  //
  // Construct VTK cells
  //

  vtkDebugMacro("Setting VTK cells..." << endl);

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  cells->SetNumberOfCells(numLevels*this->mesh2D.numFaces);

  vtkDataArray * cellConnect = cells->GetData();
  cellConnect->SetNumberOfTuples(numLevels*this->mesh2D.numFaces*
                                 (2*this->mesh2D.numVertsPerFace+1));

  SetConnectivityWorker connectWorker(face_nodes, numLevels, this->mesh2D.numFaces,
                                      this->mesh2D.numVertsPerFace, numNodesCurrent);
  typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Integrals> cellsDispatcher;
  if (!cellsDispatcher::Execute(cellConnect, connectWorker))
  {
    connectWorker(cellConnect);
  }
  grid->SetCells(VTK_HEXAHEDRON, cells);

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
        for (size_t iface = 0; iface < this->mesh2D.numFaces; iface++)
        {
          ghosts->SetValue(cellId, vtkDataSetAttributes::DUPLICATECELL);
          cellId++;
        }
      }
      else
      {
        cellId += this->mesh2D.numFaces;
      }
    }
  }

  this->UpdateProgress(1.0);

  vtkDebugMacro("Finished CreateVTKGrid" << endl);

  return 1;
}

//----------------------------------------------------------------------------
int vtkNetCDFLFRicReader::CreateVTKPoints(netCDFLFRicFile& inputFile, vtkUnstructuredGrid *grid,
                                          const size_t startLevel, const size_t numLevels,
                                          const size_t numGhostsAbove, const size_t numGhostsBelow)
{
  vtkDebugMacro("Entering CreateVTKPoints..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure is not available." << endl);
    return 0;
  }

  this->SetProgressText("Creating VTK Points");
  this->UpdateProgress(0.0);

  std::vector<double> edge_coords_x = inputFile.GetVarDouble(this->mesh2D.edgeCoordXVar,
                                                             {0}, {this->mesh2D.numEdges});

  std::vector<double> edge_coords_y = inputFile.GetVarDouble(this->mesh2D.edgeCoordYVar,
                                                             {0}, {this->mesh2D.numEdges});

  //
  // Determine vertical vertex heights
  //

  std::vector<double> levels;
  if (this->UseIndexAsVertCoord or this->zAxis.numLevels == 1)
  {
    levels.resize(numLevels);
    for (size_t idx = 0; idx < numLevels; idx++)
    {
      levels[idx] = static_cast<double>(startLevel + idx);
      levels[idx] = this->VerticalScale*(levels[idx] + this->VerticalBias);
    }
  }
  else if (this->mesh2D.meshType == fullLevelFaceMesh)
  {
    levels = inputFile.GetVarDouble(this->zAxis.axisVar,
                                    {startLevel}, {numLevels+1});

    // Compute half-level heights for points visualisation,
    // last vector element is not needed
    for (size_t idx = 0; idx < numLevels; idx++)
    {
      levels[idx] = 0.5*(levels[idx]+levels[idx+1]);
      levels[idx] = this->VerticalScale*(levels[idx] + this->VerticalBias);
    }
    levels.resize(numLevels);
  }
  else
  {
    vtkErrorMacro("Cannot determine vertical vertex heights." << endl);
    return 0;
  }

  this->UpdateProgress(0.5);

  //
  // Construct VTK grid points
  //

  vtkDebugMacro("Setting VTK points..." << endl);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(this->mesh2D.numEdges*numLevels);
  vtkDataArray * pointLocs = points->GetData();

  SetPointLocationWorker pointsWorker(edge_coords_x, edge_coords_y, levels,
                                      this->UseCartCoords);
  typedef vtkArrayDispatch::DispatchByValueType<vtkArrayDispatch::Reals> pointsDispatcher;
  if (!pointsDispatcher::Execute(pointLocs, pointsWorker))
  {
    pointsWorker(pointLocs);
  }

  grid->SetPoints(points);

  // Mark ghost points as "duplicate point"
  if (numGhostsAbove > 0 || numGhostsBelow > 0)
  {
    grid->AllocatePointGhostArray();
    vtkUnsignedCharArray * ghosts = grid->GetPointGhostArray();
    vtkIdType pointId = 0;
    for (size_t ilevel = 0; ilevel < numLevels; ilevel++)
    {
      if ((ilevel < numGhostsBelow) || (ilevel > (numLevels-numGhostsAbove-1)))
      {
        for (size_t iEdge = 0; iEdge < this->mesh2D.numEdges; iEdge++)
        {
          ghosts->SetValue(pointId, vtkDataSetAttributes::DUPLICATEPOINT);
          pointId++;
        }
      }
      else
      {
        pointId += this->mesh2D.numEdges;
      }
    }
  }

  this->UpdateProgress(1.0);

  vtkDebugMacro("Finished CreateVTKPoints" << endl);

  return 1;
}

//----------------------------------------------------------------------------
// Read field data from netCDF file and add to the VTK grid
int vtkNetCDFLFRicReader::LoadFields(netCDFLFRicFile & inputFile, vtkUnstructuredGrid * grid,
                                     const std::map<std::string, DataField> & fields,
                                     const size_t timestep, const size_t startLevel,
                                     const size_t numLevels)
{
  vtkDebugMacro("Entering LoadFields..." << endl);

  if (grid == nullptr)
  {
    vtkErrorMacro("Grid data structure not available." << endl);
    return 0;
  }

  this->SetProgressText("Loading Field Data");
  this->UpdateProgress(0.0);

  int fieldCount = 0;
  for (std::pair<std::string, DataField> const &field : fields)
  {
    // Load variable?
    if (field.second.active)
    {
      const std::string & fieldName = field.first;
      const DataField & fieldSpec = field.second;

      vtkDebugMacro("Reading variable " << fieldName << endl);

      std::vector<double> read_buffer;
      std::vector<size_t> start;
      std::vector<size_t> count;

      const size_t numFieldDims = fieldSpec.dims.size();
      start.resize(numFieldDims);
      count.resize(numFieldDims);

      // Number of vertical levels is imposed by numLevels
      // Number of timesteps to be read is always 1
      size_t numHorizontal = 1;
      size_t numComponents = 1;

      // Need strides for flexible dimension ordering
      size_t componentStride = 1;
      size_t horizontalStride = 1;
      size_t verticalStride = 1;

      // Retrieve field dimensions and strides
      // Horizontal and component dimensions are always read in full
      for (size_t dimIdx = 0; dimIdx < numFieldDims; dimIdx++)
      {
        switch (fieldSpec.dims[dimIdx].dimType)
	{
          case horizontalAxisDim:
            numHorizontal = fieldSpec.dims[dimIdx].dimLength;
            horizontalStride = fieldSpec.dims[dimIdx].dimStride;
            start[dimIdx] = 0;
            count[dimIdx] = numHorizontal;
            break;

          case componentAxisDim:
            numComponents = fieldSpec.dims[dimIdx].dimLength;
            componentStride = fieldSpec.dims[dimIdx].dimStride;
            start[dimIdx] = 0;
            count[dimIdx] = numComponents;
            break;

          case verticalAxisDim:
            verticalStride = fieldSpec.dims[dimIdx].dimStride;
            start[dimIdx] = startLevel;
            count[dimIdx] = numLevels;
            // Full-level meshes are defined on cell interfaces,
            // need to read extra level in that case
            if (fieldSpec.meshType == fullLevelFaceMesh)
	    {
              count[dimIdx]++;
            }
            break;

          case timeAxisDim:
            start[dimIdx] = timestep;
            count[dimIdx] = 1;
            break;

          default:
            vtkErrorMacro("Unknown field dimension." << endl);
            return 0;
        }
        vtkDebugMacro("Dim " << dimIdx << " of " << numFieldDims << ": load " <<
                      count[dimIdx] << " of " << fieldSpec.dims[dimIdx].dimLength <<
                      " elements at index " << start[dimIdx] << "\n");
      }

      // Read data if field has a vertical dimension, or if this piece
      // contains the surface level
      if (fieldSpec.hasVerticalDim or (startLevel == 0))
      {
        read_buffer = inputFile.GetVarDouble(fieldName, start, count);
      }

      // Fill remaining space with NaNs if vertical dimension not present (2D fields)
      if (not fieldSpec.hasVerticalDim)
      {
        // Vertical dimension is outermost
        verticalStride = numHorizontal*numComponents;

        if (fieldSpec.meshType == halfLevelFaceMesh)
        {
          read_buffer.resize(numLevels*numHorizontal*numComponents,
                             std::numeric_limits<double>::quiet_NaN());
        }
        else if (fieldSpec.meshType == fullLevelFaceMesh)
        {
          read_buffer.resize((numLevels+1)*numHorizontal*numComponents,
                             std::numeric_limits<double>::quiet_NaN());
        }
      }

      // Disable component dimension if not present by setting stride to full
      // (possibly resized) buffer size
      if (not fieldSpec.hasComponentDim) componentStride = read_buffer.size();

      vtkDebugMacro("Setting vtkDoubleArray for this field..." << endl);

      // Create vtkDoubleArray for field data, vector components are stored separately
      vtkSmartPointer<vtkDoubleArray> dataField = vtkSmartPointer<vtkDoubleArray>::New();
      dataField->SetNumberOfComponents(numComponents);
      dataField->SetName(fieldName.c_str());

      // Sort data into VTK data structure
      switch(fieldSpec.meshType)
      {
        case halfLevelEdgeMesh :
          vtkDebugMacro("half level edge mesh: no interpolation needed" << endl);
          dataField->SetNumberOfTuples(numLevels*this->mesh2D.numEdges);
          for (size_t bufferIdx = 0; bufferIdx < read_buffer.size(); bufferIdx++)
          {
            const size_t iHorizontal = bufferIdx/horizontalStride % numHorizontal;
            const size_t iVertical = bufferIdx/verticalStride % numLevels;
            const size_t iComponent = bufferIdx/componentStride % numComponents;
            const vtkIdType iCell = static_cast<vtkIdType>(iVertical*numHorizontal + iHorizontal);
            dataField->SetComponent(iCell, iComponent, read_buffer[bufferIdx]);
          }
          break;
        case halfLevelFaceMesh :
          vtkDebugMacro("half level face mesh: no interpolation needed" << endl);
          dataField->SetNumberOfTuples(numLevels*this->mesh2D.numFaces);
          for (size_t bufferIdx = 0; bufferIdx < read_buffer.size(); bufferIdx++)
          {
            const size_t iHorizontal = bufferIdx/horizontalStride % numHorizontal;
            const size_t iVertical = bufferIdx/verticalStride % numLevels;
            const size_t iComponent = bufferIdx/componentStride % numComponents;
            const vtkIdType iCell = static_cast<vtkIdType>(iVertical*numHorizontal + iHorizontal);
            dataField->SetComponent(iCell, iComponent, read_buffer[bufferIdx]);
          }
          break;
        case fullLevelFaceMesh :
          vtkDebugMacro("full level face mesh: averaging top and bottom faces" << endl);
          dataField->SetNumberOfTuples(numLevels*this->mesh2D.numFaces);
          const size_t numLevelsFull = numLevels+1;
          for (size_t bufferIdx = 0; bufferIdx < read_buffer.size(); bufferIdx++)
          {
            const size_t iHorizontal = bufferIdx/horizontalStride % numHorizontal;
            const size_t iVertical = bufferIdx/verticalStride % numLevelsFull;
            const size_t iComponent = bufferIdx/componentStride % numComponents;
            const vtkIdType iCell = static_cast<vtkIdType>(iVertical*numHorizontal + iHorizontal);
            // Skip the last full-level as we are averaging onto half-levels
            if (iVertical < numLevels)
            {
              const double fieldValue = 0.5*(read_buffer[bufferIdx]+
                                             read_buffer[bufferIdx+verticalStride]);
              dataField->SetComponent(iCell, iComponent, fieldValue);
            }
          }
          break;
      }

      if (fieldSpec.location == cellFieldLoc)
      {
        grid->GetCellData()->AddArray(dataField);
      }
      else if (fieldSpec.location == pointFieldLoc)
      {
        grid->GetPointData()->AddArray(dataField);
      }
      else
      {
        vtkDebugMacro("Unknown location for field " << fieldName << endl);
      }

      fieldCount++;
      this->UpdateProgress(static_cast<float>(fieldCount)/
                           static_cast<float>(fields.size()));
    }
  }

  vtkDebugMacro("Finished LoadFields" << endl);

  return 1;

}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetUseIndexAsVertCoord(const int status)
{
  this->UseIndexAsVertCoord = status;
  // Notify pipeline
  this->Modified();
}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetUseCartCoords(const int status)
{
  this->UseCartCoords = status;
  this->Modified();
}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetVerticalScale(const double value)
{
  this->VerticalScale = value;
  this->Modified();
}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetVerticalBias(const double value)
{
  this->VerticalBias = value;
  this->Modified();
}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetOutputMode(const int mode)
{
  this->OutputMode = mode;
  this->Modified();
}

//----------------------------------------------------------------------------

int vtkNetCDFLFRicReader::GetNumberOfCellArrays()
{
  return this->CellFields.size();
}

//----------------------------------------------------------------------------

const char* vtkNetCDFLFRicReader::GetCellArrayName(const int index)
{
  // Using a map to store fields is convenient for most purposes, but
  // ParaView also wants us to provide array names by index. Input files
  // shouldn't contain too many arrays, so we'll simply count map keys.
  int counter = 0;
  std::map<std::string, DataField>::const_iterator it;
  for (it = this->CellFields.begin(); it != this->CellFields.end();
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

//----------------------------------------------------------------------------

int vtkNetCDFLFRicReader::GetCellArrayStatus(const char* name)
{
  int status;
  std::map<std::string, DataField>::const_iterator it = this->CellFields.find(name);
  if (it != this->CellFields.end())
  {
    status = static_cast<int>(it->second.active);
  }
  else
  {
    status = 0;
    vtkDebugMacro("GetCellArrayStatus: no array with name=" << name << endl);
  }
  return status;
}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetCellArrayStatus(const char* name, const int status)
{
  std::map<std::string, DataField>::iterator it = this->CellFields.find(name);
  if (it != this->CellFields.end())
  {
    it->second.active = static_cast<bool>(status);
    this->Modified();
  }
  else
  {
    // Ignore unknown names
    vtkDebugMacro("SetCellArrayStatus: no array with name=" << name << endl);
  }
}

//----------------------------------------------------------------------------

int vtkNetCDFLFRicReader::GetNumberOfPointArrays()
{
  return this->PointFields.size();
}

//----------------------------------------------------------------------------

const char* vtkNetCDFLFRicReader::GetPointArrayName(const int index)
{
  int counter = 0;
  std::map<std::string, DataField>::const_iterator it;
  for (it = this->PointFields.begin(); it != this->PointFields.end();
       it++, counter++)
  {
    if (counter == index)
    {
      return it->first.c_str();
    }
  }
  vtkDebugMacro("GetPointArrayName: out of range index=" << index << endl);
  return nullptr;

}

//----------------------------------------------------------------------------

int vtkNetCDFLFRicReader::GetPointArrayStatus(const char* name)
{
  int status;
  std::map<std::string, DataField>::const_iterator it = this->PointFields.find(name);
  if (it != this->PointFields.end())
  {
    status = static_cast<int>(it->second.active);
  }
  else
  {
    status = 0;
    vtkDebugMacro("GetPointArrayStatus: no array with name=" << name << endl);
  }
  return status;
}

//----------------------------------------------------------------------------

void vtkNetCDFLFRicReader::SetPointArrayStatus(const char* name, const int status)
{
  std::map<std::string, DataField>::iterator it = this->PointFields.find(name);
  if (it != this->PointFields.end())
  {
    it->second.active = static_cast<bool>(status);
    this->Modified();
  }
  else
  {
    // Ignore unknown names
    vtkDebugMacro("SetPointArrayStatus: no array with name=" << name << endl);
  }
}
