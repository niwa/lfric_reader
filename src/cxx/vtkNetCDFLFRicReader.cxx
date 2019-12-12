#include "vtkNetCDFLFRicReader.h"
#include "netCDFLFRicReaderUtils.h"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
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
  this->Fields.clear();
  this->TimeSteps.clear();
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
  if (mesh2D.numTopologies == 0)
  {
    vtkErrorMacro("Failed to determine 2D UGRID mesh description." << endl);
    return 0;
  }

  // Look for vertical axis
  this->zAxis = inputFile.GetZAxisDescription(this->mesh2D.isLFRicXIOSFile,
                                              this->mesh2D.meshType);

  // Look for time axis
  this->tAxis = inputFile.GetTAxisDescription();

  //
  // Look for data fields
  //

  // Get field variable names, ignoring UGRID mesh definitions
  // and other variables
  for (std::string const &varName : inputFile.GetVarNames())
  {
    vtkDebugMacro("Request Information: Considering variable name=" << varName);

    // There is no attribute that uniquely identifies fields, so
    // check if this combination of attributes exists
    bool hasFieldAtts = (inputFile.VarHasAtt(varName, "standard_name") or
			 inputFile.VarHasAtt(varName, "long_name")) and
                        inputFile.VarHasAtt(varName, "mesh");

    // At the moment, variables can have up to 3 dimensions (time, vertical dimension,
    // and horizontal unstructured dimension)
    bool hasValidNumDims = inputFile.GetVarNumDims(varName) < 4;

    if (hasFieldAtts and hasValidNumDims)
    {
      // If field is not in list, insert and default to "don't load"
      std::map<std::string,bool>::const_iterator it = this->Fields.find(varName);
      if (it == this->Fields.end())
      {
        this->Fields.insert(it, std::pair<std::string,bool>(varName,false));
        vtkDebugMacro("=> inserted field" << endl);
      }
    }
  }
  vtkDebugMacro("Number of data arrays found=" << this->Fields.size() << endl);

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
                                                             {0}, {mesh2D.numNodes});
  std::vector<double> node_coords_y = inputFile.GetVarDouble(this->mesh2D.nodeCoordYVar,
                                                             {0}, {mesh2D.numNodes});
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
  // Determine number of edges
  //

  if (inputFile.VarHasAtt(this->mesh2D.meshTopologyVar, "edge_coordinates"))
  {
    const std::vector<std::string> edgeCoordVarNames =
      inputFile.GetAttTextSplit(this->mesh2D.meshTopologyVar, "edge_coordinates");

    // Assume that number of edges is the first dimension
    const std::string edgeDimName = inputFile.GetVarDimName(edgeCoordVarNames[0], 0);
    this->NumberOfEdges2D = inputFile.GetDimLen(edgeDimName);
    vtkDebugMacro("edgeDimName=" << edgeDimName << " NumberOfEdges2D=" <<
                  this->NumberOfEdges2D << endl);
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
  else if (this->mesh2D.meshType == fullLevelFace)
  {
    // Load vertex heights from file
    levels = inputFile.GetVarDouble(this->zAxis.axisVar,
                                    {startLevel}, {numLevels+1});
  }
  else if (this->mesh2D.meshType == halfLevelFace)
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

  if (!this->UseCartCoords)
  {
    vtkDebugMacro("Resolving grid periodicity..." << endl);

    resolvePeriodicGrid(node_coords_x, node_coords_y, face_nodes,
                        this->mesh2D.numFaces, this->mesh2D.numVertsPerFace);
    // Reset nnodes to new count
    this->mesh2D.numNodes = node_coords_x.size();
  }

  this->UpdateProgress(0.5);

  //
  // Construct VTK grid points
  //

  vtkDebugMacro("Setting VTK points..." << endl);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(this->mesh2D.numNodes*(numLevels+1));
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
                                      this->mesh2D.numVertsPerFace, this->mesh2D.numNodes);
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

      std::vector<double> read_buffer;
      std::vector<size_t> start;
      std::vector<size_t> count;

      read_buffer.clear();
      start.clear();
      count.clear();

      // Add time dimension to slice arrays if it exists
      if (inputFile.VarHasDim(varName, this->tAxis.axisDim))
      {
        if (inputFile.GetVarDimName(varName, 0) == this->tAxis.axisDim)
	{
          start.push_back(timestep);
          count.push_back(1);
          vtkDebugMacro("Found time dimension" << endl);
        }
        else
	{
          vtkWarningMacro("Time dimension must be first. Skipping this variable..." << endl);
        }
      }
      else
      {
        vtkDebugMacro("Found NO time dimension" << endl);
      }

      // Add vertical dimension to slice arrays if it exists
      // Note that there may be several vertical axes, so
      // these names need to be partially hard-coded
      bool hasVertDim = false;
      if (inputFile.VarHasDim(varName, "half_levels") or
          inputFile.VarHasDim(varName, "full_levels") or
          inputFile.VarHasDim(varName, this->zAxis.axisDim))
      {
        start.push_back(startLevel);
        count.push_back(numLevels);
        hasVertDim = true;
        vtkDebugMacro("Found vertical dimension" << endl);
      }
      else
      {
        vtkDebugMacro("Found NO vertical dimension" << endl);
      }

      // Find out which mesh type is used for this variable, we already
      // know that this variable has the "mesh" attribute
      mesh_types mesh_type;
      const std::string mesh_name = inputFile.GetAttText(varName, "mesh");

      // Assume half-level mesh by default, unless "Mesh2d_full_level" mesh is
      // present, in which case there will be several mesh types
      if ( this->mesh2D.meshTopologyVar != "Mesh2d_full_levels" or
           mesh_name == "Mesh2d_half_levels" )
      {
        mesh_type = half_level_face;
        vtkDebugMacro("Field is defined on half level face mesh" << endl);

        // Add horizontal index dimension to slice arrays - always exists
        start.push_back(0);
        count.push_back(this->mesh2D.numFaces);

        // Only call read method if we expect any data for this piece (partition)
        // Always the case for a 3D field, and if piece contains surface layer for
        // a 2D field
        if (hasVertDim or (startLevel == 0))
        {
          read_buffer = inputFile.GetVarDouble(varName, start, count);
        }

        // Fill remaining space with NaNs, this will only apply to 2D fields
        if (not hasVertDim)
        {
          read_buffer.resize(numLevels*this->mesh2D.numFaces,
                             std::numeric_limits<double>::quiet_NaN());
        }
      }
      else if ( mesh_name == "Mesh2d_edge_half_levels" )
      {
        mesh_type = half_level_edge;
        vtkDebugMacro("Field is defined on half level edge mesh" << endl);

        start.push_back(0);
        count.push_back(this->NumberOfEdges2D);

        if (hasVertDim or (startLevel == 0))
        {
          read_buffer = inputFile.GetVarDouble(varName, start, count);
        }

        if (not hasVertDim)
        {
          read_buffer.resize(numLevels*this->NumberOfEdges2D,
                             std::numeric_limits<double>::quiet_NaN());
        }
      }
      else if ( mesh_name == "Mesh2d_full_levels" )
      {
        mesh_type = full_level_face;
        vtkDebugMacro("Field is defined on full level face mesh" << endl);

        // Need to read extra level if field is 3D
        if (hasVertDim)
        {
          count.back() += 1;
        }

        start.push_back(0);
        count.push_back(this->mesh2D.numFaces);

        if (hasVertDim or (startLevel == 0))
        {
          read_buffer = inputFile.GetVarDouble(varName, start, count);
        }

        if (not hasVertDim)
        {
          read_buffer.resize((numLevels+1)*this->mesh2D.numFaces,
                             std::numeric_limits<double>::quiet_NaN());
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
      dataField->SetNumberOfTuples(numLevels*this->mesh2D.numFaces);
      dataField->SetName(varName.c_str());

      switch(mesh_type)
      {
        case half_level_face :
          vtkDebugMacro("half level face mesh: no interpolation needed" << endl);
          for (size_t i = 0; i < numLevels*this->mesh2D.numFaces; i++)
          {
            dataField->SetComponent(static_cast<vtkIdType>(i), 0, read_buffer[i]);
          }
          break;
        case full_level_face :
          vtkDebugMacro("full level face mesh: averaging top and bottom faces" << endl);
          for (size_t i = 0; i < numLevels*this->mesh2D.numFaces; i++)
          {
            const double bottomval = read_buffer[i];
            const double topval = read_buffer[i+this->mesh2D.numFaces];
            dataField->SetComponent(static_cast<vtkIdType>(i), 0, 0.5*(bottomval+topval));
          }
          break;
        case half_level_edge:
          vtkWarningMacro("WARNING: edge-centered fields cannot be handled correctly at the moment" << endl);
          vtkDebugMacro("half level edge mesh: averaging four edges" << endl);
          dataField->Fill(std::numeric_limits<double>::quiet_NaN());

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
          dataField->Fill(std::numeric_limits<double>::quiet_NaN());
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

int vtkNetCDFLFRicReader::GetNumberOfCellArrays()
{
  return this->Fields.size();
}

//----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------

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
