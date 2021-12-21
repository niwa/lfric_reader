#include "netCDFLFRicFile.h"
#include <vtk_netcdf.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <limits>

#define ncErrorMacro(e) if (e != NC_NOERR) {std::cerr << "NetCDF Error: " << nc_strerror(e) << "\n";}

#define errorMacro(x) std::cerr << x
#ifdef DEBUG
#define debugMacro(x) std::cerr << x
#else
#define debugMacro(x)
#endif

//----------------------------------------------------------------------------

netCDFLFRicFile::netCDFLFRicFile(const char* fileName)
{
  this->FileName = fileName;
  const int err = nc_open(fileName, NC_NOWRITE, &this->ncId);
  ncErrorMacro(err);
  if (err != NC_NOERR)
  {
    this->fileOpen = false;
  }
  else
  {
    this->fileOpen = true;
  }
}

//----------------------------------------------------------------------------

netCDFLFRicFile::~netCDFLFRicFile()
{
  if (this->fileOpen)
  {
    ncErrorMacro(nc_close(this->ncId));
    this->fileOpen = false;
  }
  this->FileName.clear();
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::IsFileOpen()
{
  return this->fileOpen;
}

//----------------------------------------------------------------------------

const char* netCDFLFRicFile::GetFileName()
{
  return this->FileName.c_str();
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::HasDim(const std::string& dimName)
{
  int dimId;
  bool retval = false;
  const int result = nc_inq_dimid(this->ncId, dimName.c_str(), &dimId);
  if (result == NC_NOERR)
  {
    retval=true;
  }
  else if (result != NC_EBADDIM)
  {
    ncErrorMacro(result);
  }
  return retval;
}

//----------------------------------------------------------------------------

size_t netCDFLFRicFile::GetDimLen(const int dimId)
{
  size_t dimLen;
  ncErrorMacro(nc_inq_dimlen(this->ncId, dimId, &dimLen));
  return dimLen;
}

//----------------------------------------------------------------------------

int netCDFLFRicFile::GetDimId(const std::string& dimName)
{
  int dimId;
  ncErrorMacro(nc_inq_dimid(this->ncId, dimName.c_str(), &dimId));
  return dimId;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetDimName(const int dimId)
{
  // NetCDF dimensions can use up to NC_MAX_NAME bytes
  char *dimName = new char[NC_MAX_NAME+1];
  ncErrorMacro(nc_inq_dimname(this->ncId, dimId, dimName));
  const std::string dimNameStr(dimName);
  delete[] dimName;
  return dimNameStr;
}

//----------------------------------------------------------------------------

int netCDFLFRicFile::GetVarId(const std::string& varName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));
  return varId;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetVarName(const int varId)
{
  // Valid netCDF IDs start at 0
  if (varId >= 0)
  {
    // NetCDF variables can use up to NC_MAX_NAME bytes
    char *varName = new char[NC_MAX_NAME+1];
    ncErrorMacro(nc_inq_varname(this->ncId, varId, varName));
    const std::string varNameStr(varName);
    delete[] varName;
    return varNameStr;
  }
  else
  {
    const std::string varNameStr("invalid varId");
    return varNameStr;
  }
}

//----------------------------------------------------------------------------

size_t netCDFLFRicFile::GetVarNumDims(const int varId)
{
  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));
  if (numDims < 0)
  {
    errorMacro("netCDFLFRicFile::GetVarNumDims: NetCDF returned invalid " <<
               "number of dimensions: " << numDims);
    return 0;
  }
  else
  {
    return static_cast<size_t>(numDims);
  }
}

//----------------------------------------------------------------------------

int netCDFLFRicFile::GetVarDimId(const int varId,
                                 const size_t dim)
{
  const size_t numDims = this->GetVarNumDims(varId);

  // Return undefined ID if out of range
  int dimId = -1;
  if (dim < numDims)
  {
    std::vector<int> dimIds;
    dimIds.resize(numDims);
    ncErrorMacro(nc_inq_vardimid(this->ncId, varId, dimIds.data()));
    dimId = dimIds[dim];
  }

  return dimId;
}

//----------------------------------------------------------------------------

int netCDFLFRicFile::GetAttInt(const int varId, const std::string& attName)
{
  int attInt;
  ncErrorMacro(nc_get_att_int(this->ncId, varId, attName.c_str(), &attInt));
  return attInt;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetAttText(const int varId, const std::string& attName)
{
  size_t attTextLen;
  ncErrorMacro(nc_inq_attlen(this->ncId, varId, attName.c_str(), &attTextLen));

  char *attText = new char[attTextLen+1];
  ncErrorMacro(nc_get_att_text(this->ncId, varId, attName.c_str(), attText));

  // We cannot assume that attText is always null-terminated, so copy only attTextLen bytes
  std::string attTextStr(attText, attTextLen);
  delete[] attText;

  return attTextStr;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetAttText(const std::string& varName, const std::string& attName)
{
  const int varId = this->GetVarId(varName);
  return this->GetAttText(varId, attName);
}

//----------------------------------------------------------------------------

std::vector<std::string> netCDFLFRicFile::GetAttTextSplit(const int varId,
                                                          const std::string& attName)
{
  const std::string attText = this->GetAttText(varId, attName);
  std::istringstream attStream(attText);
  std::vector<std::string> textSplit(std::istream_iterator<std::string>{attStream},
                                     std::istream_iterator<std::string>());
  return textSplit;
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::HasVar(const std::string& varName)
{
  int varId;
  bool retval = false;
  const int result = nc_inq_varid(this->ncId, varName.c_str(), &varId);
  if (result == NC_NOERR)
  {
    retval=true;
  }
  else if (result != NC_ENOTVAR)
  {
    ncErrorMacro(result);
  }
  return retval;
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::VarHasAtt(const int varId, const std::string& attName)
{
  int attId;
  const int result = nc_inq_attid(this->ncId, varId, attName.c_str(), &attId);
  if (result == NC_NOERR)
  {
    return true;
  }
  else if (result == NC_ENOTATT)
  {
    return false;
  }
  else
  {
    ncErrorMacro(result);
  }
  return false;
}

//----------------------------------------------------------------------------

size_t netCDFLFRicFile::GetNumVars()
{
  int numVars;
  ncErrorMacro(nc_inq_nvars(this->ncId, &numVars));
  if (numVars < 0)
  {
    errorMacro("netCDFLFRicFile::GetNumVars: NetCDF returned invalid " <<
               "number of variables: " << numVars);
    return 0;
  }
  else
  {
    return static_cast<size_t>(numVars);
  }
}

//----------------------------------------------------------------------------

void netCDFLFRicFile::LoadVarDouble(const int varId,
                                    const std::vector<size_t>& start,
                                    const std::vector<size_t>& count,
                                    double* buffer)
{
  // Check if number of dimensions matches netCDF variable
  const size_t numDims = this->GetVarNumDims(varId);
  if (numDims != start.size() || numDims != count.size())
  {
    std::cerr << "netCDFLFRicFile::LoadVarDouble: number of dimensions does not match netCDF variable.\n";
  }
  else
  {
    // NetCDF will automatically convert non-double numeric data into double
    // This function will also check index ranges automatically
    ncErrorMacro(nc_get_vara_double(this->ncId, varId, start.data(),
                                    count.data(), buffer));
  }

}

//----------------------------------------------------------------------------

std::vector<double> netCDFLFRicFile::GetVarDouble(
                    const int varId,
                    const std::vector<size_t>& start,
                    const std::vector<size_t>& count)
{
  // Compute total number of elements to read
  size_t size = 1;
  for (size_t const n : count)
  {
    size *= n;
  }

  std::vector<double> varData;
  varData.resize(size);

  // Check if number of dimensions matches netCDF variable
  const size_t numDims = this->GetVarNumDims(varId);
  if (numDims != start.size() || numDims != count.size())
  {
    std::cerr << "netCDFLFRicFile::GetVarDouble: number of dimensions does not match netCDF variable.\n";
    std::fill(varData.begin(), varData.end(), std::numeric_limits<double>::quiet_NaN());
  }
  else
  {
    // NetCDF will automatically convert non-double numeric data into double
    // This function will also check index ranges automatically
    ncErrorMacro(nc_get_vara_double(this->ncId, varId, start.data(),
                                    count.data(), varData.data()));
  }

  return varData;
}

//----------------------------------------------------------------------------

std::vector<long long> netCDFLFRicFile::GetVarLongLong(
                       const int varId,
                       const std::vector<size_t>& start,
                       const std::vector<size_t>& count)
{
  size_t size = 1;
  for (size_t const n : count)
  {
    size *= n;
  }

  std::vector<long long> varData;
  varData.resize(size);

  // Check if number of dimensions matches netCDF variable
  const size_t numDims = this->GetVarNumDims(varId);
  if (numDims != start.size() || numDims != count.size())
  {
    std::cerr << "netCDFLFRicFile::GetVarLongLong: number of dimensions does not match netCDF variable.\n";
    std::fill(varData.begin(), varData.end(), 0);
  }
  else
  {
    // netCDF will automatically convert non-double numeric data into long long
    ncErrorMacro(nc_get_vara_longlong(this->ncId, varId, start.data(),
                                      count.data(), varData.data()));
  }

  return varData;
}

//----------------------------------------------------------------------------

UGRIDMeshDescription netCDFLFRicFile::GetMesh2DDescription()
{
  UGRIDMeshDescription mesh2D = UGRIDMeshDescription();

  // LFRic output files can have a separate "half-level edge" mesh
  bool hasHalfLevelEdgeMesh = false;
  int meshTopologyVarEdgeId = -1;

  // Look for UGRID 2D mesh description dummy variables,
  // must have attribute cf_role=mesh_topology
  // There can be several mesh descriptions in a file, but only LFRic output
  // files are permitted to have multiple horizontal meshes
  for (int varId = 0; varId < static_cast<int>(this->GetNumVars()); varId++)
  {
    const std::string varName = this->GetVarName(varId);

    if (this->VarHasAtt(varId, "cf_role"))
    {
      if (this->GetAttText(varId, "cf_role") == "mesh_topology")
      {
        mesh2D.numTopologies++;

        // Mesh descriptions in LFRic output files start with "Mesh2d"
        // and can be separate, or partially or fully unified
        // topology_dimension=2 means faces
        if ((varName == "Mesh2d_full_levels" ||
             varName == "Mesh2d_face" ||
             (varName == "Mesh2d" && this->HasDim("full_levels"))) &&
            this->GetAttInt(varId, "topology_dimension") == 2)
        {
          mesh2D.meshTopologyVarId = varId;
          mesh2D.meshType = fullLevelFaceMesh;
          mesh2D.isLFRicXIOSFile = true;
        }
        // LFRic output files can come with separate edge mesh, still has
        // topology_dimension=2
        else if ((varName == "Mesh2d_edge_half_levels" ||
                  varName == "Mesh2d_edge") &&
                 this->GetAttInt(varId, "topology_dimension") == 2)
        {
          hasHalfLevelEdgeMesh = true;
          meshTopologyVarEdgeId = varId;
        }
        // Assume that mesh is half-level type otherwise (2D faces represent 3D volumes)
        // Adopt first mesh description
        else if (mesh2D.meshTopologyVarId < 0 &&
                 this->GetAttInt(varId, "topology_dimension") == 2)
        {
          mesh2D.meshTopologyVarId = varId;
          mesh2D.meshType = halfLevelFaceMesh;
          mesh2D.isLFRicXIOSFile = false;
        }
      }
    }
  }

  // Must have at least one face mesh
  if (mesh2D.meshTopologyVarId < 0 || mesh2D.numTopologies == 0)
  {
    std::cerr << "GetMesh2DDescription: At least one UGRID face mesh is required.\n";
    return UGRIDMeshDescription();
  }

  // Require that a multi-mesh file is an LFRic file - the current implementation
  // only supports a single UGRID mesh per file otherwise
  if (!mesh2D.isLFRicXIOSFile && mesh2D.numTopologies > 1)
  {
    std::cerr << "GetMesh2DDescription: Only LFRic output files can have multiple UGRID meshes.\n";
    return UGRIDMeshDescription();
  }

  // Accept up to 3 meshes in an LFRic file
  if (mesh2D.isLFRicXIOSFile && mesh2D.numTopologies > 3)
  {
    std::cerr << "GetMesh2DDescription: LFRic output files can only have up to 3 UGRID meshes.\n";
    return UGRIDMeshDescription();
  }

  //
  // Get node coordinate variables and their dimensions
  //

  if (!this->VarHasAtt(mesh2D.meshTopologyVarId, "node_coordinates"))
  {
    std::cerr << "GetMesh2DDescription: UGRID topology variable must have node_coordinates attribute.\n";
    return UGRIDMeshDescription();
  }

  std::vector<std::string> nodeCoordVarNames =
    this->GetAttTextSplit(mesh2D.meshTopologyVarId, "node_coordinates");

  if (nodeCoordVarNames.size() != 2)
  {
    std::cerr << "GetMesh2DDescription: Expected 2 node coordinate variables but received " <<
                  nodeCoordVarNames.size() << ".\n";
    return UGRIDMeshDescription();
  }

  // Longitude must be the "x axis" to detect cubed-sphere grids
  if (this->GetAttText(nodeCoordVarNames[0], "standard_name") == "latitude" ||
      this->GetAttText(nodeCoordVarNames[0], "standard_name") == "projection_y_coordinate")
  {
    nodeCoordVarNames[0].swap(nodeCoordVarNames[1]);
  }

  if ((this->GetAttText(nodeCoordVarNames[0], "standard_name") != "longitude" &&
       this->GetAttText(nodeCoordVarNames[0], "standard_name") != "projection_x_coordinate") ||
      (this->GetAttText(nodeCoordVarNames[1], "standard_name") != "latitude" &&
       this->GetAttText(nodeCoordVarNames[1], "standard_name") != "projection_y_coordinate"))
  {
    std::cerr << "GetMesh2DDescription: Node coord variables must be named latitude/longitude " <<
      "or projection_x_coordinate/projection_y_coordinate.\n";
    return UGRIDMeshDescription();
  }

  mesh2D.nodeCoordXVarId = this->GetVarId(nodeCoordVarNames[0]);
  mesh2D.nodeCoordYVarId = this->GetVarId(nodeCoordVarNames[1]);

  // Planar LAMs use "m" as location units
  if (this->GetAttText(mesh2D.nodeCoordXVarId, "units") == "m")
  {
    mesh2D.isPlanarLAM = true;
  }

  // Get node dimension name and length
  mesh2D.nodeDimId = this->GetVarDimId(mesh2D.nodeCoordXVarId, 0);
  mesh2D.numNodes = this->GetDimLen(mesh2D.nodeDimId);

  //
  // Get face-node connectivity variable and its dimensions
  //

  if (!this->VarHasAtt(mesh2D.meshTopologyVarId, "face_node_connectivity"))
  {
    std::cerr << "GetMesh2DDescription: UGRID topology variable must have face_node_connectivity attribute.\n";
    return UGRIDMeshDescription();
  }

  const std::string faceNodeConnVar = this->GetAttText(mesh2D.meshTopologyVarId,
                                      "face_node_connectivity");
  mesh2D.faceNodeConnVarId = this->GetVarId(faceNodeConnVar);

  // Number of faces must be the first dimension, although
  // the UGRID conventions are more lenient
  mesh2D.faceDimId = this->GetVarDimId(mesh2D.faceNodeConnVarId, 0);

  if (this->VarHasAtt(mesh2D.meshTopologyVarId, "face_dimension"))
  {
    if (mesh2D.faceDimId !=
        this->GetDimId(this->GetAttText(mesh2D.meshTopologyVarId, "face_dimension")))
    {
      std::cerr << "GetMesh2DDescription: face_node_connectivity must have face dimension first.\n";
      return UGRIDMeshDescription();
    }
  }

  mesh2D.numFaces = this->GetDimLen(mesh2D.faceDimId);

  // Assume that number of vertices per face is the second dimension
  mesh2D.vertsPerFaceDimId = this->GetVarDimId(mesh2D.faceNodeConnVarId, 1);
  mesh2D.numVertsPerFace = this->GetDimLen(mesh2D.vertsPerFaceDimId);

  // Correction for non-zero start index
  if (this->VarHasAtt(mesh2D.faceNodeConnVarId, "start_index"))
  {
    mesh2D.faceNodeStartIdx = static_cast<long long>(
      this->GetAttInt(mesh2D.faceNodeConnVarId, "start_index"));
  }

  // Alternative face dimension ID for old LFRic output files, required
  // for identifying field variables
  if (mesh2D.isLFRicXIOSFile && this->HasDim("nMesh2d_half_levels_face"))
  {
    mesh2D.faceDimIdAlt = this->GetDimId("nMesh2d_half_levels_face");
  }

  //
  // Get edge coordinate variables and their dimensions (if available)
  //

  if (hasHalfLevelEdgeMesh || this->VarHasAtt(mesh2D.meshTopologyVarId, "face_edge_connectivity"))
  {
    // Prioritise half-level edge mesh if it exists - the face mesh may have face-edge connectivity, too,
    // but it will not be useable for edge-centered data if the meshes are not unified in LFRic output
    const int lookupMeshTopologyVarId = hasHalfLevelEdgeMesh ? meshTopologyVarEdgeId : mesh2D.meshTopologyVarId;

    if (!this->VarHasAtt(lookupMeshTopologyVarId, "edge_coordinates"))
    {
      std::cerr << "GetMesh2DDescription: UGRID edge topology variable must have edge_coordinates attribute.\n";
      return UGRIDMeshDescription();
    }

    std::vector<std::string> edgeCoordVarNames =
      this->GetAttTextSplit(lookupMeshTopologyVarId, "edge_coordinates");

    if (edgeCoordVarNames.size() != 2)
    {
      std::cerr << "GetMesh2DDescription: Expected 2 edge coordinate variables but received " <<
                   edgeCoordVarNames.size() << ".\n";
      return UGRIDMeshDescription();
    }

    // Longitude must be the "x axis" for consistency with face mesh
    if (this->GetAttText(edgeCoordVarNames[0], "standard_name") == "latitude")
    {
      edgeCoordVarNames[0].swap(edgeCoordVarNames[1]);
    }

  if ((this->GetAttText(edgeCoordVarNames[0], "standard_name") != "longitude" &&
       this->GetAttText(edgeCoordVarNames[0], "standard_name") != "projection_x_coordinate") ||
      (this->GetAttText(edgeCoordVarNames[1], "standard_name") != "latitude" &&
       this->GetAttText(edgeCoordVarNames[1], "standard_name") != "projection_y_coordinate"))
    {
      std::cerr << "GetMesh2DDescription: Edge coord variables must be named latitude/longitude " <<
      "or projection_x_coordinate/projection_y_coordinate.\n";
      return UGRIDMeshDescription();
    }

    mesh2D.edgeCoordXVarId = this->GetVarId(edgeCoordVarNames[0]);
    mesh2D.edgeCoordYVarId = this->GetVarId(edgeCoordVarNames[1]);

    // Get edge dimension name and length
    mesh2D.edgeDimId = this->GetVarDimId(mesh2D.edgeCoordXVarId, 0);
    mesh2D.numEdges = this->GetDimLen(mesh2D.edgeDimId);

    // Get face-edge connectivity if the mesh is unified
    if (!hasHalfLevelEdgeMesh)
    {
      const std::string faceEdgeConnVar = this->GetAttText(mesh2D.meshTopologyVarId,
                                                           "face_edge_connectivity");
      mesh2D.faceEdgeConnVarId = this->GetVarId(faceEdgeConnVar);

      if (mesh2D.faceDimId != this->GetVarDimId(mesh2D.faceEdgeConnVarId, 0))
      {
        std::cerr << "GetMesh2DDescription: face_edge_connectivity must have face dimension first.\n";
        return UGRIDMeshDescription();
      }

      mesh2D.edgesPerFaceDimId = this->GetVarDimId(mesh2D.faceEdgeConnVarId, 1);
      mesh2D.numEdgesPerFace = this->GetDimLen(mesh2D.edgesPerFaceDimId);

      if (this->VarHasAtt(mesh2D.faceEdgeConnVarId, "start_index"))
      {
        mesh2D.faceEdgeStartIdx = static_cast<long long>(
          this->GetAttInt(mesh2D.faceEdgeConnVarId, "start_index"));
      }
    }
  }

  return mesh2D;
}

//----------------------------------------------------------------------------

std::map<std::string, CFAxis> netCDFLFRicFile::GetZAxisDescription(const bool isLFRicXIOSFile,
                                                                   const mesh2DTypes meshType)
{
  std::map<std::string, CFAxis> verticalAxes;

  // Workaround for LFRic XIOS output files where vertical axes do not have
  // attributes required by CF convention.
  // LFRic files typically have 2 vertical axes, "full_levels" (cell interfaces) and
  // "half_levels" (cell mid-levels). We need to store these axes for accessing the data,
  // and create a third axis ("vtk") that counts the number vertical levels, for creating
  // the VTK grid.
  if (isLFRicXIOSFile && (meshType == fullLevelFaceMesh || meshType == halfLevelFaceMesh))
  {
    if (this->HasVar("full_levels"))
    {
      CFAxis levels = CFAxis();
      levels.axisVarId = this->GetVarId("full_levels");
      levels.axisDimId = this->GetVarDimId(levels.axisVarId, 0);
      levels.axisLength = this->GetDimLen(levels.axisDimId);
      verticalAxes.insert(std::pair<std::string, CFAxis>("full_levels",levels));
      if (meshType == fullLevelFaceMesh)
      {
        verticalAxes.insert(std::pair<std::string, CFAxis>("vtk",levels));
        // Need the number of cells here, "full_levels" are interfaces between cells
        verticalAxes.at("vtk").axisLength -= 1;
      }
    }
    if (this->HasVar("half_levels"))
    {
      CFAxis levels = CFAxis();
      levels.axisVarId = this->GetVarId("half_levels");
      levels.axisDimId = this->GetVarDimId(levels.axisVarId, 0);
      levels.axisLength = this->GetDimLen(levels.axisDimId);
      verticalAxes.insert(std::pair<std::string, CFAxis>("half_levels",levels));
      if (meshType == halfLevelFaceMesh)
      {
        verticalAxes.insert(std::pair<std::string, CFAxis>("vtk",levels));
      }
    }
  }
  // Assume that other input files are CF-compliant and use "positive" attribute
  // Restrict choices to "level_height" and "model_level_number"
  else
  {
    CFAxis levels = CFAxis();
    for (int varId = 0; varId < static_cast<int>(this->GetNumVars()); varId++)
    {
      if (this->VarHasAtt(varId, "positive"))
      {
        // Prefer to use level_height variable if it exists
        const std::string varName = this->GetVarName(varId);
	if (varName == "level_height" ||
            (levels.axisVarId < 0 && varName == "model_level_number"))
        {
          levels.axisVarId = varId;
          levels.axisDimId = this->GetVarDimId(levels.axisVarId, 0);
          levels.axisLength = this->GetDimLen(levels.axisDimId);
          // Treat as "half-levels" axis by default
          verticalAxes.insert(std::pair<std::string, CFAxis>("half_levels",levels));
          verticalAxes.insert(std::pair<std::string, CFAxis>("vtk",levels));
        }
      }
    }
  }

  // Assume 2D-only file and set vertical axis to single level if no axis found,
  // the data fields will not have a vertical dimension
  if (verticalAxes.count("vtk") == 0)
  {
    CFAxis levels = CFAxis();
    levels.axisLength = 1;
    verticalAxes.insert(std::pair<std::string, CFAxis>("vtk",levels));
  }

  return verticalAxes;
}

//----------------------------------------------------------------------------

CFAxis netCDFLFRicFile::GetTAxisDescription()
{
  CFAxis time = CFAxis();

  // Look for variable with standard_name = time, even though this is not
  // strictly required by CF conventions
  for (int varId = 0; varId < static_cast<int>(this->GetNumVars()); varId++)
  {
    if (this->VarHasAtt(varId, "standard_name"))
    {
      if (this->GetAttText(varId, "standard_name") == "time")
      {
        time.axisVarId = varId;
        time.axisDimId = this->GetVarDimId(time.axisVarId, 0);
        time.axisLength = this->GetDimLen(time.axisDimId);
      }
    }
  }

  return time;
}

//----------------------------------------------------------------------------

// Discovers data fields in netCDF file using CF attributes and adds them to
// cell data and point data maps if all variable dimensions can be identified
void netCDFLFRicFile::UpdateFieldMaps(const UGRIDMeshDescription & mesh2D,
                                      const std::map<std::string, CFAxis> & zAxes,
                                      const CFAxis & tAxis,
                                      std::map<std::string, DataField> & CellFields,
                                      std::map<std::string, DataField> & PointFields)
{
  debugMacro("Entering netCDFLFRicFile::UpdateFieldMaps...\n");

  // There can be up to two vertical dimenions
  const int halfLevelsDimId = zAxes.count("half_levels") > 0 ? zAxes.at("half_levels").axisDimId : -1;
  const int fullLevelsDimId = zAxes.count("full_levels") > 0 ? zAxes.at("full_levels").axisDimId : -1;

  for (int varId = 0; varId < static_cast<int>(this->GetNumVars()); varId++)
  {
    debugMacro("Considering variable " << this->GetVarName(varId) << "\n");

    // Require these CF-netCDF attributes to distinguish fields from UGRID
    // variables and other data
    bool valid = (this->VarHasAtt(varId, "standard_name") ||
                  this->VarHasAtt(varId, "long_name")) &&
                  this->VarHasAtt(varId, "mesh") &&
                  this->VarHasAtt(varId, "coordinates") &&
                  this->VarHasAtt(varId, "location");

    // Identify and record variable dimensions, and add field to map if all
    // dimensions have been identified (allowing for up to one unidentified
    // one for multi-component fields)
    if (valid)
    {
      debugMacro("Variable has required netCDF attributes\n");

      DataField fieldSpec = DataField();

      const size_t numDims = this->GetVarNumDims(varId);
      fieldSpec.dims.resize(numDims);

      // Determine stride in flat data array for each dimension to
      // flexibly recover the data in any dimension order
      size_t stride = 1;

      // Try to identify all dimensions - up to one dimension may remain
      // unidentified and will be treated as a component dimension
      size_t numDimsIdentified = 0;

      // Loop over variable dims backwards to simplify computing strides
      for (size_t iDim = numDims-1; iDim < numDims; iDim--)
      {
        const int thisDimId = this->GetVarDimId(varId, iDim);
        bool isIdentified = false;

        // Try to identify dim
        if (thisDimId == mesh2D.edgeDimId)
        {
          fieldSpec.hasHorizontalDim = true;
          fieldSpec.meshType = halfLevelEdgeMesh;
          fieldSpec.dims[iDim].dimType = horizontalAxisDim;
          isIdentified = true;
          debugMacro("Found horizontal edge dim\n");
        }
        // Older LFRic files can have two different face meshes
        else if (thisDimId == mesh2D.faceDimId || thisDimId == mesh2D.faceDimIdAlt)
        {
          fieldSpec.hasHorizontalDim = true;
          // Assume half-level faces for now, if meshType has not been determined yet
          if (fieldSpec.meshType == unknownMesh) fieldSpec.meshType = halfLevelFaceMesh;
          fieldSpec.dims[iDim].dimType = horizontalAxisDim;
          isIdentified = true;
          debugMacro("Found horizontal face dim\n");
        }
        else if (thisDimId == fullLevelsDimId || thisDimId == halfLevelsDimId)
        {
          fieldSpec.hasVerticalDim = true;
          // Only face-type fields can be defined at full level (cell interfaces)
          if (thisDimId == fullLevelsDimId) fieldSpec.meshType = fullLevelFaceMesh;
          fieldSpec.dims[iDim].dimType = verticalAxisDim;
          isIdentified = true;
          debugMacro("Found vertical dim\n");
        }
        else if (thisDimId == tAxis.axisDimId)
        {
          fieldSpec.hasTimeDim = true;
          fieldSpec.dims[iDim].dimType = timeAxisDim;
          isIdentified = true;
          debugMacro("Found time dim\n");
        }
        // If dimension is not known, assume it is a component dimension if it
        // has up to 9 components and no other component dimension has been set
        else if (this->GetDimLen(thisDimId) < 10 && !fieldSpec.hasComponentDim)
        {
          fieldSpec.hasComponentDim = true;
          fieldSpec.dims[iDim].dimType = componentAxisDim;
          isIdentified = true;
          debugMacro("Found component dim\n");
        }
        else
        {
          debugMacro("Could not identify dim " << this->GetDimName(thisDimId) <<
                     " and it cannot be accepted as component dim\n");
        }

        // Record dimension length and stride for this dimension
        if (isIdentified)
        {
          fieldSpec.dims[iDim].dimLength = this->GetDimLen(thisDimId);
          fieldSpec.dims[iDim].dimStride = stride;
          stride *= fieldSpec.dims[iDim].dimLength;
          numDimsIdentified++;
          debugMacro("Dim has length " << fieldSpec.dims[iDim].dimLength <<
                     " and stride " << fieldSpec.dims[iDim].dimStride << "\n");
        }
      }

      debugMacro("Identified " << numDimsIdentified << " of " <<
                 numDims << " dimensions\n");

      // Fields must have a horizontal dimension and all dimensions must be
      // accounted for
      valid &= fieldSpec.hasHorizontalDim;
      valid &= (numDimsIdentified == numDims);

      if (valid)
      {
        const std::string varName = this->GetVarName(varId);

        // Make W2 fields, which are defined on half-level edge meshes
        // available as point data
        if (fieldSpec.meshType == halfLevelEdgeMesh)
        {
          fieldSpec.location = pointFieldLoc;
          // If field is not in list, insert and default to "don't load"
          std::map<std::string, DataField>::const_iterator it = PointFields.find(varName);
          if (it == PointFields.end())
          {
            debugMacro("Field is valid and is not in PointFields list, inserting...\n");
            PointFields.insert(it, std::pair<std::string, DataField>(varName, fieldSpec));
          }
        }
        // Make all fields available as cell data
        // W2 fields require face-edge connectivity data
        if (fieldSpec.meshType != halfLevelEdgeMesh || mesh2D.faceEdgeConnVarId > -1)
        {
          fieldSpec.location = cellFieldLoc;
          std::map<std::string, DataField>::const_iterator it = CellFields.find(varName);
          if (it == CellFields.end())
          {
            debugMacro("Field is valid and is not in CellFields list, inserting...\n");
            CellFields.insert(it, std::pair<std::string, DataField>(varName, fieldSpec));
          }
        }
      }
      else
      {
        debugMacro("Ignoring field without horizontal dim or unidentified dims\n");
      }
    }
  }
  debugMacro("Finished netCDFLFRicFile::UpdateFieldMaps.\n");
}
