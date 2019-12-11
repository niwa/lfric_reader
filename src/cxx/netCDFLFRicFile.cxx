#include "netCDFLFRicFile.h"
#include <vtk_netcdf.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <limits>

#define ncErrorMacro(e) if (e != NC_NOERR) {std::cerr << "NetCDF Error: " << nc_strerror(e) << "\n";}

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
  const int result = nc_inq_dimid(this->ncId, dimName.c_str(), &dimId);
  if (result == NC_NOERR)
  {
    return true;
  }
  else if (result == NC_EBADDIM)
  {
    return false;
  }
  else
  {
    ncErrorMacro(result);
  }
}

//----------------------------------------------------------------------------

size_t netCDFLFRicFile::GetDimLen(const std::string& dimName)
{
  int dimId;
  size_t dimLen;
  ncErrorMacro(nc_inq_dimid(this->ncId, dimName.c_str(), &dimId));
  ncErrorMacro(nc_inq_dimlen(this->ncId, dimId, &dimLen));
  return dimLen;
}

//----------------------------------------------------------------------------

size_t netCDFLFRicFile::GetVarNumDims(const std::string& varName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));

  // We rely on netCDF and assume that numDims >= 0
  return static_cast<size_t>(numDims);
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetVarDimName(const std::string& varName,
                                           const size_t dim)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));

  std::string dimNameStr;
  dimNameStr.clear();

  // Return empty string if out of range
  if (dim >= 0 and dim < numDims)
  {
    std::vector<int> dimIds;
    dimIds.resize(numDims);
    ncErrorMacro(nc_inq_vardimid(this->ncId, varId, dimIds.data()));

    char *dimName = new char[NC_MAX_NAME+1];
    ncErrorMacro(nc_inq_dimname(this->ncId, dimIds[dim], dimName));

    // NetCDF docs say that dimension names are guaranteed to be null-terminated
    dimNameStr = dimName;

    delete[] dimName;
  }

  return dimNameStr;
}

//----------------------------------------------------------------------------

int netCDFLFRicFile::GetAttInt(const std::string& varName, const std::string& attName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int attInt;
  ncErrorMacro(nc_get_att_int(this->ncId, varId, attName.c_str(), &attInt));

  return attInt;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetAttText(const std::string& varName, const std::string& attName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

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

std::vector<std::string> netCDFLFRicFile::GetAttTextSplit(const std::string& varName,
                                                          const std::string& attName)
{
  std::string attText = this->GetAttText(varName, attName);
  std::istringstream attStream(attText);
  std::vector<std::string> textSplit(std::istream_iterator<std::string>{attStream},
                                     std::istream_iterator<std::string>());
  return textSplit;
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::HasVar(const std::string& varName)
{
  int varId;
  const int result = nc_inq_varid(this->ncId, varName.c_str(), &varId);
  if (result == NC_NOERR)
  {
    return true;
  }
  else if (result == NC_ENOTVAR)
  {
    return false;
  }
  else
  {
    ncErrorMacro(result);
  }
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::VarHasDim(const std::string& varName, const std::string& dimName)
{
  // If dimName does not exist, the variable won't have it...
  if (not this->HasDim(dimName))
  {
    return false;
  }

  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));

  std::vector<int> dimIds;
  dimIds.resize(numDims);
  ncErrorMacro(nc_inq_vardimid(this->ncId, varId, dimIds.data()));

  int dimIdOfInterest;
  ncErrorMacro(nc_inq_dimid(this->ncId, dimName.c_str(), &dimIdOfInterest));

  bool dimFound = false;
  for (int const dimId : dimIds)
  {
    if (dimId == dimIdOfInterest) dimFound = true;
  }

  return dimFound;
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::VarHasAtt(const std::string& varName, const std::string& attName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

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
}

//----------------------------------------------------------------------------

std::vector<std::string> netCDFLFRicFile::GetVarNames()
{
  std::vector<std::string> varNames;
  int numVars;
  ncErrorMacro(nc_inq_nvars(this->ncId, &numVars));
  varNames.resize(numVars);

  // NetCDF variables can use up to NC_MAX_NAME bytes
  char *varName = new char[NC_MAX_NAME+1];
  for (int iVar = 0; iVar < numVars; iVar++)
  {
    ncErrorMacro(nc_inq_varname(this->ncId, iVar, varName));
    // NetCDF docs say that variable names are guaranteed to be null-terminated
    varNames[iVar] = varName;
  }
  delete[] varName;
  return varNames;
}

//----------------------------------------------------------------------------

std::vector<double> netCDFLFRicFile::GetVarDouble(
                    const std::string& varName,
                    const std::vector<size_t> start,
                    const std::vector<size_t> count)
{
  // Find variable by name
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  // Compute total number of elements to read
  size_t size = 1;
  for (size_t const n : count)
  {
    size *= n;
  }

  std::vector<double> varData;
  varData.resize(size);

  // Check if number of dimensions matches netCDF variable
  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));
  if (numDims != start.size() or numDims != count.size())
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
                       const std::string& varName,
                       const std::vector<size_t> start,
                       const std::vector<size_t> count)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  size_t size = 1;
  for (size_t const n : count)
  {
    size *= n;
  }

  std::vector<long long> varData;
  varData.resize(size);

  // Check if number of dimensions matches netCDF variable
  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));
  if (numDims != start.size() or numDims != count.size())
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

  // Look for UGRID 2D mesh description dummy variable,
  // must have attribute cf_role=mesh_topology
  // There can be several mesh descriptions in a file
  mesh2D.meshTopologyVar.clear();
  for (std::string const &varName : this->GetVarNames())
  {
    if (this->VarHasAtt(varName, "cf_role"))
    {
      if (this->GetAttText(varName, "cf_role") == "mesh_topology")
      {
        mesh2D.numTopologies++;

        // Prefer LFRic full-level face mesh which matches VTK grids
        // topology_dimension=2 means faces
        if (varName == "Mesh2d_full_levels" and
            this->GetAttInt(varName, "topology_dimension") == 2)
	{
          mesh2D.meshTopologyVar = varName;
          mesh2D.meshType = fullLevelFace;
          // Workaround for non-CF-compliant vertical axis
          mesh2D.isLFRicXIOSFile = true;
        }
        // Assume that mesh is half-level type otherwise
        else if (mesh2D.meshTopologyVar.empty() and
                 this->GetAttInt(varName, "topology_dimension") == 2)
        {
          mesh2D.meshTopologyVar = varName;
          mesh2D.meshType = halfLevelFace;
          mesh2D.isLFRicXIOSFile = false;
        }
      }
    }
  }

  // Must have at least one face mesh
  if (mesh2D.meshTopologyVar.empty() or mesh2D.numTopologies == 0)
  {
    std::cerr << "GetMesh2DDescription: At least one UGRID face mesh is required.\n";
    return UGRIDMeshDescription();
  }

  // Require that a multi-mesh file is an LFRic file - we will not be able make
  // sense of multiple meshes otherwise
  if (not mesh2D.isLFRicXIOSFile and mesh2D.numTopologies > 1)
  {
    std::cerr << "GetMesh2DDescription: Only LFRic output files can have multiple UGRID meshes.\n";
    return UGRIDMeshDescription();
  }

  // Accept up to 3 meshes in an LFRic file
  if (mesh2D.isLFRicXIOSFile and mesh2D.numTopologies > 3)
  {
    std::cerr << "GetMesh2DDescription: LFRic output files can only have up to 3 UGRID meshes.\n";
    return UGRIDMeshDescription();
  }

  //
  // Get node coordinate variables and their dimensions
  //

  if (not this->VarHasAtt(mesh2D.meshTopologyVar, "node_coordinates"))
  {
    std::cerr << "GetMesh2DDescription: UGRID topology variable must have node_coordinates attribute.\n";
    return UGRIDMeshDescription();
  }

  std::vector<std::string> nodeCoordVarNames =
    this->GetAttTextSplit(mesh2D.meshTopologyVar, "node_coordinates");

  if (nodeCoordVarNames.size() != 2)
  {
    std::cerr << "GetMesh2DDescription: Expected 2 node coordinate variables but received " <<
                  nodeCoordVarNames.size() << ".\n";
    return UGRIDMeshDescription();
  }

  // Longitude must be the "x axis" to detect cubed-sphere grids
  if (this->GetAttText(nodeCoordVarNames[0], "standard_name") == "latitude")
  {
    nodeCoordVarNames[0].swap(nodeCoordVarNames[1]);
  }

  if (this->GetAttText(nodeCoordVarNames[0], "standard_name") != "longitude" or
      this->GetAttText(nodeCoordVarNames[1], "standard_name") != "latitude")
  {
    std::cerr << "GetMesh2DDescription: Node coord variables must be named latitude/longitude.\n";
    return UGRIDMeshDescription();
  }

  mesh2D.nodeCoordXVar = nodeCoordVarNames[0];
  mesh2D.nodeCoordYVar = nodeCoordVarNames[1];

  // Get node dimension name and length
  mesh2D.nodeDim = this->GetVarDimName(mesh2D.nodeCoordXVar, 0);
  mesh2D.numNodes = this->GetDimLen(mesh2D.nodeDim);

  //
  // Get face-node connectivity variable and its dimensions
  //

  if (not this->VarHasAtt(mesh2D.meshTopologyVar, "face_node_connectivity"))
  {
    std::cerr << "GetMesh2DDescription: UGRID topology variable must have face_node_connectivity attribute.\n";
    return UGRIDMeshDescription();
  }

  mesh2D.faceNodeConnVar = this->GetAttText(mesh2D.meshTopologyVar,
                           "face_node_connectivity");

  // Number of faces must be the first dimension, although
  // the UGRID conventions are more lenient
  mesh2D.faceDim = this->GetVarDimName(mesh2D.faceNodeConnVar, 0);

  if (this->VarHasAtt(mesh2D.meshTopologyVar, "face_dimension"))
  {
    if (mesh2D.faceDim !=
        this->GetAttText(mesh2D.meshTopologyVar, "face_dimension"))
    {
      std::cerr << "GetMesh2DDescription: face_node_connectivity must have face_dimension first.\n";
      return UGRIDMeshDescription();
    }
  }

  mesh2D.numFaces = this->GetDimLen(mesh2D.faceDim);

  // Assume that number of vertices per face is the second dimension
  mesh2D.vertDim = this->GetVarDimName(mesh2D.faceNodeConnVar, 1);
  mesh2D.numVertsPerFace = this->GetDimLen(mesh2D.vertDim);

  // Correction for non-zero start index
  if (this->VarHasAtt(mesh2D.faceNodeConnVar, "start_index"))
  {
    mesh2D.faceNodeStartIdx = static_cast<long long>(
      this->GetAttInt(mesh2D.faceNodeConnVar, "start_index"));
  }

  return mesh2D;
}
