/*
 * @class   netCDFLFRicFile
 * @brief   Utility class for netCDF file handling
 *
 * Reads metadata and data from a netCDF file using
 * convenient utility functions.
*/

#ifndef netCDFLFRicFile_h
#define netCDFLFRicFile_h
#include <vector>
#include <string>

// The following LFRic mesh 2D types are currently supported
enum mesh2DTypes {halfLevelFace, fullLevelFace};

// Holds metadata for unstructured part of VTK grid
struct UGRIDMeshDescription
{
  // LFRic XIOS output files require special treatment
  bool isLFRicXIOSFile;

  // Up to three topologies can exist
  size_t numTopologies;

  // Type of main mesh
  mesh2DTypes meshType;

  // Mesh dimensions
  size_t numNodes;
  size_t numFaces;
  size_t numVertsPerFace;

  // NetCDF dimension names
  std::string nodeDim;
  std::string faceDim;
  std::string vertDim;

  // NetCDF variable names
  std::string meshTopologyVar;
  std::string nodeCoordXVar;
  std::string nodeCoordYVar;
  std::string faceNodeConnVar;

  // Start index for face-node connectivity
  long long faceNodeStartIdx;

  UGRIDMeshDescription() : isLFRicXIOSFile(false), numTopologies(0),
    numNodes(0), numFaces(0), numVertsPerFace(0), faceNodeStartIdx(0) {}
};

// Holds metadata for vertical axis in VTK grid
struct CFVerticalAxis
{
  // Axis length (number of cells in the vertical)
  size_t numLevels;

  // NetCDF dimension and variable name
  std::string axisDim;
  std::string axisVar;

  CFVerticalAxis(): numLevels(0) {}
};

class netCDFLFRicFile
{

public:

  netCDFLFRicFile(const char* fileName);
  ~netCDFLFRicFile();

  bool IsFileOpen();

  const char* GetFileName();

  bool HasDim(const std::string& dimName);

  size_t GetDimLen(const std::string& dimName);

  size_t GetVarNumDims(const std::string& varName);

  std::string GetVarDimName(const std::string& varName, const size_t dim);

  int GetAttInt(const std::string& varName, const std::string& attName);

  std::string GetAttText(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetAttTextSplit(const std::string& varName,
                                           const std::string& attName);

  bool HasVar(const std::string& varName);

  bool VarHasDim(const std::string& varName, const std::string& dimName);

  bool VarHasAtt(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetVarNames();

  std::vector<double> GetVarDouble(const std::string& varName,
                                   const std::vector<size_t> start,
                                   const std::vector<size_t> count);

  std::vector<long long> GetVarLongLong(const std::string& varName,
                                        const std::vector<size_t> start,
                                        const std::vector<size_t> count);

  UGRIDMeshDescription GetMesh2DDescription();

  CFVerticalAxis GetZAxisDescription(const bool isLFRicXIOSFile,
                                     const mesh2DTypes meshType);

private:

  netCDFLFRicFile(const netCDFLFRicFile&) = delete;
  void operator=(const netCDFLFRicFile&) = delete;

  int ncId;
  bool fileOpen;
  std::string FileName;
};

#endif
