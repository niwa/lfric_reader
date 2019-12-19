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
#include <map>

// LFRic 2D mesh types
enum mesh2DTypes {unknownMesh, halfLevelFaceMesh, fullLevelFaceMesh};

// Dimension types for identifying field variable dimensions
enum dimTypes {unknownAxisDim, horizontalAxisDim, verticalAxisDim,
               timeAxisDim, componentAxisDim};

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
    meshType(unknownMesh), numNodes(0), numFaces(0), numVertsPerFace(0),
    nodeDim("None"), faceDim("None"), vertDim("None"),
    meshTopologyVar("None"), nodeCoordXVar("None"),
    nodeCoordYVar("None"), faceNodeConnVar("None"),
    faceNodeStartIdx(0) {}
};

// Holds metadata for vertical axis in VTK grid
struct CFVerticalAxis
{
  // Axis length (number of cells in the vertical)
  size_t numLevels;

  // NetCDF dimension and variable name
  std::string axisDim;
  std::string axisVar;

  CFVerticalAxis() : numLevels(0), axisDim("None"), axisVar("None") {}
};

// Holds metadata for time axis
struct CFTimeAxis
{
  size_t numTimeSteps;

  // NetCDF dimension and variable name
  std::string axisDim;
  std::string axisVar;

  CFTimeAxis() : numTimeSteps(0), axisDim("None"), axisVar("None") {}
};

// Holds metadata for field variable dimensions
struct fieldDim
{
  dimTypes dimType;
  size_t dimLength;
  size_t dimStride;

  fieldDim() : dimType(unknownAxisDim), dimLength(0), dimStride(1) {}
};

// Holds metadata for data fields
struct DataField
{
  bool active;
  mesh2DTypes meshType;
  bool hasComponentDim;
  bool hasHorizontalDim;
  bool hasVerticalDim;
  bool hasTimeDim;
  std::vector<fieldDim> dims;

  DataField() : active(false), meshType(unknownMesh),
    hasComponentDim(false), hasHorizontalDim(false),
    hasVerticalDim(false), hasTimeDim(false) {}
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
                                   const std::vector<size_t>& start,
                                   const std::vector<size_t>& count);

  std::vector<long long> GetVarLongLong(const std::string& varName,
                                        const std::vector<size_t>& start,
                                        const std::vector<size_t>& count);

  UGRIDMeshDescription GetMesh2DDescription();

  CFVerticalAxis GetZAxisDescription(const bool isLFRicXIOSFile,
                                     const mesh2DTypes meshType);

  CFTimeAxis GetTAxisDescription();

  void UpdateFieldMap(std::map<std::string, DataField> & fields,
                      const std::string & fieldLoc,
                      const std::string & horizontalDim,
                      const mesh2DTypes & horizontalMeshType,
                      const std::string & verticalDim,
                      const std::string & timeDim);

private:

  netCDFLFRicFile(const netCDFLFRicFile&) = delete;
  void operator=(const netCDFLFRicFile&) = delete;

  int ncId;
  bool fileOpen;
  std::string FileName;
};

#endif
