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
enum mesh2DTypes {unknownMesh, halfLevelEdgeMesh, halfLevelFaceMesh, fullLevelFaceMesh};

// Dimension types for identifying field variable dimensions
enum dimTypes {unknownAxisDim, horizontalAxisDim, verticalAxisDim,
               timeAxisDim, componentAxisDim};

// Field data location on VTK grid
enum fieldLocs {unknownFieldLoc, pointFieldLoc, cellFieldLoc};

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
  size_t numEdges;
  size_t numFaces;
  size_t numVertsPerFace;

  // NetCDF dimension IDs
  int nodeDimId;
  int edgeDimId;
  int faceDimId;
  int vertDimId;

  // NetCDF variable IDs
  int meshTopologyVarId;
  int nodeCoordXVarId;
  int nodeCoordYVarId;
  int faceNodeConnVarId;
  int edgeCoordXVarId;
  int edgeCoordYVarId;

  // Start index for face-node connectivity
  long long faceNodeStartIdx;

  UGRIDMeshDescription() : isLFRicXIOSFile(false), numTopologies(0),
    meshType(unknownMesh), numNodes(0), numEdges(0), numFaces(0), numVertsPerFace(0),
    // Initialise netCDF IDs to -1 - valid IDs are always >= 0
    nodeDimId(-1), edgeDimId(-1), faceDimId(-1), vertDimId(-1),
    meshTopologyVarId(-1), nodeCoordXVarId(-1), nodeCoordYVarId(-1),
    faceNodeConnVarId(-1), edgeCoordXVarId(-1), edgeCoordYVarId(-1),
    faceNodeStartIdx(0) {}
};

// Holds metadata for vertical and time axes
struct CFAxis
{
  size_t axisLength;

  // NetCDF dimension and variable IDs
  int axisDimId;
  int axisVarId;

  CFAxis() : axisLength(0), axisDimId(-1), axisVarId(-1) {}
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
  fieldLocs location;
  bool hasComponentDim;
  bool hasHorizontalDim;
  bool hasVerticalDim;
  bool hasTimeDim;
  std::vector<fieldDim> dims;

  DataField() : active(false), meshType(unknownMesh),
    location(unknownFieldLoc),
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

  size_t GetDimLen(const int dimId);

  int GetDimId(const std::string& dimName);

  std::string GetDimName(const int dimId);

  int GetVarId(const std::string& varName);

  std::string GetVarName(const int varId);

  size_t GetVarNumDims(const int varId);

  int GetVarDimId(const int varId, const size_t dim);

  int GetAttInt(const int varId, const std::string& attName);

  std::string GetAttText(const int varId, const std::string& attName);
  std::string GetAttText(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetAttTextSplit(const int varId,
                                           const std::string& attName);

  bool HasVar(const std::string& varName);

  bool VarHasAtt(const int varId, const std::string& attName);

  size_t GetNumVars();

  void LoadVarDouble(const int varId,
                     const std::vector<size_t>& start,
                     const std::vector<size_t>& count,
                     double* buffer);

  std::vector<double> GetVarDouble(const int varId,
                                   const std::vector<size_t>& start,
                                   const std::vector<size_t>& count);

  std::vector<long long> GetVarLongLong(const int varId,
                                        const std::vector<size_t>& start,
                                        const std::vector<size_t>& count);

  UGRIDMeshDescription GetMesh2DDescription();

  CFAxis GetZAxisDescription(const bool isLFRicXIOSFile,
                             const mesh2DTypes meshType);

  CFAxis GetTAxisDescription();

  void UpdateFieldMap(std::map<std::string, DataField> & fields,
                      const std::string & fieldLoc,
                      const int horizontalDimId,
                      const mesh2DTypes & horizontalMeshType,
                      const int verticalDimId_1,
                      const int verticalDimId_2,
                      const int timeDimId);

private:

  netCDFLFRicFile(const netCDFLFRicFile&) = delete;
  void operator=(const netCDFLFRicFile&) = delete;

  int ncId;
  bool fileOpen;
  std::string FileName;
};

#endif
