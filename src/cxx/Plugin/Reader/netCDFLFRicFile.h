/**
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
  size_t numEdgesPerFace;

  // NetCDF dimension IDs
  int nodeDimId;
  int edgeDimId;
  int faceDimId;
  int vertsPerFaceDimId;
  int edgesPerFaceDimId;

  // Alternative face dimension ID for old LFRic output files
  int faceDimIdAlt;

  // NetCDF variable IDs
  int meshTopologyVarId;
  int nodeCoordXVarId;
  int nodeCoordYVarId;
  int faceNodeConnVarId;
  int faceEdgeConnVarId;
  int edgeCoordXVarId;
  int edgeCoordYVarId;

  // Start indices for face-node and edge-node connectivity
  long long faceNodeStartIdx;
  long long faceEdgeStartIdx;

  UGRIDMeshDescription() : isLFRicXIOSFile(false), numTopologies(0),
    meshType(unknownMesh), numNodes(0), numEdges(0), numFaces(0),
    numVertsPerFace(0), numEdgesPerFace(0),
    // Initialise netCDF IDs to -1 - valid IDs are always >= 0
    nodeDimId(-1), edgeDimId(-1), faceDimId(-1), faceDimIdAlt(-1),
    vertsPerFaceDimId(-1), edgesPerFaceDimId(-1),
    meshTopologyVarId(-1), nodeCoordXVarId(-1), nodeCoordYVarId(-1),
    edgeCoordXVarId(-1), edgeCoordYVarId(-1), faceNodeConnVarId(-1), 
    faceNodeStartIdx(0), faceEdgeConnVarId(-1), faceEdgeStartIdx(0) {}
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

  /**
   * Construct new netCDFLFRicFile object for netCDF file fileName.
   * The file is opened on construction and will be closed by the
   * destructor.
   */
  netCDFLFRicFile(const char* fileName);
  ~netCDFLFRicFile();

  /**
   * Check if the class constructor managed to open the netCDF file.
   */
  bool IsFileOpen();

  /**
   * Return name of the netCDF file.
   */
  const char* GetFileName();

  /**
   * Check if netCDF dimension dimName exists.
   */
  bool HasDim(const std::string& dimName);

  /**
   * Return length of netCDF dimension with ID dimId.
   */
  size_t GetDimLen(const int dimId);

  /**
   * Return netCDF dimension ID for dimension dimName.
   */
  int GetDimId(const std::string& dimName);

  /**
   * Return netCDF dimension name for ID dimId.
   */
  std::string GetDimName(const int dimId);

  /**
   * Return netCDF variable ID for variable varName.
   */
  int GetVarId(const std::string& varName);

  /**
   * Return netCDF variable name for ID varID.
   */
  std::string GetVarName(const int varId);

  /**
   * Return the number of dimensions for netCDF variable with ID varId.
   */
  size_t GetVarNumDims(const int varId);

  /**
   * Return ID of the dimth dimension of netCDF variable with ID varId.
   */
  int GetVarDimId(const int varId, const size_t dim);

  /**
   * Return integer attribute attName of netCDF variable with ID varId.
   */
  int GetAttInt(const int varId, const std::string& attName);

  ///@{
  /**
   * Return string attribute attName of netCDF variable with ID varId
   * or name varName.
   */
  std::string GetAttText(const int varId, const std::string& attName);
  std::string GetAttText(const std::string& varName, const std::string& attName);
  ///@}

  /**
   * Return string attribute attName of netCDF variable with ID varId
   * split into a vector of strings.
   */
  std::vector<std::string> GetAttTextSplit(const int varId,
                                           const std::string& attName);

  /**
   * Check if netCDF variable varName exists.
   */
  bool HasVar(const std::string& varName);

  /**
   * Check if netCDF variable with ID varId has attribute attName.
   */
  bool VarHasAtt(const int varId, const std::string& attName);

  /**
   * Return the total number of netCDF variables in the file.
   */
  size_t GetNumVars();

  /**
   * Load data for double precision netCDF variable with ID varId
   * into user-provided memory buffer, using vectors start and count
   * for subsetting.
   *
   * @note Make sure that buffer is sufficiently larget to hold the
   * data!
   */
  void LoadVarDouble(const int varId,
                     const std::vector<size_t>& start,
                     const std::vector<size_t>& count,
                     double* buffer);

  ///@{
  /**
   * Return double precision or long long vectors with data for
   * netCDF variable with ID varID, using vectors start and count
   * for subsetting.
   */
  std::vector<double> GetVarDouble(const int varId,
                                   const std::vector<size_t>& start,
                                   const std::vector<size_t>& count);

  std::vector<long long> GetVarLongLong(const int varId,
                                        const std::vector<size_t>& start,
                                        const std::vector<size_t>& count);
  ///@}

  UGRIDMeshDescription GetMesh2DDescription();

  std::map<std::string, CFAxis> GetZAxisDescription(const bool isLFRicXIOSFile,
                                                    const mesh2DTypes meshType);

  CFAxis GetTAxisDescription();

  void UpdateFieldMaps(const UGRIDMeshDescription & mesh2D,
                       const std::map<std::string, CFAxis> & zAxes,
                       const CFAxis & tAxis,
                       std::map<std::string, DataField> & CellFields,
                       std::map<std::string, DataField> & PointFields);

private:

  netCDFLFRicFile(const netCDFLFRicFile&) = delete;
  void operator=(const netCDFLFRicFile&) = delete;

  int ncId;
  bool fileOpen;
  std::string FileName;
};

#endif
