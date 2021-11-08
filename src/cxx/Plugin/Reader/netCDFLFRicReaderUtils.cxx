#include "netCDFLFRicReaderUtils.h"

#include <unordered_map>
#include <algorithm>

#define errorMacro(x) std::cerr << x
#ifdef DEBUG
#define debugMacro(x) std::cerr << x
#else
#define debugMacro(x)
#endif

// Some local area models cross the dateline with wrap-around cells at +-180 degree longitude, but
// do not cover the entire 360 degree longitude range. Mesh node longitudes need to be shifted by
// 360 degrees in such cases, to avoid a fractured mesh. This routine looks for gaps in mesh node
// longitudes that are larger than a given threshold and applies the shift to the western mesh fraction.
void resolveLongitudeGap(std::vector<double> & nodeCoordsLon, const double lonMin, const double lonMax,
                         const double lonGapSizeThreshold)
{
  // Generate a histogram of point longitudes with 1 degree worst-case resolution
  const size_t numLonBins = 360;
  const double lonBinWidth = (lonMax-lonMin)/static_cast<double>(numLonBins);
  std::vector<long long> lonHist(numLonBins);
  const size_t numNodeCoords = nodeCoordsLon.size();
  for (size_t idx = 0; idx < numNodeCoords; idx++)
  {
    long binId = static_cast<long>(std::floor((nodeCoordsLon[idx]-lonMin)/lonBinWidth));
    binId = std::min(std::max(binId,0l), static_cast<long>(numLonBins-1));
    lonHist[static_cast<size_t>(binId)]++;
  }

  // Find the largest gap (consecutive number of empty bins) in the histogram
  size_t zeroCounter = 0;
  size_t zeroCounterMax = 0;
  size_t lastGapIdx = 0;
  for (size_t idx = 0; idx < numLonBins; idx++)
  {
    if (lonHist[idx] == 0)
    {
      zeroCounter++;
    }
    else
    {
      zeroCounter = 0;
    }
    if (zeroCounter > zeroCounterMax)
    {
      zeroCounterMax = zeroCounter;
      lastGapIdx = idx;
    }
  }
  const double lonGapSize = static_cast<double>(zeroCounterMax)*lonBinWidth;

  debugMacro("resolveLongitudeGap: Largest gap in longitudes is " << zeroCounterMax <<
             " bins, which corresponds to " << lonGapSize << " degrees" << endl);

  // Assume that a longitude gap that is larger than this threshold means that
  // the LAM crosses the dateline, in which case its western part will be
  // shifted by 360 degrees
  if (lonGapSize > lonGapSizeThreshold)
  {
    const double lonShiftThreshold = lonMin+static_cast<double>(lastGapIdx)*lonBinWidth;
    for (size_t idx = 0; idx < numNodeCoords; idx++)
    {
      if (nodeCoordsLon[idx] < lonShiftThreshold)
      {
        nodeCoordsLon[idx] += 360.0;
      }
    }
    debugMacro("resolveLongitudeGap: Shifted longitudes below " << lonShiftThreshold <<
               " by 360 degrees eastward" << endl);
  }
  else
  {
    debugMacro("resolveLongitudeGap: Gap is below threshold " << lonGapSizeThreshold <<
               ", assuming that the mesh is not fractured" << endl);
  }
}

// Computing total solid angle covered by the mesh, excluding wrap-around cells
int computeSolidAngle(const std::vector<double> & nodeCoordsLon,
                      const std::vector<double> & nodeCoordsLat,
                      const std::vector<long long> & faceNodeConnectivity,
                      const size_t numFaces,
                      const size_t numVertsPerFace,
                      double & solidAngle,
                      bool & hasWrapAroundCells)
{
  solidAngle = 0.0;
  hasWrapAroundCells = false;

  if (numVertsPerFace != 4)
  {
    errorMacro("computeSolidAngle: ERROR: Only quadrilateral cells are currently " <<
               "supported in the horizontal domain." << endl);
    return 1;
  }

  // Traverse grid and compute total solid angle of the horizontal domain on the sphere
  double lon[4];
  double lat[4];
  const double deg2rad = vtkMath::Pi()/180.0;
  for (size_t iFace = 0; iFace < numFaces; iFace++)
  {
    // Face base index in connectivity array
    const size_t faceBaseIdx = iFace*numVertsPerFace;

    for (size_t iVertex = 0; iVertex < numVertsPerFace; iVertex++)
    {
      const long long nodeId = faceNodeConnectivity[faceBaseIdx + iVertex];
      lon[iVertex] = nodeCoordsLon[nodeId]*deg2rad;
      lat[iVertex] = nodeCoordsLat[nodeId]*deg2rad;
    }

    // Compute solid angle of cell on the sphere using the cross-product of cell
    // diagonal vectors using coordinates mu=sin(lat), nu=lon/pi
    const double cellSolidAngle = 0.5*((lon[2]-lon[0])*(sin(lat[3])-sin(lat[1])) - \
                                       (lon[3]-lon[1])*(sin(lat[2])-sin(lat[0])));

    // Cells that wrap around at a periodic boundary will have negative solid angle;
    // ignore these in the calculation
    if (cellSolidAngle > 0.0)
    {
      solidAngle += cellSolidAngle;
    }
    else
    {
      hasWrapAroundCells = true;
    }
  }

  debugMacro("computeSolidAngle: Mesh covers a total solid angle of " << solidAngle << endl);
  if (hasWrapAroundCells)
  {
    debugMacro("computeSolidAngle: Mesh covers a total solid angle of " << solidAngle << endl);
  }
  else
  {
    debugMacro("computeSolidAngle: Mesh covers a total solid angle of " << solidAngle << endl);
  }
  return 0;
}

//----------------------------------------------------------------------------

void resolvePeriodicGrid(std::vector<double> & nodeCoordsLon,
                         std::vector<double> & nodeCoordsLat,
                         std::vector<long long> & faceNodeConnectivity,
                         const size_t numFaces,
                         const size_t numVertsPerFace,
                         const bool globalModel,
                         const double latMin,
                         const double latMax,
                         const double lonMin,
                         const double lonMax)
{
  debugMacro("Entering resolvePeriodicGrid..." << endl);

  double lonOffset;
  double latOffset;
  double cellThreshold;
  bool mirrorFromWest;

  if (globalModel)
  {
    // Assume cubed-sphere grid with range (newer output files use -180..180 degrees)
    if (lonMin >= -180.0 and lonMax <= 180.0)
    {
      lonOffset = 0.0;
      latOffset = 0.0;
      cellThreshold = 0.6;
      mirrorFromWest = false;
      debugMacro("resolvePeriodicGrid: Assuming cubed-sphere grid with -180..180 lon range, setting lonOffset=" <<
                 lonOffset << " latOffset=" << latOffset << endl);
    }
    else if (lonMin >= 0.0 and lonMax <= 360.0)
    {
      lonOffset = 180.0;
      latOffset = 0.0;
      cellThreshold = 0.5;
      mirrorFromWest = true;
      debugMacro("resolvePeriodicGrid: Assuming cubed-sphere grid with 0..360 lon range, setting lonOffset=" <<
                 lonOffset << " latOffset=" << latOffset << endl);
    }
    else
    {
      errorMacro("resolvePeriodicGrid: ERROR: Unexpected longitude range " << lonMin <<
                 ".." << lonMax << endl);
      return;
    }
  }
  else
  {
    lonOffset = 0.5*(lonMin + lonMax);
    latOffset = 0.5*(latMin + latMax);
    cellThreshold = 0.5;
    mirrorFromWest = true;
    debugMacro("Detected LAM grid, setting lonOffset=" <<
               lonOffset << " latOffset=" << latOffset << endl);
  }

  // Compute lon-lat grid size
  const double gridDlon = lonMax-lonMin;
  const double gridDlat = latMax-latMin;

  // Need to keep track of nodes that have already been duplicated,
  // to avoid degenerate nodes
  std::unordered_map<long long, long long> mirrorNodesLon;
  std::unordered_map<long long, long long> mirrorNodesLat;
  std::unordered_map<long long, long long> mirrorNodesLonLat;
  std::unordered_map<long long, long long>::const_iterator mirrorNodesIt;

  // Search all faces in 2D grid
  for (size_t iFace = 0; iFace < numFaces; iFace++)
  {
    // Face base index in connectivity array
    const size_t faceBaseIdx = iFace*numVertsPerFace;

    // Compute xy face size
    long long nodeId = faceNodeConnectivity[faceBaseIdx];
    double lonMinFace = nodeCoordsLon[nodeId];
    double lonMaxFace = lonMinFace;
    double latMinFace = nodeCoordsLat[nodeId];
    double latMaxFace = latMinFace;
    for (size_t iVertex = 1; iVertex < numVertsPerFace; iVertex++)
    {
      nodeId = faceNodeConnectivity[faceBaseIdx + iVertex];
      const double lonVertex = nodeCoordsLon[nodeId];
      const double latVertex = nodeCoordsLat[nodeId];
      if (lonVertex < lonMinFace) lonMinFace = lonVertex;
      if (lonVertex > lonMaxFace) lonMaxFace = lonVertex;
      if (latVertex < latMinFace) latMinFace = latVertex;
      if (latVertex > latMaxFace) latMaxFace = latVertex;
    }
    const double faceDlon = lonMaxFace - lonMinFace;
    const double faceDlat = latMaxFace - latMinFace;

    // Find faces that span across the grid, as well as certain "defective" cells
    const bool spanLon = faceDlon > cellThreshold*gridDlon;
    const bool spanLat = faceDlat > cellThreshold*gridDlat;

    // If face spans across, loop over vertices (nodes) and mirror
    // nodes at the left boundary to resolve grid periodicity
    if (spanLon or spanLat)
    {
      for (size_t iVertex = 0; iVertex < numVertsPerFace; iVertex++)
      {
        nodeId = faceNodeConnectivity[faceBaseIdx + iVertex];

        // Offset node coordinates to compute their mirror locations
        const double nodeCoordsOffsetLon = nodeCoordsLon[nodeId] - lonOffset;
        const double nodeCoordsOffsetLat = nodeCoordsLat[nodeId] - latOffset;

        // Mirror corner nodes
        if (spanLon and spanLat and nodeCoordsOffsetLon < 0 and nodeCoordsOffsetLat < 0)
        {
          // Keep track of mirrored node to avoid degeneracy; insert a new node if
          // no mirror node has been created yet
          mirrorNodesIt = mirrorNodesLonLat.find(nodeId);
          if (mirrorNodesIt == mirrorNodesLonLat.end())
          {
            // Add node and register it in face-node connectivity array
            nodeCoordsLon.push_back(-nodeCoordsOffsetLon + lonOffset);
            nodeCoordsLat.push_back(-nodeCoordsOffsetLat + latOffset);
            const long long newNodeId = nodeCoordsLon.size() - 1;
            mirrorNodesLonLat.insert({nodeId, newNodeId});
            faceNodeConnectivity[faceBaseIdx + iVertex] = newNodeId;
          }
          else
          {
            // Use previously added node
            faceNodeConnectivity[faceBaseIdx + iVertex] = mirrorNodesIt->second;
          }
        }
        // Mirror nodes on left or right domain boundary
        else if (spanLon and ((nodeCoordsOffsetLon < 0 and mirrorFromWest) or \
                              (nodeCoordsOffsetLon > 0 and not mirrorFromWest)))
        {
          mirrorNodesIt = mirrorNodesLon.find(nodeId);
          if (mirrorNodesIt == mirrorNodesLon.end())
          {
            nodeCoordsLon.push_back(-nodeCoordsOffsetLon + lonOffset);
            nodeCoordsLat.push_back( nodeCoordsOffsetLat + latOffset);
            const long long newNodeId = nodeCoordsLon.size() - 1;
            mirrorNodesLon.insert({nodeId, newNodeId});
            faceNodeConnectivity[faceBaseIdx + iVertex] = newNodeId;
          }
          else
          {
            faceNodeConnectivity[faceBaseIdx + iVertex] = mirrorNodesIt->second;
          }
        }
        // Mirror nodes on bottom domain boundary
        else if (spanLat and nodeCoordsOffsetLat < 0)
        {
          mirrorNodesIt = mirrorNodesLat.find(nodeId);
          if (mirrorNodesIt == mirrorNodesLat.end())
          {
            nodeCoordsLon.push_back( nodeCoordsOffsetLon + lonOffset);
            nodeCoordsLat.push_back(-nodeCoordsOffsetLat + latOffset);
            const long long newNodeId = nodeCoordsLon.size() - 1;
            mirrorNodesLat.insert({nodeId, newNodeId});
            faceNodeConnectivity[faceBaseIdx + iVertex] = newNodeId;
          }
          else
          {
            faceNodeConnectivity[faceBaseIdx + iVertex] = mirrorNodesIt->second;
          }
        }
      }
    }
  }

  debugMacro("resolvePeriodicGrid: mirrorNodesLon: " << mirrorNodesLon.size() << endl);
  debugMacro("resolvePeriodicGrid: mirrorNodesLat: " << mirrorNodesLat.size() << endl);
  debugMacro("resolvePeriodicGrid: mirrorNodesLonLat: " << mirrorNodesLonLat.size() << endl);

  debugMacro("Finished resolvePeriodicGrid" << endl);
}

//----------------------------------------------------------------------------

void prepareGrid(std::vector<double> & nodeCoordsX,
                 std::vector<double> & nodeCoordsY,
                 std::vector<long long> & faceNodeConnectivity,
                 const size_t numFaces,
                 const size_t numVertsPerFace,
                 const bool isPlanarLAM)
{
  debugMacro("Entering prepareGrid..." << endl);

  //
  // Determine mesh dimensions
  //

  // Determine horizontal coordinate ranges
  const auto xMinMaxIt = std::minmax_element(nodeCoordsX.begin(), nodeCoordsX.end());
  const double xMin = *xMinMaxIt.first;
  const double xMax = *xMinMaxIt.second;

  const auto yMinMaxIt = std::minmax_element(nodeCoordsY.begin(), nodeCoordsY.end());
  const double yMin = *yMinMaxIt.first;
  const double yMax = *yMinMaxIt.second;

  debugMacro("prepareGrid: Found X coord range of [" <<
             xMin << "," << xMax << "] and Y coord range of [" <<
             yMin << "," << yMax << "]" << endl);

  //
  // Detect and possibly resolve gap in longitudes for non-planar-LAMs that cross the dateline
  //

  if (not isPlanarLAM)
  {
    const double lonGapSizeThreshold = 30.0;
    resolveLongitudeGap(nodeCoordsX, xMin, xMax, lonGapSizeThreshold);
  }

  //
  // Detect grid type (global or LAM), if not planar LAM
  //

  double solidAngle = 0.0;
  bool hasWrapAroundCells = false;
  if (not isPlanarLAM)
  {
    if (computeSolidAngle(nodeCoordsX, nodeCoordsY, faceNodeConnectivity,
                          numFaces, numVertsPerFace, solidAngle, hasWrapAroundCells))
    {
      errorMacro("prepareGrid: ERROR: computeSolidAngle returned error" << endl);
      return;
    }
  }

  const double solidAngleThreshold = 3.0*vtkMath::Pi();
  bool globalModel = false;
  if (solidAngle > solidAngleThreshold and not isPlanarLAM)
  {
    globalModel = true;
    debugMacro("prepareGrid: Assuming mesh is global." << endl);
  }
  else
  {
    debugMacro("prepareGrid: Assuming mesh is LAM." << endl);
  }

  // computeSolidAngle scans for wrap-around cells, but is only executed when isPlanarLAM == false
  if (hasWrapAroundCells or isPlanarLAM)
  {
    resolvePeriodicGrid(nodeCoordsX, nodeCoordsY, faceNodeConnectivity,
                        numFaces, numVertsPerFace, globalModel, yMin,
                        yMax, xMin, xMax);
  }

  debugMacro("Finished prepareGrid" << endl);
}
