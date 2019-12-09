#include "netCDFLFRicReaderUtils.h"

#include <unordered_map>
#include <algorithm>

void resolvePeriodicGrid(std::vector<double> & nodeCoordsX,
                         std::vector<double> & nodeCoordsY,
                         std::vector<long long> & faceNodeConnectivity,
                         const size_t numFaces,
                         const size_t numVertsPerFace)
{
  vtkDebugMacro("Entering resolvePeriodicGrid..." << endl);

  //
  // Detect grid type
  //

  // Work out horizontal grid ranges to distinguish cubed-sphere and biperiodic grids
  // This could be based on netCDF variable attributes if available
  const double xMin = *std::min_element(nodeCoordsX.begin(), nodeCoordsX.end());
  const double xMax = *std::max_element(nodeCoordsX.begin(), nodeCoordsX.end());
  const double yMin = *std::min_element(nodeCoordsY.begin(), nodeCoordsY.end());
  const double yMax = *std::max_element(nodeCoordsY.begin(), nodeCoordsY.end());

  // These offsets are needed to mirror nodes into the correct location
  double xOffset;
  double yOffset;

  // We expect lon between 0..360 and lat between -90..90 for a cubed-sphere grid
  if (xMin >= 0.0 and xMax <= 360.0 and yMin >= -90.0 and yMax <= 90.0)
  {
    xOffset = 180.0;
    yOffset = 0.0;
    vtkDebugMacro("Detected cubed-sphere grid, setting xOffset=" <<
                  xOffset << " yOffset=" << yOffset << endl);
  }
  else
  {
    xOffset = 0.5*(xMin + xMax);
    yOffset = 0.5*(yMin + yMax);
    vtkDebugMacro("Detected grid other than cubed-sphere, setting xOffset=" <<
                  xOffset << " yOffset=" << yOffset << endl);
  }

  // Compute xy grid size
  const double gridDx = xMax-xMin;
  const double gridDy = yMax-yMin;

  // Need to keep track of nodes that have already been duplicated,
  // to avoid degenerate nodes
  std::unordered_map<long long, long long> mirrorNodesX;
  std::unordered_map<long long, long long> mirrorNodesY;
  std::unordered_map<long long, long long> mirrorNodesXY;
  std::unordered_map<long long, long long>::const_iterator mirrorNodesIt;

  // Search all faces in 2D grid
  for (size_t iFace = 0; iFace < numFaces; iFace++)
  {
    // Face base index in connectivity array
    const size_t faceBaseIdx = iFace*numVertsPerFace;

    // Compute xy face size
    long long nodeId = faceNodeConnectivity[faceBaseIdx];
    double xMinFace = nodeCoordsX[nodeId];
    double xMaxFace = xMinFace;
    double yMinFace = nodeCoordsY[nodeId];
    double yMaxFace = yMinFace;
    for (size_t iVertex = 1; iVertex < numVertsPerFace; iVertex++)
    {
      nodeId = faceNodeConnectivity[faceBaseIdx + iVertex];
      const double xVertex = nodeCoordsX[nodeId];
      const double yVertex = nodeCoordsY[nodeId];
      if (xVertex < xMinFace) xMinFace = xVertex;
      if (xVertex > xMaxFace) xMaxFace = xVertex;
      if (yVertex < yMinFace) yMinFace = yVertex;
      if (yVertex > yMaxFace) yMaxFace = yVertex;
    }
    const double faceDx = xMaxFace - xMinFace;
    const double faceDy = yMaxFace - yMinFace;

    // Find faces that span across the grid
    const bool spanX = faceDx > 0.5*gridDx;
    const bool spanY = faceDy > 0.5*gridDy;

    // If face spans across, loop over vertices (nodes) and mirror
    // nodes at the left boundary to resolve grid periodicity
    if (spanX or spanY)
    {
      for (size_t iVertex = 0; iVertex < numVertsPerFace; iVertex++)
      {
	nodeId = faceNodeConnectivity[faceBaseIdx + iVertex];

        // Offset node coordinates to compute their mirror locations
        const double nodeCoordsOffsetX = nodeCoordsX[nodeId] - xOffset;
        const double nodeCoordsOffsetY = nodeCoordsY[nodeId] - yOffset;

        // Mirror corner nodes
        if (spanX and spanY and nodeCoordsOffsetX < 0 and nodeCoordsOffsetY < 0)
        {
          // Keep track of mirrored node to avoid degeneracy; insert a new node if
          // no mirror node has been created yet
          mirrorNodesIt = mirrorNodesXY.find(nodeId);
          if (mirrorNodesIt == mirrorNodesXY.end())
          {
            // Add node and register it in face-node connectivity array
            nodeCoordsX.push_back(-nodeCoordsOffsetX + xOffset);
            nodeCoordsY.push_back(-nodeCoordsOffsetY + yOffset);
            const long long newNodeId = nodeCoordsX.size() - 1;
            mirrorNodesXY.insert({nodeId, newNodeId});
            faceNodeConnectivity[faceBaseIdx + iVertex] = newNodeId;
          }
          else
          {
            // Use previously added node
            faceNodeConnectivity[faceBaseIdx + iVertex] = mirrorNodesIt->second;
          }
        }
        // Mirror nodes on left domain boundary
        else if (spanX && nodeCoordsOffsetX < 0)
        {
          mirrorNodesIt = mirrorNodesX.find(nodeId);
          if (mirrorNodesIt == mirrorNodesX.end())
          {
            nodeCoordsX.push_back(-nodeCoordsOffsetX + xOffset);
            nodeCoordsY.push_back( nodeCoordsOffsetY + yOffset);
            const long long newNodeId = nodeCoordsX.size() - 1;
            mirrorNodesX.insert({nodeId, newNodeId});
            faceNodeConnectivity[faceBaseIdx + iVertex] = newNodeId;
          }
          else {
            faceNodeConnectivity[faceBaseIdx + iVertex] = mirrorNodesIt->second;
          }
        }
        // Mirror nodes on bottom domain boundary
        else if (spanY && nodeCoordsOffsetY < 0)
        {
          mirrorNodesIt = mirrorNodesY.find(nodeId);
          if (mirrorNodesIt == mirrorNodesY.end())
          {
            nodeCoordsX.push_back( nodeCoordsOffsetX + xOffset);
            nodeCoordsY.push_back(-nodeCoordsOffsetY + yOffset);
            const long long newNodeId = nodeCoordsX.size() - 1;
            mirrorNodesY.insert({nodeId, newNodeId});
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
 
  vtkDebugMacro("mirrorNodesX: " << mirrorNodesX.size() << endl);
  vtkDebugMacro("mirrorNodesY: " << mirrorNodesY.size() << endl);
  vtkDebugMacro("mirrorNodesXY: " << mirrorNodesXY.size() << endl);

  vtkDebugMacro("Finished resolvePeriodicGrid" << endl);
}
