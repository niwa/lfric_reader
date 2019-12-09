#ifndef netCDFLFRicReaderUtils_h
#define netCDFLFRicReaderUtils_h

#include <vtkMath.h>
#include <vtkAssume.h>
#include <vtkDataArrayAccessor.h>

#include <vector>

// Functor for setting a vtkDataArray with point coordinates
// in layer-first ordering
struct SetPointLocationWorker
{
public:

  SetPointLocationWorker(const std::vector<double> & X,
                         const std::vector<double> & Y,
                         const std::vector<double> & Z,
                         const bool Cartesian) :
    X(X), Y(Y), Z(Z), cartCoords(Cartesian) {}

  template <typename CoordArray>
  void operator()(CoordArray * pointLocs)
  {
    VTK_ASSUME(pointLocs->GetNumberOfComponents() == 3);
    vtkDataArrayAccessor<CoordArray> p(pointLocs);

    const size_t numNodes = X.size();
    const size_t numLevels = Z.size();
    vtkIdType pointBaseIdx = 0;

    if (cartCoords)
    {
      const double deg2rad = vtkMath::Pi()/180.0;
      for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++)
      {	
	// Convert from lon-lat-rad to Cartesian coordinates
	const double sinLon = sin(X[nodeIdx]*deg2rad);
	const double cosLon = cos(X[nodeIdx]*deg2rad);
	const double sinLat = sin(Y[nodeIdx]*deg2rad);
	const double cosLat = cos(Y[nodeIdx]*deg2rad);

	const double normX = cosLon*cosLat;
	const double normY = sinLon*cosLat;
	const double normZ = sinLat;

        // Compute column-wise to avoid recomputing trigonometric functions
	for (size_t levelIdx = 0; levelIdx < numLevels; levelIdx++)
	{
          const vtkIdType pointIdx = pointBaseIdx+
            static_cast<vtkIdType>(levelIdx*numNodes);
	  p.Set(pointIdx, 0, Z[levelIdx]*normX);
	  p.Set(pointIdx, 1, Z[levelIdx]*normY);
	  p.Set(pointIdx, 2, Z[levelIdx]*normZ);
	}
	pointBaseIdx++;
      }
    }
    else
    {
      for (size_t levelIdx = 0; levelIdx < numLevels; levelIdx++)
      {
	for (size_t nodeIdx = 0; nodeIdx < numNodes; nodeIdx++)
	{
	  p.Set(pointBaseIdx, 0, X[nodeIdx]);
	  p.Set(pointBaseIdx, 1, Y[nodeIdx]);
	  p.Set(pointBaseIdx, 2, Z[levelIdx]);
	  pointBaseIdx++;
	}
      }
    }
  }

private:

  const std::vector<double> & X;
  const std::vector<double> & Y;
  const std::vector<double> & Z;
  const bool cartCoords;
};

// Functor for setting a vtkDataArray with cell connectivities
// in layer-first ordering
struct SetConnectivityWorker
{
public:

  SetConnectivityWorker(const std::vector<long long> & faceNode,
                        const size_t nLevels,
                        const size_t nFaces,
                        const size_t nVertsPerFace,
                        const size_t nNodes) :
    faceNodeConnect(faceNode), numLevels(nLevels), numFaces(nFaces),
    numVertsPerFace(nVertsPerFace), numNodes(nNodes) {}

  template <typename ConnectArray>
  void operator()(ConnectArray * cellConnect)
  {
    VTK_ASSUME(cellConnect->GetNumberOfComponents() == 1);
    vtkDataArrayAccessor<ConnectArray> c(cellConnect);

    // Cell definitions are a flat list with format
    // Number of points n, point ID 1, ..., point ID n
    vtkIdType cellIdx = 0;
    for (size_t iLevel = 0; iLevel < numLevels; iLevel++)
    {
      const long long levelOffset = static_cast<long long>(iLevel*numNodes);
      for (size_t iFace = 0; iFace < numFaces; iFace++)
      {
        c.Set(cellIdx, 0, static_cast<vtkIdType>(2*numVertsPerFace));
        const size_t faceBaseIdx = iFace*numVertsPerFace;
    	cellIdx++;
    	for (size_t iVertex = 0; iVertex < numVertsPerFace; iVertex++)
    	{
          const vtkIdType pointId = static_cast<vtkIdType>(
            faceNodeConnect[faceBaseIdx+iVertex] + levelOffset);
          c.Set(cellIdx, 0, pointId);
          c.Set(cellIdx+static_cast<vtkIdType>(numVertsPerFace), 0,
                pointId+static_cast<vtkIdType>(numNodes));
    	  cellIdx++;
    	}
        cellIdx += static_cast<vtkIdType>(numVertsPerFace);
      }
    }
  }

private:

  const std::vector<long long> & faceNodeConnect;
  const size_t numLevels;
  const size_t numFaces;
  const size_t numVertsPerFace;
  const size_t numNodes;
};

void resolvePeriodicGrid(std::vector<double> & nodeCoordsX,
                         std::vector<double> & nodeCoordsY,
                         std::vector<long long> & faceNodeConnectivity,
                         const size_t numFaces,
                         const size_t numVertsPerFace);

#endif
