#ifndef netCDFLFRicReaderUtils_h
#define netCDFLFRicReaderUtils_h

#include <vtkMath.h>
#include <vtkAssume.h>
#include <vtkDataArrayAccessor.h>
#include <vector>

// MSVC compiler requires explicit symbol import/export for DLLs
// Use mechanism provided by VTK/ParaView build system to handle
// this automatically
#include "vtkNetCDFLFRicReaderModule.h"
#define LFRICREADERUTILS_EXPORT VTKNETCDFLFRICREADER_EXPORT

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

LFRICREADERUTILS_EXPORT void resolveLongitudeGap(std::vector<double> & nodeCoordsLon,
                                                 const double lonMin,
                                                 const double lonMax,
                                                 const double lonGapSizeThreshold);

LFRICREADERUTILS_EXPORT int computeSolidAngle(const std::vector<double> & nodeCoordsLon,
                                              const std::vector<double> & nodeCoordsLat,
                                              const std::vector<long long> & faceNodeConnectivity,
                                              const size_t numFaces,
                                              const size_t numVertsPerFace,
                                              double & solidAngle,
                                              bool & hasWrapAroundCells);

LFRICREADERUTILS_EXPORT int computeEdgeLength(const std::vector<double> & nodeCoordsX,
                                              const std::vector<double> & nodeCoordsY,
                                              const std::vector<long long> & faceNodeConnectivity,
                                              const size_t numFaces,
                                              const size_t numVertsPerFace,
                                              double & edgeLengthX,
                                              double & edgeLengthY);

LFRICREADERUTILS_EXPORT void resolvePeriodicGrid(std::vector<double> & nodeCoordsX,
                                                 std::vector<double> & nodeCoordsY,
                                                 std::vector<long long> & faceNodeConnectivity,
                                                 const size_t numFaces,
                                                 const size_t numVertsPerFace,
                                                 const bool globalModel,
                                                 const double latMin,
                                                 const double latMax,
                                                 const double lonMin,
                                                 const double lonMax);

void prepareGrid(std::vector<double> & nodeCoordsX,
                 std::vector<double> & nodeCoordsY,
                 std::vector<long long> & faceNodeConnectivity,
                 const size_t numFaces,
                 const size_t numVertsPerFace,
                 const bool isPlanarLAM);
#endif
