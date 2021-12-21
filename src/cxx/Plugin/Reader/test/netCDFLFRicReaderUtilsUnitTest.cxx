#include "Catch2/catch.hpp"
#include "netCDFLFRicReaderUtils.h"

#include <vector>

// ------------------------------------------------------------------------------------------

TEST_CASE( "ResolveLongitudeGap Test", "[basic]" )
{
  SECTION( "No GAP Detected In LAM Without Gap" )
  {
    const std::vector<double> lonOriginal({10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0});
    const auto lonMinMaxIt = std::minmax_element(lonOriginal.begin(), lonOriginal.end());
    const double lonMin = *lonMinMaxIt.first;
    const double lonMax = *lonMinMaxIt.second;

    double lonGapSizeThreshold = 10.0;
    std::vector<double> lon = lonOriginal;
    resolveLongitudeGap(lon, lonMin, lonMax, lonGapSizeThreshold);

    // Output should be bit-identical
    REQUIRE( lon == lonOriginal );
  }

  SECTION( "No GAP Detected In Global Model Without Gap" )
  {
    const std::vector<double> lonOriginal({-180.0,-135.0,-90.0,-45.0,0.0,45.0,90.0,135.0});
    const auto lonMinMaxIt = std::minmax_element(lonOriginal.begin(), lonOriginal.end());
    const double lonMin = *lonMinMaxIt.first;
    const double lonMax = *lonMinMaxIt.second;

    double lonGapSizeThreshold = 50.0;
    std::vector<double> lon = lonOriginal;
    resolveLongitudeGap(lon, lonMin, lonMax, lonGapSizeThreshold);

    // Output should be bit-identical
    REQUIRE( lon == lonOriginal );
  }

  SECTION( "GAP Detected And Resolved In LAM With Gap" )
  {
    const std::vector<double> lonOriginal({-175.0,-170.0,-165.0,-160.0,165.0,170.0,175.0,180.0});
    const auto lonMinMaxIt = std::minmax_element(lonOriginal.begin(), lonOriginal.end());
    const double lonMin = *lonMinMaxIt.first;
    const double lonMax = *lonMinMaxIt.second;

    double lonGapSizeThreshold = 30.0;
    std::vector<double> lon = lonOriginal;
    resolveLongitudeGap(lon, lonMin, lonMax, lonGapSizeThreshold);

    for(size_t i = 0; i < lon.size(); i++)
    {
      if(lonOriginal[i] < 0.0)
      {
        REQUIRE( lon[i] == Approx(lonOriginal[i]+360.0) );
      }
      else
      {
        REQUIRE( lon[i] == Approx(lonOriginal[i]) );
      }
    }
  }
}

// ------------------------------------------------------------------------------------------

TEST_CASE( "ComputeSolidAngle Test", "[basic]" )
{
  SECTION( "Fails For Triangular Cells" )
  {
    const std::vector<double> lon({-180.0,180.0,-180.0});
    const std::vector<double> lat({ -90.0,-90.0,  90.0});
    const std::vector<long long> faceNodeConnectivity({0,1,2});
    const size_t numFaces = 1;
    const size_t numVertsPerFace = 3;

    double solidAngle = -1.0;
    bool hasWrapAroundCell = true;
    
    REQUIRE( computeSolidAngle(lon, lat, faceNodeConnectivity, numFaces,
			       numVertsPerFace, solidAngle, hasWrapAroundCell) );
    REQUIRE( solidAngle == Approx(0.0) );
    REQUIRE( !hasWrapAroundCell );
  }

  SECTION( "Correct Results for Single Global Cell" )
  {
    const std::vector<double> lon({-180.0,180.0,-180.0,180.0});
    const std::vector<double> lat({ -90.0,-90.0,  90.0, 90.0});
    const std::vector<long long> faceNodeConnectivity({0,1,3,2});
    const size_t numFaces = 1;
    const size_t numVertsPerFace = 4;

    double solidAngle = -1.0;
    bool hasWrapAroundCell = true;
    
    REQUIRE( !computeSolidAngle(lon, lat, faceNodeConnectivity, numFaces,
			           numVertsPerFace, solidAngle, hasWrapAroundCell) );
    REQUIRE( solidAngle == Approx(4.0*vtkMath::Pi()) );
    REQUIRE( !hasWrapAroundCell );
  }

  SECTION( "Correct Results for Global Model With Wrap-Around Cells" )
  {
    // Construct a "global model" with a smaller wrap-around cell
    const std::vector<double> lon({-180.0,-90.0,  0.0, 90.0,175.0,-180.0,-90.0, 0.0,90.0,175.0});
    const std::vector<double> lat({ -90.0,-90.0,-90.0,-90.0,-90.0,  90.0, 90.0,90.0,90.0, 90.0});
    const std::vector<long long> faceNodeConnectivity({0,1,6,5,1,2,7,6,2,3,8,7,3,4,9,8,4,5,0,9});
    const size_t numFaces = 5;
    const size_t numVertsPerFace = 4;

    double solidAngle = -1.0;
    bool hasWrapAroundCell = true;
    
    REQUIRE( !computeSolidAngle(lon, lat, faceNodeConnectivity, numFaces,
			           numVertsPerFace, solidAngle, hasWrapAroundCell) );

    // Compute gap left by the wrap-around cell
    const double deg2rad = vtkMath::Pi()/180.0;
    const double gap = (sin(90.0*deg2rad)-sin(-90.0*deg2rad))*((180.0-175.0)*deg2rad);
    REQUIRE( solidAngle == Approx(4.0*vtkMath::Pi()-gap) );
    REQUIRE( hasWrapAroundCell );
  }
}

// ------------------------------------------------------------------------------------------

TEST_CASE( "ResolvePeriodicGrid Test", "[basic]" )
{
  // Define 3 quads along the horizontal axis
  std::vector<double> lonOriginal({0.0, 1.0, 2.0, 3.0,
                                   0.0, 1.0, 2.0, 3.0,
                                   0.0, 1.0, 2.0, 3.0,
                                   0.0, 1.0, 2.0, 3.0});
  std::vector<double> latOriginal({0.0, 0.0, 0.0, 0.0,
                                   1.0, 1.0, 1.0, 1.0,
                                   2.0, 2.0, 2.0, 2.0,
                                   3.0, 3.0, 3.0, 3.0});
  std::vector<long long> faceNodeOriginal({0,1, 5, 4, 1, 2, 6, 5,  2, 3, 7, 6,
                                           4,5, 9, 8, 5, 6,10, 9,  6, 7,11,10,
                                           8,9,13,12, 9,10,14,13, 10,11,15,14});
  const double lonMin = *std::min_element(lonOriginal.begin(),lonOriginal.end());
  const double lonMax = *std::max_element(lonOriginal.begin(),lonOriginal.end());
  const double latMin = *std::min_element(latOriginal.begin(),latOriginal.end());
  const double latMax = *std::max_element(latOriginal.begin(),latOriginal.end());
  const bool globalModel = false;
  const size_t numVertsPerFace = 4;
  size_t numFaces = 9;

  // Copy vectors for manipulation
  std::vector<double> lon = lonOriginal;
  std::vector<double> lat = latOriginal;
  std::vector<long long> faceNode = faceNodeOriginal;

  SECTION( "Non-Periodic Planar Grid Is Not Mirrored" )
  {
    resolvePeriodicGrid(lon, lat, faceNode, numFaces, numVertsPerFace, globalModel, lonMin, lonMax, latMin, latMax);

    REQUIRE( lon.size() == lonOriginal.size() );
    REQUIRE( lat.size() == latOriginal.size() );
    REQUIRE( faceNode.size() == faceNodeOriginal.size() );

    REQUIRE( lon == lonOriginal );
    REQUIRE( lat == latOriginal );
    REQUIRE( faceNode == faceNodeOriginal );
  }

  SECTION( "Lon-Periodic Planar Grid Is Mirrored" )
  {
    // Insert wrap-around cell and update copy
    std::vector<long long> wrapAroundCell({3,0,4,7});
    faceNodeOriginal.insert(faceNodeOriginal.begin()+12,wrapAroundCell.begin(),wrapAroundCell.end());
    faceNode = faceNodeOriginal;
    numFaces += 1;

    resolvePeriodicGrid(lon, lat, faceNode, numFaces, numVertsPerFace, globalModel, lonMin, lonMax, latMin, latMax);

    // Require two extra points and the same number of cells
    REQUIRE( lon.size() == lonOriginal.size()+2 );
    REQUIRE( lat.size() == latOriginal.size()+2 );
    REQUIRE( faceNode.size() == faceNodeOriginal.size() );

    // Add mirrored points - for a planar grid, these coincide with the points with
    // highest longitude - and edit connectivity accordingly
    lonOriginal.insert(lonOriginal.end(), 2, 3.0);
    latOriginal.insert(latOriginal.end(), 0.0);
    latOriginal.insert(latOriginal.end(), 1.0);
    faceNodeOriginal.at(13) = 16;
    faceNodeOriginal.at(14) = 17;

    REQUIRE( lon == lonOriginal );
    REQUIRE( lat == latOriginal );
    REQUIRE( faceNode == faceNodeOriginal );
  }

  SECTION( "Lat-Periodic Planar Grid Is Mirrored" )
  {
    // Insert wrap-around cell and update copy
    std::vector<long long> wrapAroundCell({12,13,1,0});
    faceNodeOriginal.insert(faceNodeOriginal.begin(),wrapAroundCell.begin(),wrapAroundCell.end());
    faceNode = faceNodeOriginal;
    numFaces += 1;

    resolvePeriodicGrid(lon, lat, faceNode, numFaces, numVertsPerFace, globalModel, lonMin, lonMax, latMin, latMax);

    // Require two extra points and the same number of cells
    REQUIRE( lon.size() == lonOriginal.size()+2 );
    REQUIRE( lat.size() == latOriginal.size()+2 );
    REQUIRE( faceNode.size() == faceNodeOriginal.size() );

    // Add mirrored points and edit connectivity accordingly
    lonOriginal.insert(lonOriginal.end(), 1.0);
    lonOriginal.insert(lonOriginal.end(), 0.0);
    latOriginal.insert(latOriginal.end(), 2, 3.0);
    faceNodeOriginal.at(2) = 16;
    faceNodeOriginal.at(3) = 17;

    REQUIRE( lon == lonOriginal );
    REQUIRE( lat == latOriginal );
    REQUIRE( faceNode == faceNodeOriginal );
  }

  SECTION( "Lat-Lon-Periodic Planar Grid Is Mirrored" )
  {
    // Insert wrap-around cell and update copy
    std::vector<long long> wrapAroundCell({15,12,0,3});
    faceNodeOriginal.insert(faceNodeOriginal.end(),wrapAroundCell.begin(),wrapAroundCell.end());
    faceNode = faceNodeOriginal;
    numFaces += 1;

    resolvePeriodicGrid(lon, lat, faceNode, numFaces, numVertsPerFace, globalModel, lonMin, lonMax, latMin, latMax);

    // Require two extra points and the same number of cells
    REQUIRE( lon.size() == lonOriginal.size()+3 );
    REQUIRE( lat.size() == latOriginal.size()+3 );
    REQUIRE( faceNode.size() == faceNodeOriginal.size() );

    // Add mirrored points and edit connectivity accordingly
    lonOriginal.insert(lonOriginal.end(), 3, 3.0);
    latOriginal.insert(latOriginal.end(), 3, 3.0);
    faceNodeOriginal.at(37) = 16;
    faceNodeOriginal.at(38) = 17;
    faceNodeOriginal.at(39) = 18;

    REQUIRE( lon == lonOriginal );
    REQUIRE( lat == latOriginal );
    REQUIRE( faceNode == faceNodeOriginal );
  }
}
