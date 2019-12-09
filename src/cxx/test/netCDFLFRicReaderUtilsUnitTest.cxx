#include "Catch2/catch.hpp"
#include "netCDFLFRicReaderUtils.h"

#include <vector>

// ------------------------------------------------------------------------------------------

TEST_CASE( "ResolvePeriodicGrid Test", "[basic]" )
{

  SECTION( "Non-Periodic Grid With Three Equal Cells Is Not Mirrored" )
  { 
    const std::vector<double> XOriginal({0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0});
    const std::vector<double> YOriginal({0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0});
    const std::vector<long long> faceNodeOriginal({0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6});

    std::vector<double> X = XOriginal;
    std::vector<double> Y = YOriginal;
    std::vector<long long> faceNode = faceNodeOriginal;

    const size_t numFaces = 3;
    const size_t numVertsPerFace = 4;

    resolvePeriodicGrid(X, Y, faceNode, numFaces, numVertsPerFace);

    REQUIRE( X.size() == 8 );
    REQUIRE( Y.size() == 8 );
    REQUIRE( faceNode.size() == numFaces*numVertsPerFace );

    REQUIRE( X == XOriginal );
    REQUIRE( Y == YOriginal );
    REQUIRE( faceNode == faceNodeOriginal );
  }

}
