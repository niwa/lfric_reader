#include "Catch2/catch.hpp"
#include "netCDFLFRicFile.h"

// ------------------------------------------------------------------------------------------

TEST_CASE( "netCDFLFRicFile Basic Class Tests", "[basic]" )
{

  SECTION( "Construct With Nonexistent File" )
  { 
    netCDFLFRicFile ncFile("nonexistentfile");
    REQUIRE( ncFile.IsFileOpen() == false );
  }

  SECTION( "Construct With Existing File" )
  { 
    netCDFLFRicFile ncFile("testdata_valid.nc");
    REQUIRE( ncFile.IsFileOpen() == true );
  }

  SECTION( "GetFileName Works" )
  { 
    netCDFLFRicFile ncFile("testdata_valid.nc");
    std::string result = ncFile.GetFileName();
    REQUIRE( result == "testdata_valid.nc" );
  }

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "NetCDF Metadata Tests", "[metadata]" )
{
  netCDFLFRicFile ncFile("testdata_valid.nc");

  SECTION( "GetDimLen Works" )
  {
    REQUIRE( ncFile.GetDimLen("Two") == 2 );
  }

  SECTION( "GetVarDimName With Out Of Range Dim Index" )
  {
    REQUIRE( ncFile.GetVarDimName("var1", -1) == "" );
    REQUIRE( ncFile.GetVarDimName("var1", 100) == "" );
  }

  SECTION( "GetVarDimName With Valid Index" )
  {
    REQUIRE( ncFile.GetVarDimName("var1", 0) == "full_levels" );
  }

  SECTION( "GetAttText Returns Correct Attribute" )
  {
    REQUIRE( ncFile.GetAttText("var1", "long_name") == "var1" );
  }

  SECTION( "GetAttTextSplit Returns Correct Attributes" )
  {
    std::vector<std::string> result =
      ncFile.GetAttTextSplit("Mesh2d_full_levels", "node_coordinates");

    REQUIRE( result.size() == 2 );
    REQUIRE( result[0] == "Mesh2d_full_levels_node_x" );
    REQUIRE( result[1] == "Mesh2d_full_levels_node_y" );
  }

  SECTION( "VarHasDim Returns False If Dimension Nonexistent" )
  {
    REQUIRE( ncFile.VarHasDim("var1", "nonexistentdim") == false );
  }

  SECTION( "VarHasDim Returns True If Dimension Exists" )
  {
    REQUIRE( ncFile.VarHasDim("var1", "full_levels") == true );
  }

  SECTION( "VarHasDim Returns False If Attribute Nonexistent" )
  {
    REQUIRE( ncFile.VarHasAtt("var1", "nonexistentatt") == false );
  }

  SECTION( "VarHasDim Returns True If Attribute Exists" )
  {
    REQUIRE( ncFile.VarHasAtt("var1", "long_name") == true );
  }

  SECTION( "GetVarNames Returns Correct Names" )
  {
    std::vector<std::string> result = ncFile.GetVarNames();
    // Just check that there is more than 1 var name as
    // the number may change, and look at the first name
    REQUIRE( result.size() > 1 );
    REQUIRE( result[0] == "Mesh2d_full_levels" );
  }

  SECTION( "GetVarDouble Returns Correct Result" )
  {
    std::vector<double> result = ncFile.GetVarDouble("var1", {0,0}, {1,1});
    REQUIRE( result[0] == Approx(1.0) );
  }

  SECTION( "GetVarLongLong Returns Correct Result" )
  {
    std::vector<long long> result =
      ncFile.GetVarLongLong("Mesh2d_full_levels_edge_nodes", {0,0}, {1,1});
    REQUIRE( result[0] == 1 );
  }
  
}
