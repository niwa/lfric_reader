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

  SECTION( "HasDim With Valid Dim Works" )
  {
    REQUIRE( ncFile.HasDim("Two") == true );
  }

  SECTION( "HasDim With Invalid Dim Works" )
  {
    REQUIRE( ncFile.HasDim("TwoTwo") == false );
  }

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

  SECTION( "GetAttInt Returns Correct Attribute" )
  {
    REQUIRE( ncFile.GetAttInt("var1", "index") == 7 );
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

  SECTION( "VarHasDim Returns False If Existing Dimension Not Used" )
  {
    REQUIRE( ncFile.HasDim("Four") == true );
    REQUIRE( ncFile.VarHasDim("var1", "Four") == false );
  }

  SECTION( "VarHasDim Returns True If Dimension Exists And Is Used" )
  {
    REQUIRE( ncFile.HasDim("full_levels") == true );
    REQUIRE( ncFile.VarHasDim("var1", "full_levels") == true );
  }

  SECTION( "VarHasAtt Returns False If Attribute Nonexistent" )
  {
    REQUIRE( ncFile.VarHasAtt("var1", "nonexistentatt") == false );
  }

  SECTION( "VarHasAtt Returns True If Attribute Exists" )
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

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "NetCDF Data Tests", "[data]" )
{
  netCDFLFRicFile ncFile("testdata_valid.nc");

  SECTION( "GetVarDouble Returns Correct Result" )
  {
    std::vector<double> result = ncFile.GetVarDouble("var1", {0,0}, {1,1});
    REQUIRE( result[0] == Approx(1.0) );
  }

  SECTION( "GetVarDouble Returns NaN With Incorrect Number Of Dimensions" )
  {
    std::vector<double> result = ncFile.GetVarDouble("var1", {0}, {1});
    REQUIRE( std::isnan(result[0]) == true );
    result = ncFile.GetVarDouble("var1", {0,0,0}, {1,1,1});
    REQUIRE( std::isnan(result[0]) == true );
  }

  SECTION( "GetVarLongLong Returns Correct Result" )
  {
    std::vector<long long> result =
      ncFile.GetVarLongLong("Mesh2d_full_levels_edge_nodes", {0,0}, {1,1});
    REQUIRE( result[0] == 1 );
  }

  SECTION( "GetVarLongLong Returns Zero With Incorrect Number Of Dimensions" )
  {
    std::vector<long long> result =
      ncFile.GetVarLongLong("Mesh2d_full_levels_edge_nodes", {0}, {1});
    REQUIRE( result[0] == 0 );
    result =
      ncFile.GetVarLongLong("Mesh2d_full_levels_edge_nodes", {0,0,0}, {1,1,1});
    REQUIRE( result[0] == 0 );
  }

}
