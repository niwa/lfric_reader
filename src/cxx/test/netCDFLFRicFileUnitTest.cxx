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

  SECTION( "GetVarNumDims Works" )
  {
    REQUIRE( ncFile.GetVarNumDims("var1") == 2 );
  }

  SECTION( "GetVarDimName With Out Of Range Dim Index" )
  {
    REQUIRE( ncFile.GetVarDimName("var1", -1) == "" );
    REQUIRE( ncFile.GetVarDimName("var1", 100) == "" );
  }

  SECTION( "GetVarDimName With Valid Index" )
  {
    REQUIRE( ncFile.GetVarDimName("var1", 0) == "half_levels" );
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

  SECTION( "HasVar With Valid Var Works" )
  {
    REQUIRE( ncFile.HasVar("var1") == true );
  }

  SECTION( "HasVar With Invalid Var Works" )
  {
    REQUIRE( ncFile.HasVar("notavalidvar") == false );
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
    REQUIRE( ncFile.HasDim("half_levels") == true );
    REQUIRE( ncFile.VarHasDim("var1", "half_levels") == true );
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
    std::vector<double> result = ncFile.GetVarDouble("var2", {0,0}, {1,1});
    REQUIRE( result[0] == Approx(0.0) );
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

// ------------------------------------------------------------------------------------------

TEST_CASE( "UGRID, Vertical Axis, And Time Axis Tests", "[metadata]" )
{
  netCDFLFRicFile ncFile("testdata_valid.nc");

  SECTION( "GetMesh2DDescription Returns Correct Result" )
  {
    UGRIDMeshDescription result = ncFile.GetMesh2DDescription();
    REQUIRE( result.isLFRicXIOSFile == true );
    REQUIRE( result.numTopologies == 2 );
    REQUIRE( result.meshType == fullLevelFaceMesh );
    REQUIRE( result.numNodes == 56 );
    REQUIRE( result.numEdges == 108 );
    REQUIRE( result.numFaces == 54 );
    REQUIRE( result.numVertsPerFace == 4 );
    REQUIRE( result.meshTopologyVar == "Mesh2d_full_levels" );
    REQUIRE( result.nodeCoordXVar == "Mesh2d_full_levels_node_x" );
    REQUIRE( result.nodeCoordYVar == "Mesh2d_full_levels_node_y" );
    REQUIRE( result.faceNodeConnVar == "Mesh2d_full_levels_face_nodes" );
    // Edge-half-level-mesh is currently used here to due inconsistency with
    // full-level mesh
    REQUIRE( result.edgeCoordXVar == "Mesh2d_edge_half_levels_edge_x" );
    REQUIRE( result.edgeCoordYVar == "Mesh2d_edge_half_levels_edge_y" );
    REQUIRE( result.faceNodeStartIdx == 1 );
  }

  SECTION( "GetZAxisDescription Returns Correct Result For Full-Level LFRic File" )
  {
    CFVerticalAxis result = ncFile.GetZAxisDescription(true, fullLevelFaceMesh);
    REQUIRE( result.numLevels == 3 );
    REQUIRE( result.axisDim == "full_levels" );
    REQUIRE( result.axisVar == "full_levels" );
  }

  SECTION( "GetZAxisDescription Returns Correct Result For Half-Level LFRic File" )
  {
    CFVerticalAxis result = ncFile.GetZAxisDescription(true, halfLevelFaceMesh);
    REQUIRE( result.numLevels == 3 );
    REQUIRE( result.axisDim == "half_levels" );
    REQUIRE( result.axisVar == "half_levels" );
  }

  SECTION( "GetZAxisDescription Returns No Axis For Non-LFRic File" )
  {
    CFVerticalAxis result = ncFile.GetZAxisDescription(false, halfLevelFaceMesh);
    REQUIRE( result.numLevels == 1 );
    REQUIRE( result.axisDim == "None" );
    REQUIRE( result.axisVar == "None" );
  }

  SECTION( "GetZAxisDescription Returns No Axis For Unknown Mesh" )
  {
    CFVerticalAxis result = ncFile.GetZAxisDescription(false, unknownMesh);
    REQUIRE( result.numLevels == 1 );
    REQUIRE( result.axisDim == "None" );
    REQUIRE( result.axisVar == "None" );
    result = ncFile.GetZAxisDescription(true, unknownMesh);
    REQUIRE( result.numLevels == 1 );
    REQUIRE( result.axisDim == "None" );
    REQUIRE( result.axisVar == "None" );
  }

  SECTION( "GetTAxisDescription Returns No Axis" )
  {
    CFTimeAxis result = ncFile.GetTAxisDescription();
    REQUIRE( result.numTimeSteps == 0 );
    REQUIRE( result.axisDim == "None" );
    REQUIRE( result.axisVar == "None" );
  }

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "Field Tests", "[metadata]" )
{
  netCDFLFRicFile ncFile("testdata_valid.nc");

  SECTION( "UpdateFieldMap Returns Correct Result For Full-Levels Face Fields" )
  {
    std::map<std::string, DataField> fields;
    ncFile.UpdateFieldMap(fields, "face", "nMesh2d_full_levels_face",
                          fullLevelFaceMesh, "full_levels", "None");
    REQUIRE( fields.size() == 1 );

    // var3 has horizontal and vertical dimensions
    REQUIRE( fields.count("var3") == 1 );
    DataField & field = fields["var3"];
    REQUIRE( field.active == false );
    REQUIRE( field.meshType == fullLevelFaceMesh );
    REQUIRE( field.location == cellFieldLoc );
    REQUIRE( field.hasComponentDim == false );
    REQUIRE( field.hasVerticalDim == true );
    REQUIRE( field.hasTimeDim == false );
    REQUIRE( field.dims.size() == 2 );
    REQUIRE( field.dims[0].dimType == horizontalAxisDim );
    REQUIRE( field.dims[0].dimLength == 54 );
    REQUIRE( field.dims[0].dimStride == 4 );
    REQUIRE( field.dims[1].dimType == verticalAxisDim );
    REQUIRE( field.dims[1].dimLength == 4 );
    REQUIRE( field.dims[1].dimStride == 1 );
  }

  SECTION( "UpdateFieldMap Returns Correct Result For Half-Levels Face Fields" )
  {
    std::map<std::string, DataField> fields;
    ncFile.UpdateFieldMap(fields, "face", "nMesh2d_half_levels_face",
                          halfLevelFaceMesh, "half_levels", "None");
    REQUIRE( fields.size() == 2 );

    // var1 has horizontal and vertical dimensions
    REQUIRE( fields.count("var1") == 1 );
    DataField & field = fields["var1"];
    REQUIRE( field.active == false );
    REQUIRE( field.meshType == halfLevelFaceMesh );
    REQUIRE( field.location == cellFieldLoc );
    REQUIRE( field.hasComponentDim == false );
    REQUIRE( field.hasVerticalDim == true );
    REQUIRE( field.hasTimeDim == false );
    REQUIRE( field.dims.size() == 2 );
    REQUIRE( field.dims[0].dimType == verticalAxisDim );
    REQUIRE( field.dims[0].dimLength == 3 );
    REQUIRE( field.dims[0].dimStride == 54 );
    REQUIRE( field.dims[1].dimType == horizontalAxisDim );
    REQUIRE( field.dims[1].dimLength == 54 );
    REQUIRE( field.dims[1].dimStride == 1 );

    // var1 has horizontal and component dimensions
    REQUIRE( fields.count("var2") == 1 );
    field = fields["var2"];
    REQUIRE( field.active == false );
    REQUIRE( field.meshType == halfLevelFaceMesh );
    REQUIRE( field.location == cellFieldLoc );
    REQUIRE( field.hasComponentDim == true );
    REQUIRE( field.hasVerticalDim == false );
    REQUIRE( field.hasTimeDim == false );
    REQUIRE( field.dims.size() == 2 );
    REQUIRE( field.dims[0].dimType == horizontalAxisDim );
    REQUIRE( field.dims[0].dimLength == 54 );
    REQUIRE( field.dims[0].dimStride == 3 );
    REQUIRE( field.dims[1].dimType == componentAxisDim );
    REQUIRE( field.dims[1].dimLength == 3 );
    REQUIRE( field.dims[1].dimStride == 1 );
  }

  SECTION( "UpdateFieldMap Returns No Result For Undefined Location" )
  {
    std::map<std::string, DataField> fields;
    ncFile.UpdateFieldMap(fields, "nowhere", "nMesh2d_half_levels_face",
                          halfLevelFaceMesh, "half_levels", "None");
    REQUIRE( fields.size() == 0 );
  }

  SECTION( "UpdateFieldMap Returns No Result For Two Unidentified Dimensions" )
  {
    std::map<std::string, DataField> fields;
    ncFile.UpdateFieldMap(fields, "nowhere", "nohorizontaldim",
                          halfLevelFaceMesh, "noverticaldim", "None");
    REQUIRE( fields.size() == 0 );
  }

}
