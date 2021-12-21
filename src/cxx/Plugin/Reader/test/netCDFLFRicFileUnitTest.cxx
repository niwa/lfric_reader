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
    netCDFLFRicFile ncFile("testdata_single_mesh_valid.nc");
    REQUIRE( ncFile.IsFileOpen() == true );
  }

  SECTION( "GetFileName Works" )
  { 
    netCDFLFRicFile ncFile("testdata_single_mesh_valid.nc");
    std::string result = ncFile.GetFileName();
    REQUIRE( result == "testdata_single_mesh_valid.nc" );
  }

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "NetCDF Metadata Tests", "[metadata]" )
{
  netCDFLFRicFile ncFile("testdata_single_mesh_valid.nc");

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
    const int dimId = ncFile.GetDimId("Two");
    REQUIRE( ncFile.GetDimLen(dimId) == 2 );
  }

  SECTION( "GetDimId Works" )
  {
    REQUIRE( ncFile.GetDimId("Two") == 3 );
  }

  SECTION( "GetDimName Works" )
  {
    const int dimId = ncFile.GetDimId("Two");
    REQUIRE( ncFile.GetDimName(dimId) == "Two" );
  }

  SECTION( "GetVarId Works" )
  {
    REQUIRE( ncFile.GetVarId("var1") == 13 );
  }

  SECTION( "GetVarName With Valid VarId Works" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.GetVarName(varId) == "var1" );
  }

  SECTION( "GetVarName With Invalid VarId Works" )
  {
    REQUIRE( ncFile.GetVarName(-1) == "invalid varId" );
  }

  SECTION( "GetVarNumDims Works" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.GetVarNumDims(varId) == 2 );
  }

  SECTION( "GetVarDimId With Valid Dim Works" )
  {
    const int varId = ncFile.GetVarId("var1");
    const int dimId = ncFile.GetDimId("nMesh2d_face");
    REQUIRE( ncFile.GetVarDimId(varId, 1) == dimId );
  }

  SECTION( "GetVarDimId With Invalid Dim Works" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.GetVarDimId(varId, 10) == -1 );
  }

  SECTION( "GetAttInt Returns Correct Attribute" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.GetAttInt(varId, "index") == 7 );
  }

  SECTION( "GetAttText Returns Correct Attribute" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.GetAttText(varId, "long_name") == "var1" );
    REQUIRE( ncFile.GetAttText("var1", "long_name") == "var1" );
  }

  SECTION( "GetAttTextSplit Returns Correct Attributes" )
  {
    const int varId = ncFile.GetVarId("Mesh2d");
    std::vector<std::string> result =
      ncFile.GetAttTextSplit(varId, "node_coordinates");

    REQUIRE( result.size() == 2 );
    REQUIRE( result[0] == "Mesh2d_node_x" );
    REQUIRE( result[1] == "Mesh2d_node_y" );
  }

  SECTION( "HasVar With Valid Var Works" )
  {
    REQUIRE( ncFile.HasVar("var1") == true );
  }

  SECTION( "HasVar With Invalid Var Works" )
  {
    REQUIRE( ncFile.HasVar("notavalidvar") == false );
  }

  SECTION( "VarHasAtt Returns False If Attribute Nonexistent" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.VarHasAtt(varId, "nonexistentatt") == false );
  }

  SECTION( "VarHasAtt Returns True If Attribute Exists" )
  {
    const int varId = ncFile.GetVarId("var1");
    REQUIRE( ncFile.VarHasAtt(varId, "long_name") == true );
  }

  SECTION( "GetNumVars Returns Correct Result" )
  {
    REQUIRE( ncFile.GetNumVars() == 17 );
  }

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "NetCDF Data Tests", "[data]" )
{
  netCDFLFRicFile ncFile("testdata_single_mesh_valid.nc");

  SECTION( "LoadVarDouble Returns Correct Result" )
  {
    const int varId = ncFile.GetVarId("var2");
    std::vector<double> result({-1.0});
    ncFile.LoadVarDouble(varId, {0,0}, {1,1}, result.data());
    REQUIRE( result[0] == Approx(0.0) );
  }

  SECTION( "LoadVarDouble Does Not Modify Buffer With Incorrect Number Of Dims" )
  {
    const int varId = ncFile.GetVarId("var2");
    std::vector<double> result({-1.0});
    ncFile.LoadVarDouble(varId, {0}, {1}, result.data());
    REQUIRE( result[0] == Approx(-1.0) );
  }

  SECTION( "GetVarDouble Returns Correct Result" )
  {
    const int varId = ncFile.GetVarId("var2");
    std::vector<double> result = ncFile.GetVarDouble(varId, {0,0}, {1,1});
    REQUIRE( result[0] == Approx(0.0) );
  }

  SECTION( "GetVarDouble Returns NaN With Incorrect Number Of Dimensions" )
  {
    const int varId = ncFile.GetVarId("var1");
    std::vector<double> result = ncFile.GetVarDouble(varId, {0}, {1});
    REQUIRE( std::isnan(result[0]) == true );
    result = ncFile.GetVarDouble(varId, {0,0,0}, {1,1,1});
    REQUIRE( std::isnan(result[0]) == true );
  }

  SECTION( "GetVarLongLong Returns Correct Result" )
  {
    const int varId = ncFile.GetVarId("Mesh2d_edge_nodes");
    std::vector<long long> result = ncFile.GetVarLongLong(varId, {0,0}, {1,1});
    REQUIRE( result[0] == 1 );
  }

  SECTION( "GetVarLongLong Returns Zero With Incorrect Number Of Dimensions" )
  {
    const int varId = ncFile.GetVarId("Mesh2d_edge_nodes");
    std::vector<long long> result = ncFile.GetVarLongLong(varId, {0}, {1});
    REQUIRE( result[0] == 0 );
    result = ncFile.GetVarLongLong(varId, {0,0,0}, {1,1,1});
    REQUIRE( result[0] == 0 );
  }

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "UGRID, Vertical Axis, And Time Axis Tests", "[metadata]" )
{
  netCDFLFRicFile ncFile("testdata_single_mesh_valid.nc");

  SECTION( "GetMesh2DDescription Returns Correct Result" )
  {
    UGRIDMeshDescription result = ncFile.GetMesh2DDescription();
    REQUIRE( result.isLFRicXIOSFile == true );
    REQUIRE( result.numTopologies == 1 );
    REQUIRE( result.meshType == fullLevelFaceMesh );
    REQUIRE( result.numNodes == 56 );
    REQUIRE( result.numEdges == 108 );
    REQUIRE( result.numFaces == 54 );
    REQUIRE( result.numVertsPerFace == 4 );
    REQUIRE( result.numEdgesPerFace == 4 );
    REQUIRE( result.nodeDimId > -1 );
    REQUIRE( result.edgeDimId > -1 );
    REQUIRE( result.faceDimId > -1 );
    REQUIRE( result.vertsPerFaceDimId > -1 );
    REQUIRE( result.edgesPerFaceDimId > -1 );
    REQUIRE( result.faceDimIdAlt == -1 );
    REQUIRE( result.meshTopologyVarId > -1 );
    REQUIRE( result.nodeCoordXVarId > -1 );
    REQUIRE( result.nodeCoordYVarId > -1 );
    REQUIRE( result.faceNodeConnVarId > -1 );
    REQUIRE( result.faceEdgeConnVarId > -1 );
    REQUIRE( result.edgeCoordXVarId > -1 );
    REQUIRE( result.edgeCoordYVarId > -1 );
    const std::string meshTopologyVar = ncFile.GetVarName(result.meshTopologyVarId);
    REQUIRE( meshTopologyVar == "Mesh2d" );
    const std::string nodeCoordXVar = ncFile.GetVarName(result.nodeCoordXVarId);
    REQUIRE( nodeCoordXVar == "Mesh2d_node_x" );
    const std::string nodeCoordYVar = ncFile.GetVarName(result.nodeCoordYVarId);
    REQUIRE( nodeCoordYVar == "Mesh2d_node_y" );
    const std::string faceNodeConnVar = ncFile.GetVarName(result.faceNodeConnVarId);
    REQUIRE( faceNodeConnVar == "Mesh2d_face_nodes" );
    // Edge-half-level-mesh is currently used here to due inconsistency with
    // full-level mesh
    const std::string edgeCoordXVar = ncFile.GetVarName(result.edgeCoordXVarId);
    REQUIRE( edgeCoordXVar == "Mesh2d_edge_x" );
    const std::string edgeCoordYVar = ncFile.GetVarName(result.edgeCoordYVarId);
    REQUIRE( edgeCoordYVar == "Mesh2d_edge_y" );
    REQUIRE( result.faceNodeStartIdx == 1 );
  }

  SECTION( "GetZAxisDescription Returns Correct Result For Full-Level LFRic File" )
  {
    std::map<std::string, CFAxis> result = ncFile.GetZAxisDescription(true, fullLevelFaceMesh);

    // Map needs to contain these three entries
    REQUIRE( result.size() == 3 );
    REQUIRE( result.count("vtk") == 1 );
    REQUIRE( result.count("half_levels") == 1 );
    REQUIRE( result.count("full_levels") == 1 );

    // VTK axis has "full_levels" level heights, but count the number of cells
    REQUIRE( result.at("vtk").axisLength == 3 );
    std::string axisDim = ncFile.GetDimName(result.at("vtk").axisDimId);
    REQUIRE( axisDim == "full_levels" );
    std::string axisVar = ncFile.GetVarName(result.at("vtk").axisVarId);
    REQUIRE( axisVar == "full_levels" );

    // Check "full_levels" axis
    REQUIRE( result.at("full_levels").axisLength == 4 );
    axisDim = ncFile.GetDimName(result.at("full_levels").axisDimId);
    REQUIRE( axisDim == "full_levels" );
    axisVar = ncFile.GetVarName(result.at("full_levels").axisVarId);
    REQUIRE( axisVar == "full_levels" );

    // Check "half_levels" axis
    REQUIRE( result.at("half_levels").axisLength == 3 );
    axisDim = ncFile.GetDimName(result.at("half_levels").axisDimId);
    REQUIRE( axisDim == "half_levels" );
    axisVar = ncFile.GetVarName(result.at("half_levels").axisVarId);
    REQUIRE( axisVar == "half_levels" );
  }

  SECTION( "GetZAxisDescription Returns Correct Result For Half-Level LFRic File" )
  {
    std::map<std::string, CFAxis> result = ncFile.GetZAxisDescription(true, halfLevelFaceMesh);

    // Map needs to contain these three entries
    REQUIRE( result.size() == 3 );
    REQUIRE( result.count("vtk") == 1 );
    REQUIRE( result.count("half_levels") == 1 );
    REQUIRE( result.count("full_levels") == 1 );

    // VTK axis must be "half_levels" axis, counting the number of cells
    REQUIRE( result.at("vtk").axisLength == 3 );
    std::string axisDim = ncFile.GetDimName(result.at("vtk").axisDimId);
    REQUIRE( axisDim == "half_levels" );
    std::string axisVar = ncFile.GetVarName(result.at("vtk").axisVarId);
    REQUIRE( axisVar == "half_levels" );

    // Only "half_levels" axis is required
    REQUIRE( result.at("half_levels").axisLength == 3 );
    axisDim = ncFile.GetDimName(result.at("half_levels").axisDimId);
    REQUIRE( axisDim == "half_levels" );
    axisVar = ncFile.GetVarName(result.at("half_levels").axisVarId);
    REQUIRE( axisVar == "half_levels" );
  }

  SECTION( "GetZAxisDescription Returns No Axis For Non-LFRic File" )
  {
    std::map<std::string, CFAxis> result = ncFile.GetZAxisDescription(false, halfLevelFaceMesh);

    REQUIRE( result.size() == 1 );
    REQUIRE( result.count("vtk") == 1 );
    REQUIRE( result.at("vtk").axisLength == 1 );
    REQUIRE( result.at("vtk").axisDimId == -1 );
    REQUIRE( result.at("vtk").axisVarId == -1 );
  }

  SECTION( "GetZAxisDescription Returns No Axis For Unknown Mesh" )
  {
    std::map<std::string, CFAxis> result = ncFile.GetZAxisDescription(false, unknownMesh);

    REQUIRE( result.size() == 1 );
    REQUIRE( result.count("vtk") == 1 );
    REQUIRE( result.at("vtk").axisLength == 1 );
    REQUIRE( result.at("vtk").axisDimId == -1 );
    REQUIRE( result.at("vtk").axisVarId == -1 );

    result = ncFile.GetZAxisDescription(true, unknownMesh);

    REQUIRE( result.size() == 1 );
    REQUIRE( result.count("vtk") == 1 );
    REQUIRE( result.at("vtk").axisLength == 1 );
    REQUIRE( result.at("vtk").axisDimId == -1 );
    REQUIRE( result.at("vtk").axisVarId == -1 );
  }

  SECTION( "GetTAxisDescription Returns No Axis" )
  {
    CFAxis result = ncFile.GetTAxisDescription();
    REQUIRE( result.axisLength == 0 );
    REQUIRE( result.axisDimId == -1 );
    REQUIRE( result.axisVarId == -1 );
  }

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "Field Tests", "[metadata]" )
{
  netCDFLFRicFile ncFile("testdata_single_mesh_valid.nc");

  // Collect mesh metadata
  UGRIDMeshDescription mesh2D = ncFile.GetMesh2DDescription();
  CFAxis tAxis = ncFile.GetTAxisDescription();
  std::map<std::string, CFAxis> zAxes = ncFile.GetZAxisDescription(mesh2D.isLFRicXIOSFile,
                                                                   fullLevelFaceMesh);
  std::map<std::string, DataField> cellfields;
  std::map<std::string, DataField> pointfields;

  SECTION( "UpdateFieldMap Returns Correct Results For LFRic Files" )
  {
    ncFile.UpdateFieldMaps(mesh2D, zAxes, tAxis, cellfields, pointfields);

    // Expect both cell data and point data
    REQUIRE( cellfields.size() == 4 );
    REQUIRE( pointfields.size() == 1 );

    // var1 has horizontal and vertical dimensions
    REQUIRE( cellfields.count("var1") == 1 );
    DataField & field = cellfields["var1"];
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

    // var2 has horizontal and component dimensions
    REQUIRE( cellfields.count("var2") == 1 );
    field = cellfields["var2"];
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

    // var3 has horizontal and vertical dimensions
    REQUIRE( cellfields.count("var3") == 1 );
    field = cellfields["var3"];
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

    // var4 has horizontal and vertical dimensions (cell field)
    REQUIRE( cellfields.count("var4") == 1 );
    field = cellfields["var4"];
    REQUIRE( field.active == false );
    REQUIRE( field.meshType == halfLevelEdgeMesh );
    REQUIRE( field.location == cellFieldLoc );
    REQUIRE( field.hasComponentDim == false );
    REQUIRE( field.hasVerticalDim == true );
    REQUIRE( field.hasTimeDim == false );
    REQUIRE( field.dims.size() == 2 );
    REQUIRE( field.dims[0].dimType == horizontalAxisDim );
    REQUIRE( field.dims[0].dimLength == 108 );
    REQUIRE( field.dims[0].dimStride == 3 );
    REQUIRE( field.dims[1].dimType == verticalAxisDim );
    REQUIRE( field.dims[1].dimLength == 3 );
    REQUIRE( field.dims[1].dimStride == 1 );

    // var4 has horizontal and vertical dimensions (point field)
    REQUIRE( pointfields.count("var4") == 1 );
    field = pointfields["var4"];
    REQUIRE( field.active == false );
    REQUIRE( field.meshType == halfLevelEdgeMesh );
    REQUIRE( field.location == pointFieldLoc );
    REQUIRE( field.hasComponentDim == false );
    REQUIRE( field.hasVerticalDim == true );
    REQUIRE( field.hasTimeDim == false );
    REQUIRE( field.dims.size() == 2 );
    REQUIRE( field.dims[0].dimType == horizontalAxisDim );
    REQUIRE( field.dims[0].dimLength == 108 );
    REQUIRE( field.dims[0].dimStride == 3 );
    REQUIRE( field.dims[1].dimType == verticalAxisDim );
    REQUIRE( field.dims[1].dimLength == 3 );
    REQUIRE( field.dims[1].dimStride == 1 );
  }

  SECTION( "UpdateFieldMaps Treats Unknown Full_Levels Axis As Component Dim" )
  {
    // Remove full_levels axis
    zAxes.erase("full_levels");

    ncFile.UpdateFieldMaps(mesh2D, zAxes, tAxis, cellfields, pointfields);

    // Expect all variables
    REQUIRE( cellfields.size() == 4 );
    REQUIRE( pointfields.size() == 1 );
  }

  SECTION( "UpdateFieldMaps Treats Unknown Half_Levels Axis As Component Dim" )
  {
    // Remove half_levels axis
    zAxes.erase("half_levels");

    ncFile.UpdateFieldMaps(mesh2D, zAxes, tAxis, cellfields, pointfields);

    // Expect all variables
    REQUIRE( cellfields.size() == 4 );
    REQUIRE( pointfields.size() == 1 );
  }

  SECTION( "UpdateFieldMaps Ignores Vars With Unidentified Dimensions" )
  {
    mesh2D.faceDimId = -1;

    ncFile.UpdateFieldMaps(mesh2D, zAxes, tAxis, cellfields, pointfields);

    // Only edge field must be returned
    REQUIRE( cellfields.size() == 1 );
    REQUIRE( pointfields.size() == 1 );
  }

  SECTION( "UpdateDataFields Returns No Fields For Undefined Metadata" )
  {
    // Reset metadata
    mesh2D = UGRIDMeshDescription();
    tAxis = CFAxis();
    zAxes.clear();

    // File is treated as a general UGRID file, which requires metadata to
    // identify all field dimensions and accept field variables
    ncFile.UpdateFieldMaps(mesh2D, zAxes, tAxis, cellfields, pointfields);
    REQUIRE( cellfields.size() == 0 );
    REQUIRE( pointfields.size() == 0 );
  }
}
