#include <vtk_netcdf.h>
#include <iostream>
#include <string>
#include <vector>

#define ncErrorMacro(e) if (e != NC_NOERR) {std::cerr << "NetCDF Error: " << nc_strerror(e) << "\n";}

void generate_testfile(const bool valid)
{
  // Function can generate "invalid" test files that will not be accepted by LFRic reader
  std::string fileName = "testdata";
  if (valid)
  {
    fileName += "_valid.nc";
  }
  else
  {
    fileName += "_invalid.nc";
  }
  int ncId;
  ncErrorMacro(nc_create(fileName.c_str(), NC_CLOBBER, &ncId));

  // ------------------------------------------------------------------------------------------

  // Define dimensions

  // UGRID

  const size_t nMesh2dFullNevelsNodeLen = 56;
  int nMesh2dFullLevelsNodeDimId;
  ncErrorMacro(nc_def_dim(ncId, "nMesh2d_full_levels_node", nMesh2dFullNevelsNodeLen,
                          &nMesh2dFullLevelsNodeDimId));

  const size_t nMesh2dFullLevelsEdgeLen = 108;
  int nMesh2dFullLevelsEdgeDimId;
  ncErrorMacro(nc_def_dim(ncId, "nMesh2d_full_levels_edge", nMesh2dFullLevelsEdgeLen,
                          &nMesh2dFullLevelsEdgeDimId));

  const size_t nMesh2dFullLevelsFaceLen = 54;
  int nMesh2dFullLevelsFaceDimId;
  ncErrorMacro(nc_def_dim(ncId, "nMesh2d_full_levels_face", nMesh2dFullLevelsFaceLen,
                          &nMesh2dFullLevelsFaceDimId));

  const size_t TwoLen = 2;
  int TwoDimId;
  ncErrorMacro(nc_def_dim(ncId, "Two", TwoLen, &TwoDimId));

  const size_t FourLen = 4;
  int FourDimId;
  ncErrorMacro(nc_def_dim(ncId, "Four", FourLen, &FourDimId));

  // Time and vertical levels

  const size_t timeCounterLen = 0;
  int timeCounterDimId;
  ncErrorMacro(nc_def_dim(ncId, "time_counter", timeCounterLen, &timeCounterDimId));

  const size_t levelsLen = 2;
  int levelsDimId;
  ncErrorMacro(nc_def_dim(ncId, "full_levels", levelsLen, &levelsDimId));

  // ------------------------------------------------------------------------------------------

  // Define variables

  // Define UGRID mesh description variables

  // Valid file needs to have dummy variable "Mesh2d_full_levels" with UGRID metadata
  int Mesh2dFullLevelsId;
  std::string meshVarName = "Mesh2d_full_levels";
  ncErrorMacro(nc_def_var(ncId, meshVarName.c_str(), NC_INT, 0,
                          0, &Mesh2dFullLevelsId));

  const int Mesh2dFullLevelsFaceNodesDims[] = {nMesh2dFullLevelsFaceDimId, FourDimId};
  int Mesh2dFullLevelsFaceNodesId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_face_nodes", NC_INT, 2,
                          Mesh2dFullLevelsFaceNodesDims,
                          &Mesh2dFullLevelsFaceNodesId));

  const int Mesh2dFullLevelsEdgeNodesDims[] = {nMesh2dFullLevelsEdgeDimId, TwoDimId};
  int Mesh2dFullLevelsEdgeNodesId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_edge_nodes", NC_INT, 2,
                          Mesh2dFullLevelsEdgeNodesDims,
                          &Mesh2dFullLevelsEdgeNodesId));

  const int Mesh2dFullLevelsFaceEdgesDims[] = {nMesh2dFullLevelsFaceDimId, FourDimId};
  int Mesh2dFullLevelsFaceEdgesId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_face_edges", NC_INT, 2,
                          Mesh2dFullLevelsFaceEdgesDims,
                          &Mesh2dFullLevelsFaceEdgesId));

  const int Mesh2dFullLevelsFaceLinksDims[] = {nMesh2dFullLevelsFaceDimId, FourDimId};
  int Mesh2dFullLevelsFaceLinksId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_face_links", NC_INT, 2,
                          Mesh2dFullLevelsFaceLinksDims,
                          &Mesh2dFullLevelsFaceLinksId));

  const int Mesh2dFullLevelsNodeXDims[] = {nMesh2dFullLevelsNodeDimId};
  int Mesh2dFullLevelsNodeXId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_node_x", NC_DOUBLE, 1,
                          Mesh2dFullLevelsNodeXDims,
                          &Mesh2dFullLevelsNodeXId));

  const int Mesh2dFullLevelsNodeYDims[] = {nMesh2dFullLevelsNodeDimId};
  int Mesh2dFullLevelsNodeYId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_node_y", NC_DOUBLE, 1,
                          Mesh2dFullLevelsNodeYDims,
                          &Mesh2dFullLevelsNodeYId));

  const int Mesh2dFullLevelsFaceXDims[] = {nMesh2dFullLevelsFaceDimId};
  int Mesh2dFullLevelsFaceXId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_face_x", NC_DOUBLE, 1,
                          Mesh2dFullLevelsFaceXDims,
                          &Mesh2dFullLevelsFaceXId));

  const int Mesh2dFullLevelsFaceYDims[] = {nMesh2dFullLevelsFaceDimId};
  int Mesh2dFullLevelsFaceYId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_face_y", NC_DOUBLE, 1,
                          Mesh2dFullLevelsFaceYDims,
                          &Mesh2dFullLevelsFaceYId));

  const int Mesh2dFullLevelsEdgeXDims[] = {nMesh2dFullLevelsEdgeDimId};
  int Mesh2dFullLevelsEdgeXId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_edge_x", NC_DOUBLE, 1,
                          Mesh2dFullLevelsEdgeXDims,
                          &Mesh2dFullLevelsEdgeXId));

  const int Mesh2dFullLevelsEdgeYDims[] = {nMesh2dFullLevelsEdgeDimId};
  int Mesh2dFullLevelsEdgeYId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_full_levels_edge_y", NC_DOUBLE, 1,
                          Mesh2dFullLevelsEdgeYDims,
                          &Mesh2dFullLevelsEdgeYId));

  const int FullLevelsDims[] = {levelsDimId};
  int FullLevelsId;
  ncErrorMacro(nc_def_var(ncId, "full_levels", NC_FLOAT, 1, FullLevelsDims,
                          &FullLevelsId));

  // Set UGRID mesh variable attributes
  if (valid)
  {
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "cf_role", 13,
                                 "mesh_topology"));
  }
  else
  {
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "cf_role", 14,
                                 "something_else"));
  }
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "long_name", 37,
                               "Topology data of 2D unstructured mesh"));

  const int topDimAtt[] = {2};
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFullLevelsId, "topology_dimension",
                              NC_INT, 1, topDimAtt));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "node_coordinates", 51,
                               "Mesh2d_full_levels_node_x Mesh2d_full_levels_node_y"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "face_node_connectivity", 29,
                               "Mesh2d_full_levels_face_nodes"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "edge_node_connectivity", 29,
                               "Mesh2d_full_levels_edge_nodes"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "face_coordinates", 51,
                               "Mesh2d_full_levels_face_x Mesh2d_full_levels_face_y"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "face_edge_connectivity", 29,
                               "Mesh2d_full_levels_face_edges"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "face_face_connectivity", 29,
                               "Mesh2d_full_levels_face_links"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsId, "edge_coordinates", 51,
                               "Mesh2d_full_levels_edge_x Mesh2d_full_levels_edge_y"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceNodesId, "cf_role", 22,
                               "face_node_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceNodesId, "long_name", 55,
                               "Maps every quadrilateral face to its four corner nodes."));
  const int startIndexAtt[] = {1};
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFullLevelsFaceNodesId, "start_index",
                              NC_INT, 1, startIndexAtt));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsEdgeNodesId, "cf_role", 22,
                               "edge_node_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsEdgeNodesId, "long_name", 50,
                               "Maps every edge to the two nodes that it connects."));
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFullLevelsEdgeNodesId, "start_index",
                              NC_INT, 1, startIndexAtt));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceEdgesId, "cf_role", 22,
                               "face_edge_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceEdgesId, "long_name", 48,
                               "Maps every quadrilateral face to its four edges."));
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFullLevelsFaceEdgesId, "start_index",
                              NC_INT, 1, startIndexAtt));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceLinksId, "cf_role", 22,
                               "face_face_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceLinksId, "long_name", 48,
                               "Indicates which other faces neighbour each face."));
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFullLevelsFaceLinksId, "start_index",
                              NC_INT, 1, startIndexAtt));
  const int flagValuesAtt[] = {-1};
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFullLevelsFaceLinksId, "flag_values",
                              NC_INT, 1, flagValuesAtt));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceLinksId, "flag_meanings", 11,
                               "out_of_mesh"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsNodeXId, "standard_name", 9,
                               "longitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsNodeXId, "long_name", 27,
                               "longitude of 2D mesh nodes."));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsNodeXId, "units", 12,
                               "degrees_east"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsNodeYId, "standard_name", 8,
                               "latitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsNodeYId, "long_name", 26,
                               "latitude of 2D mesh nodes."));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsNodeYId, "units", 13,
                               "degrees_north"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceXId, "standard_name", 9,
                               "longitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceXId, "long_name", 28,
                               "longitude of 2D face centres"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceXId, "units", 12,
                               "degrees_east"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceYId, "standard_name", 8,
                               "latitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceYId, "long_name", 27,
                               "latitude of 2D face centres"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFullLevelsFaceYId, "units", 13,
                               "degrees_north"));

  // Test data variables and attributes

  const int var1Dims[] = {levelsDimId, nMesh2dFullLevelsFaceDimId};
  int var1Id;
  ncErrorMacro(nc_def_var(ncId, "var1", NC_DOUBLE, 2, var1Dims, &var1Id));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "long_name", 4, "var1"));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "mesh", meshVarName.size(), meshVarName.c_str()));
  const int intAttValue = 7;
  ncErrorMacro(nc_put_att_int(ncId, var1Id, "index", NC_INT, 1, &intAttValue));

  const int var2Dims[] = {levelsDimId, nMesh2dFullLevelsFaceDimId};
  int var2Id;
  ncErrorMacro(nc_def_var(ncId, "var2", NC_DOUBLE, 2, var2Dims, &var2Id));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "long_name", 4, "var2"));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "mesh", meshVarName.size(), meshVarName.c_str()));

  const int var3Dims[] = {levelsDimId, nMesh2dFullLevelsFaceDimId};
  int var3Id;
  ncErrorMacro(nc_def_var(ncId, "var3", NC_DOUBLE, 2, var3Dims, &var3Id));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "long_name", 4, "var3"));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "mesh", meshVarName.size(), meshVarName.c_str()));

  // End of define mode
  ncErrorMacro(nc_enddef(ncId));

  const size_t start = 0;
  const size_t start2[] = {0,0};
  size_t count = 0;
  size_t count2[2];

  static int Mesh2d_full_levels_data[] = {-2147483647};
  ncErrorMacro(nc_put_var1(ncId, Mesh2dFullLevelsId, &count, Mesh2d_full_levels_data));

  int Mesh2d_full_levels_face_nodes_data[] =
    {4, 5, 2, 1, 5, 6, 3, 2, 6, 13, 10, 3, 7, 8, 5, 4, 8, 9, 6, 5, 9, 16, 13, 6, 50,
     47, 8, 7, 47, 44, 9, 8, 44, 41, 16, 9, 13, 14, 11, 10, 14, 15, 12, 11, 15, 22,
     19, 12, 16, 17, 14, 13, 17, 18, 15, 14, 18, 25, 22, 15, 41, 42, 17, 16, 42, 43,
     18, 17, 43, 53, 25, 18, 23, 20, 19, 22, 24, 21, 20, 23, 31, 28, 21, 24, 26, 23,
     22, 25, 27, 24, 23, 26, 34, 31, 24, 27, 54, 26, 25, 53, 55, 27, 26, 54, 56, 34,
     27, 55, 32, 29, 28, 31, 33, 30, 29, 32, 4, 1, 30, 33, 35, 32, 31, 34, 36, 33, 32,
     35, 7, 4, 33, 36, 52, 35, 34, 56, 51, 36, 35, 52, 50, 7, 36, 51, 1, 2, 37, 30, 30,
     37, 38, 29, 29, 38, 21, 28, 2, 3, 39, 37, 37, 39, 40, 38, 38, 40, 20, 21, 3, 10,
     11, 39, 39, 11, 12, 40, 40, 12, 19, 20, 44, 45, 42, 41, 45, 46, 43, 42, 46, 54,
     53, 43, 47, 48, 45, 44, 48, 49, 46, 45, 49, 55, 54, 46, 50, 51, 48, 47, 51, 52,
     49, 48, 52, 56, 55, 49} ;
  count2[0] = 54;
  count2[1] = 4;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsFaceNodesId, start2, count2,
                           Mesh2d_full_levels_face_nodes_data));

  int Mesh2d_full_levels_edge_nodes_data[] =
    {1, 2, 4, 1, 5, 4, 2, 3, 5, 2, 6, 5, 3, 10, 6, 3, 13, 6, 7, 4, 8, 7, 8, 5, 9, 8,
     9, 6, 16, 9, 50, 7, 47, 50, 47, 8, 44, 47, 44, 9, 41, 44, 10, 11, 13, 10, 14, 13,
     11, 12, 14, 11, 15, 14, 12, 19, 15, 12, 22, 15, 16, 13, 17, 16, 17, 14, 18, 17,
     18, 15, 25, 18, 41, 16, 42, 41, 42, 17, 43, 42, 43, 18, 53, 43, 19, 20, 22, 19,
     23, 22, 20, 21, 23, 20, 24, 23, 21, 28, 24, 21, 31, 24, 25, 22, 26, 25, 26, 23,
     27, 26, 27, 24, 34, 27, 53, 25, 54, 53, 54, 26, 55, 54, 55, 27, 56, 55, 28, 29,
     31, 28, 32, 31, 29, 30, 32, 29, 33, 32, 30, 1, 33, 30, 4, 33, 34, 31, 35, 34, 35,
     32, 36, 35, 36, 33, 7, 36, 56, 34, 52, 56, 52, 35, 51, 52, 51, 36, 50, 51, 37, 2,
     38, 37, 21, 38, 39, 3, 40, 39, 20, 40, 30, 37, 29, 38, 37, 39, 38, 40, 39, 11, 40,
     12, 45, 44, 46, 45, 54, 46, 48, 47, 49, 48, 55, 49, 42, 45, 43, 46, 45, 48, 46, 49,
     48, 51, 49, 52} ;
  count2[0] = 108;
  count2[1] = 2;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsEdgeNodesId, start2, count2,
                           Mesh2d_full_levels_edge_nodes_data));

  int Mesh2d_full_levels_face_edges_data[] =
    {2, 3, 5, 1, 5, 6, 8, 4, 8, 9, 23, 7, 10, 11, 12, 3, 12, 13, 14, 6, 14, 15, 31, 9,
     16, 17, 18, 11, 18, 19, 20, 13, 20, 21, 37, 15, 23, 24, 26, 22, 26, 27, 29, 25, 29,
     30, 44, 28, 31, 32, 33, 24, 33, 34, 35, 27, 35, 36, 52, 30, 37, 38, 39, 32, 39, 40,
     41, 34, 41, 42, 58, 36, 45, 47, 43, 44, 48, 50, 46, 47, 51, 65, 49, 50, 53, 54, 45,
     52, 55, 56, 48, 54, 57, 73, 51, 56, 59, 60, 53, 58, 61, 62, 55, 60, 63, 79, 57, 62,
     66, 68, 64, 65, 69, 71, 67, 68, 72, 2, 70, 71, 74, 75, 66, 73, 76, 77, 69, 75, 78,
     10, 72, 77, 80, 81, 74, 79, 82, 83, 76, 81, 84, 16, 78, 83, 70, 1, 85, 91, 67, 91,
     86, 92, 64, 92, 87, 49, 85, 4, 88, 93, 86, 93, 89, 94, 87, 94, 90, 46, 88, 7, 22,
     95, 89, 95, 25, 96, 90, 96, 28, 43, 21, 97, 103, 38, 103, 98, 104, 40, 104, 99, 59,
     42, 19, 100, 105, 97, 105, 101, 106, 98, 106, 102, 61, 99, 17, 84, 107, 100, 107,
     82, 108, 101, 108, 80, 63, 102} ;
  count2[0] = 54;
  count2[1] = 4;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsFaceEdgesId, start2, count2,
                           Mesh2d_full_levels_face_edges_data));

  int Mesh2d_full_levels_face_links_data[] =
    {30, 4, 2, 37, 1, 5, 3, 40, 2, 6, 10, 43, 33, 7, 5, 1, 4, 8, 6, 2, 5, 9, 13, 3, 36,
     52, 8, 4, 7, 49, 9, 5, 8, 46, 16, 6, 3, 13, 11, 43, 10, 14, 12, 44, 11, 15, 19, 45,
     6, 16, 14, 10, 13, 17, 15, 11, 14, 18, 22, 12, 9, 46, 17, 13, 16, 47, 18, 14, 17,
     48, 25, 15, 22, 20, 45, 12, 23, 21, 42, 19, 24, 28, 39, 20, 25, 23, 19, 15, 26, 24,
     20, 22, 27, 31, 21, 23, 48, 26, 22, 18, 51, 27, 23, 25, 54, 34, 24, 26, 31, 29, 39,
     21, 32, 30, 38, 28, 33, 1, 37, 29, 34, 32, 28, 24, 35, 33, 29, 31, 36, 4, 30, 32, 54,
     35, 31, 27, 53, 36, 32, 34, 52, 7, 33, 35, 30, 1, 40, 38, 29, 37, 41, 39, 28, 38, 42,
     21, 37, 2, 43, 41, 38, 40, 44, 42, 39, 41, 45, 20, 40, 3, 10, 44, 41, 43, 11, 45, 42,
     44, 12, 19, 9, 49, 47, 16, 46, 50, 48, 17, 47, 51, 25, 18, 8, 52, 50, 46, 49, 53, 51,
     47, 50, 54, 26, 48, 7, 36, 53, 49, 52, 35, 54, 50, 53, 34, 27, 51};
  count2[0] = 54;
  count2[1] = 4;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsFaceLinksId, start2, count2,
                           Mesh2d_full_levels_face_links_data));

  double Mesh2d_full_levels_node_x_data[] =
    {315, 345, 15, 315, 345, 15, 315, 345, 15, 45, 75, 105, 45, 75, 105,
     45, 75, 105, 135, 165, 195, 135, 165, 195, 135, 165, 195, 225, 255,
     285, 225, 255, 285, 225, 255, 285, 315, 225, 45, 135, 45, 75, 105,
     15, 45, 135, 345, 315, 225, 315, 285, 255, 135, 165, 195, 225};
  count = 56;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsNodeXId, &start, &count,
                           Mesh2d_full_levels_node_x_data));

  double Mesh2d_full_levels_node_y_data[] =
    {35.2643896827547, 44.0070271956363, 44.0070271956363, 10.7285831216091,
     14.5108186990699, 14.5108186990699, -10.7285831216091, -14.5108186990699,
     -14.5108186990699, 35.2643896827547, 44.0070271956363, 44.0070271956363,
     10.7285831216091, 14.5108186990699, 14.5108186990699, -10.7285831216091,
     -14.5108186990699, -14.5108186990699, 35.2643896827547, 44.0070271956363,
     44.0070271956363, 10.7285831216091, 14.5108186990699, 14.5108186990699,
     -10.7285831216091, -14.5108186990699, -14.5108186990699, 35.2643896827547,
     44.0070271956363, 44.0070271956363, 10.7285831216091, 14.5108186990699,
     14.5108186990699, -10.7285831216091, -14.5108186990699, -14.5108186990699,
     69.2464290163152, 69.2464290163152, 69.2464290163152, 69.2464290163152,
     -35.2643896827547, -44.0070271956363, -44.0070271956363, -44.0070271956363,
     -69.2464290163152, -69.2464290163152, -44.0070271956363, -69.2464290163152,
     -69.2464290163152, -35.2643896827547, -44.0070271956363, -44.0070271956363,
     -35.2643896827547, -44.0070271956363, -44.0070271956363, -35.2643896827547};
  count = 56;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsNodeYId, &start, &count,
                           Mesh2d_full_levels_node_y_data));

  double Mesh2d_full_levels_face_x_data[] =
    {329.508305958902, 0, 30.4916940410979, 329.886509860388, 0, 30.1134901396116,
     329.508305958902, 2.9271427301566e-15, 30.4916940410979, 59.5083059589021,
     90, 120.491694041098, 59.8865098603884, 90, 120.113490139612, 59.5083059589021,
     90, 120.491694041098, 149.508305958902, 180, 210.491694041098, 149.886509860388,
     180, 210.113490139612, 149.508305958902, 180, 210.491694041098, 239.508305958902,
     270, 300.491694041098, 239.886509860388, 270, 300.113490139612, 239.508305958902,
     270, 300.491694041098, 315, 270, 225, 360, 296.565051177078, 180, 45, 90, 135,
     45, 90, 135, 7.5702706973617e-15, 63.434948822922, 180, 315, 270, 225};
  count = 54;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsFaceXId, &start, &count,
                           Mesh2d_full_levels_face_x_data));

  double Mesh2d_full_levels_face_y_data[] =
    {26.9038500851278, 30.1134901396116, 26.9038500851278, 2.53207598583004e-15,
     3.40125124033106e-15, 2.53207598583004e-15, -26.9038500851278, -30.1134901396116,
     -26.9038500851278, 26.9038500851278, 30.1134901396116, 26.9038500851278,
     2.53207598583004e-15, 3.40125124033106e-15, 2.53207598583004e-15, -26.9038500851278,
     -30.1134901396116, -26.9038500851278, 26.9038500851278, 30.1134901396116,
     26.9038500851278, 2.95408865013505e-15, 3.40125124033106e-15, 2.53207598583004e-15,
     -26.9038500851278, -30.1134901396116, -26.9038500851278, 26.9038500851278,
     30.1134901396116, 26.9038500851278, 2.95408865013505e-15, 3.40125124033106e-15,
     2.53207598583004e-15, -26.9038500851278, -30.1134901396116, -26.9038500851278,
     50.213843782473, 59.8865098603884, 50.213843782473, 59.8865098603884, 90,
     59.8865098603884, 50.213843782473, 59.8865098603884, 50.213843782473,
     -50.213843782473, -59.8865098603884, -50.213843782473, -59.8865098603884,
     -90, -59.8865098603884, -50.213843782473, -59.8865098603884, -50.213843782473};
  count = 54;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFullLevelsFaceYId, &start, &count,
                           Mesh2d_full_levels_face_y_data));

  float FullLevelsData[] = {0.0, 0.5};
  count = 2;
  ncErrorMacro(nc_put_vara(ncId, FullLevelsId, &start, &count, FullLevelsData));

  // Test variables

  std::vector<double> varData;
  varData.resize(nMesh2dFullLevelsFaceLen*levelsLen);

  count2[0] = levelsLen;
  count2[1] = nMesh2dFullLevelsFaceLen;

  std::fill(varData.begin(), varData.end(), 1.0);
  ncErrorMacro(nc_put_vara_double(ncId, var1Id, start2, count2, varData.data()));

  std::fill(varData.begin(), varData.end(), 2.0);
  ncErrorMacro(nc_put_vara_double(ncId, var2Id, start2, count2, varData.data()));

  std::fill(varData.begin(), varData.end(), 3.0);
  ncErrorMacro(nc_put_vara_double(ncId, var3Id, start2, count2, varData.data()));

  ncErrorMacro(nc_close(ncId));
}
