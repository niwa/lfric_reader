#include <vtk_netcdf.h>
#include <iostream>
#include <string>
#include <vector>

#define ncErrorMacro(e) if (e != NC_NOERR) {std::cerr << "NetCDF Error: " << nc_strerror(e) << "\n";}

void generate_testfile(const bool multiple_mesh, const bool valid)
{
  // Function can generate test files with single mesh or multiple meshes (to test LFRic output), as
  // well as invalid files that will not be accepted by LFRic reader
  std::string fileName = "testdata";
  if (multiple_mesh)
  {
    fileName += "_multiple_mesh";
  }
  else
  {
    fileName += "_single_mesh";
  }
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

  const size_t nMesh2dNodeLen = 56;
  int nMesh2dNodeDimId;
  ncErrorMacro(nc_def_dim(ncId, "nMesh2d_node", nMesh2dNodeLen, &nMesh2dNodeDimId));

  const size_t nMesh2dEdgeLen = 108;
  int nMesh2dEdgeDimId;
  ncErrorMacro(nc_def_dim(ncId, "nMesh2d_edge", nMesh2dEdgeLen, &nMesh2dEdgeDimId));

  const size_t nMesh2dFaceLen = 54;
  int nMesh2dFaceDimId;
  ncErrorMacro(nc_def_dim(ncId, "nMesh2d_face", nMesh2dFaceLen, &nMesh2dFaceDimId));

  // Separate edge-half-levels mesh
  int nMesh2dEdgeHalfLevelsNodeDimId = -1;
  int nMesh2dEdgeHalfLevelsEdgeDimId = -1;
  if (multiple_mesh)
  {
    ncErrorMacro(nc_def_dim(ncId, "nMesh2d_edge_half_levels_node", nMesh2dNodeLen,
                            &nMesh2dEdgeHalfLevelsNodeDimId));

    ncErrorMacro(nc_def_dim(ncId, "nMesh2d_edge_half_levels_edge", nMesh2dEdgeLen,
                            &nMesh2dEdgeHalfLevelsEdgeDimId));
  }

  // Helper dims
  const size_t TwoLen = 2;
  int TwoDimId;
  ncErrorMacro(nc_def_dim(ncId, "Two", TwoLen, &TwoDimId));

  const size_t FourLen = 4;
  int FourDimId;
  ncErrorMacro(nc_def_dim(ncId, "Four", FourLen, &FourDimId));

  // Time, vertical levels, and field components

  const size_t timeCounterLen = 0;
  int timeCounterDimId;
  ncErrorMacro(nc_def_dim(ncId, "time_counter", timeCounterLen, &timeCounterDimId));

  const size_t levelsLen = 3;
  int halfLevelsDimId;
  ncErrorMacro(nc_def_dim(ncId, "half_levels", levelsLen, &halfLevelsDimId));

  int fullLevelsDimId;
  ncErrorMacro(nc_def_dim(ncId, "full_levels", levelsLen+1, &fullLevelsDimId));

  const size_t componentLen = 3;
  int componentDimId;
  ncErrorMacro(nc_def_dim(ncId, "some_component", componentLen, &componentDimId));

  // ------------------------------------------------------------------------------------------

  // Define variables

  // Define UGRID mesh description variables

  std::string meshVarName;
  if (multiple_mesh)
  {
    meshVarName = "Mesh2d_face";
  }
  else
  {
    meshVarName = "Mesh2d";
  }

  // Faces mesh

  // Valid file needs to have dummy variable "Mesh2d" with UGRID metadata
  int Mesh2dId;
  ncErrorMacro(nc_def_var(ncId, meshVarName.c_str(), NC_INT, 0, 0, &Mesh2dId));

  // Mesh connectivity

  const int Mesh2dFaceNodesDims[] = {nMesh2dFaceDimId, FourDimId};
  int Mesh2dFaceNodesId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_face_nodes", NC_INT, 2,
                          Mesh2dFaceNodesDims, &Mesh2dFaceNodesId));

  const int Mesh2dEdgeNodesDims[] = {nMesh2dEdgeDimId, TwoDimId};
  int Mesh2dEdgeNodesId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_nodes", NC_INT, 2,
                          Mesh2dEdgeNodesDims, &Mesh2dEdgeNodesId));

  const int Mesh2dFaceEdgesDims[] = {nMesh2dFaceDimId, FourDimId};
  int Mesh2dFaceEdgesId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_face_edges", NC_INT, 2,
                          Mesh2dFaceEdgesDims, &Mesh2dFaceEdgesId));

  const int Mesh2dFaceLinksDims[] = {nMesh2dFaceDimId, FourDimId};
  int Mesh2dFaceLinksId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_face_links", NC_INT, 2,
                          Mesh2dFaceLinksDims, &Mesh2dFaceLinksId));

  // Mesh element locations node - edge - face

  const int Mesh2dNodeXDims[] = {nMesh2dNodeDimId};
  int Mesh2dNodeXId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_node_x", NC_DOUBLE, 1,
                          Mesh2dNodeXDims, &Mesh2dNodeXId));

  const int Mesh2dNodeYDims[] = {nMesh2dNodeDimId};
  int Mesh2dNodeYId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_node_y", NC_DOUBLE, 1,
                          Mesh2dNodeYDims, &Mesh2dNodeYId));
    
  const int Mesh2dEdgeXDims[] = {nMesh2dEdgeDimId};
  int Mesh2dEdgeXId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_x", NC_DOUBLE, 1,
                          Mesh2dEdgeXDims, &Mesh2dEdgeXId));
  
  const int Mesh2dEdgeYDims[] = {nMesh2dEdgeDimId};
  int Mesh2dEdgeYId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_y", NC_DOUBLE, 1,
                          Mesh2dEdgeYDims, &Mesh2dEdgeYId));

  const int Mesh2dFaceXDims[] = {nMesh2dFaceDimId};
  int Mesh2dFaceXId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_face_x", NC_DOUBLE, 1,
                          Mesh2dFaceXDims, &Mesh2dFaceXId));
  
  const int Mesh2dFaceYDims[] = {nMesh2dFaceDimId};
  int Mesh2dFaceYId;
  ncErrorMacro(nc_def_var(ncId, "Mesh2d_face_y", NC_DOUBLE, 1,
                          Mesh2dFaceYDims, &Mesh2dFaceYId));

  // Additional connectivity and element locations for separate edge mesh

  int Mesh2dEdgeHalfLevelsId = -1;
  int Mesh2dEdgeHalfLevelsEdgeNodesId = -1;
  int Mesh2dEdgeHalfLevelsNodeXId = -1;
  int Mesh2dEdgeHalfLevelsNodeYId = -1;
  int Mesh2dEdgeHalfLevelsEdgeXId = -1;
  int Mesh2dEdgeHalfLevelsEdgeYId = -1;
  if (multiple_mesh)
  {

    // Edge-half-levels
    const std::string edgeHalfLevelMeshVarName = "Mesh2d_edge_half_levels";
    ncErrorMacro(nc_def_var(ncId, edgeHalfLevelMeshVarName.c_str(), NC_INT, 0,
                            0, &Mesh2dEdgeHalfLevelsId));

    const int Mesh2dEdgeHalfLevelsEdgeNodesDims[] = {nMesh2dEdgeHalfLevelsEdgeDimId, TwoDimId};
    ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_half_levels_edge_nodes", NC_INT, 2,
                            Mesh2dEdgeHalfLevelsEdgeNodesDims,
                            &Mesh2dEdgeHalfLevelsEdgeNodesId));

    const int Mesh2dEdgeHalfLevelsNodeXDims[] = {nMesh2dEdgeHalfLevelsNodeDimId};
    ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_half_levels_node_x", NC_DOUBLE, 1,
                            Mesh2dEdgeHalfLevelsNodeXDims,
                            &Mesh2dEdgeHalfLevelsNodeXId));

    const int Mesh2dEdgeHalfLevelsNodeYDims[] = {nMesh2dEdgeHalfLevelsNodeDimId};
    ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_half_levels_node_y", NC_DOUBLE, 1,
                            Mesh2dEdgeHalfLevelsNodeYDims,
                            &Mesh2dEdgeHalfLevelsNodeYId));

    const int Mesh2dEdgeHalfLevelsEdgeXDims[] = {nMesh2dEdgeHalfLevelsEdgeDimId};
    ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_half_levels_edge_x", NC_DOUBLE, 1,
                            Mesh2dEdgeHalfLevelsEdgeXDims,
                            &Mesh2dEdgeHalfLevelsEdgeXId));

    const int Mesh2dEdgeHalfLevelsEdgeYDims[] = {nMesh2dEdgeHalfLevelsEdgeDimId};
    ncErrorMacro(nc_def_var(ncId, "Mesh2d_edge_half_levels_edge_y", NC_DOUBLE, 1,
                            Mesh2dEdgeHalfLevelsEdgeYDims,
                            &Mesh2dEdgeHalfLevelsEdgeYId));
  }

  // Vertical levels

  const int FullLevelsDims[] = {fullLevelsDimId};
  int FullLevelsId;
  ncErrorMacro(nc_def_var(ncId, "full_levels", NC_FLOAT, 1, FullLevelsDims,
                          &FullLevelsId));

  const int HalfLevelsDims[] = {halfLevelsDimId};
  int HalfLevelsId;
  ncErrorMacro(nc_def_var(ncId, "half_levels", NC_FLOAT, 1, HalfLevelsDims,
                          &HalfLevelsId));

  // Set UGRID mesh variable attributes
  std::string attStr = "Topology data of 2D unstructured mesh";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "long_name", attStr.length(),
                               attStr.c_str()));

  if (valid)
  {
    attStr = "mesh_topology";
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "cf_role", attStr.length(),
                                 attStr.c_str()));
  }
  else
  {
    attStr = "something_else";
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "cf_role", attStr.length(),
                                 attStr.c_str()));
  }

  // Topology dimension = 2 means faces
  const int topDimAtt[] = {2};
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dId, "topology_dimension",
                              NC_INT, 1, topDimAtt));
  attStr = "Mesh2d_node_x Mesh2d_node_y";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "node_coordinates", attStr.length(),
                               attStr.c_str()));
  attStr = "Mesh2d_face_nodes";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "face_node_connectivity", attStr.length(),
                               attStr.c_str()));
  attStr = "Mesh2d_edge_nodes";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "edge_node_connectivity", attStr.length(),
                               attStr.c_str()));
  attStr = "Mesh2d_face_x Mesh2d_face_y";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "face_coordinates", attStr.length(),
                               attStr.c_str()));
  attStr = "Mesh2d_face_edges";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "face_edge_connectivity", attStr.length(),
                               attStr.c_str()));
  attStr = "Mesh2d_face_links";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "face_face_connectivity", attStr.length(),
                               attStr.c_str()));
  attStr = "Mesh2d_edge_x Mesh2d_edge_y";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dId, "edge_coordinates", attStr.length(),
                               attStr.c_str()));
  attStr = "face_node_connectivity";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceNodesId, "cf_role", attStr.length(),
                               attStr.c_str()));
  attStr = "Maps every quadrilateral face to its four corner nodes.";
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceNodesId, "long_name", attStr.length(),
                               attStr.c_str()));
  const int startIndexAtt[] = {1};
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFaceNodesId, "start_index",
                              NC_INT, 1, startIndexAtt));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeNodesId, "cf_role", 22,
                               "edge_node_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeNodesId, "long_name", 50,
                               "Maps every edge to the two nodes that it connects."));
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dEdgeNodesId, "start_index",
                              NC_INT, 1, startIndexAtt));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceEdgesId, "cf_role", 22,
                               "face_edge_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceEdgesId, "long_name", 48,
                               "Maps every quadrilateral face to its four edges."));
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFaceEdgesId, "start_index",
                              NC_INT, 1, startIndexAtt));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceLinksId, "cf_role", 22,
                               "face_face_connectivity"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceLinksId, "long_name", 48,
                               "Indicates which other faces neighbour each face."));
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFaceLinksId, "start_index",
                              NC_INT, 1, startIndexAtt));
  const int flagValuesAtt[] = {-1};
  ncErrorMacro(nc_put_att_int(ncId, Mesh2dFaceLinksId, "flag_values",
                              NC_INT, 1, flagValuesAtt));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceLinksId, "flag_meanings", 11,
                               "out_of_mesh"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dNodeXId, "standard_name", 9,
                               "longitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dNodeXId, "long_name", 27,
                               "longitude of 2D mesh nodes."));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dNodeXId, "units", 12,
                               "degrees_east"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dNodeYId, "standard_name", 8,
                               "latitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dNodeYId, "long_name", 26,
                               "latitude of 2D mesh nodes."));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dNodeYId, "units", 13,
                               "degrees_north"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeXId, "standard_name", 9,
                               "longitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeXId, "long_name", 28,
                               "longitude of 2D edge centres"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeXId, "units", 12,
                               "degrees_east"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeYId, "standard_name", 8,
                               "latitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeYId, "long_name", 27,
                               "latitude of 2D edge centres"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeYId, "units", 13,
                               "degrees_north"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceXId, "standard_name", 9,
                               "longitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceXId, "long_name", 28,
                               "longitude of 2D face centres"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceXId, "units", 12,
                               "degrees_east"));

  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceYId, "standard_name", 8,
                               "latitude"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceYId, "long_name", 27,
                               "latitude of 2D face centres"));
  ncErrorMacro(nc_put_att_text(ncId, Mesh2dFaceYId, "units", 13,
                               "degrees_north"));

  // Edge-half-levels
  if (multiple_mesh)
  {
    if (valid)
    {
      ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsId, "cf_role", 13,
                                   "mesh_topology"));
    }
    else
    {
      ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsId, "cf_role", 14,
                                   "something_else"));
    }
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsId, "long_name", 37,
                                 "Topology data of 2D unstructured mesh"));

    // LFRic output files currently set topology dimension = 2 for edge mesh
    const int topDimAttEdge[] = {2};
    ncErrorMacro(nc_put_att_int(ncId, Mesh2dEdgeHalfLevelsId, "topology_dimension",
                                NC_INT, 1, topDimAttEdge));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsId, "node_coordinates", 61,
                                 "Mesh2d_edge_half_levels_node_x Mesh2d_edge_half_levels_node_y"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsId, "edge_node_connectivity", 34,
                                 "Mesh2d_edge_half_levels_edge_nodes"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsId, "edge_coordinates", 61,
                                 "Mesh2d_edge_half_levels_edge_x Mesh2d_edge_half_levels_edge_y"));
  
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeNodesId, "cf_role", 22,
                                 "edge_node_connectivity"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeNodesId, "long_name", 50,
                                 "Maps every edge to the two nodes that it connects."));
    ncErrorMacro(nc_put_att_int(ncId, Mesh2dEdgeHalfLevelsEdgeNodesId, "start_index",
                                NC_INT, 1, startIndexAtt));
  
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsNodeXId, "standard_name", 9,
                                 "longitude"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsNodeXId, "long_name", 27,
                                 "longitude of 2D mesh nodes."));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsNodeXId, "units", 12,
                                 "degrees_east"));
  
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsNodeYId, "standard_name", 8,
                                 "latitude"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsNodeYId, "long_name", 26,
                                 "latitude of 2D mesh nodes."));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsNodeYId, "units", 13,
                                 "degrees_north"));
  
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeXId, "standard_name", 9,
                                 "longitude"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeXId, "long_name", 28,
                                 "longitude of 2D edge centres"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeXId, "units", 12,
                                 "degrees_east"));
  
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeYId, "standard_name", 8,
                                 "latitude"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeYId, "long_name", 27,
                                 "latitude of 2D edge centres"));
    ncErrorMacro(nc_put_att_text(ncId, Mesh2dEdgeHalfLevelsEdgeYId, "units", 13,
                                 "degrees_north"));
  }

  // Test data variables and attributes

  const int var1Dims[] = {halfLevelsDimId, nMesh2dFaceDimId};
  int var1Id;
  ncErrorMacro(nc_def_var(ncId, "var1", NC_DOUBLE, 2, var1Dims, &var1Id));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "long_name", 4, "var1"));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "mesh", 16, "value not needed"));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "location", 4, "face"));
  ncErrorMacro(nc_put_att_text(ncId, var1Id, "coordinates", 16, "value not needed"));
  const int intAttValue = 7;
  ncErrorMacro(nc_put_att_int(ncId, var1Id, "index", NC_INT, 1, &intAttValue));

  const int var2Dims[] = {nMesh2dFaceDimId, componentDimId};
  int var2Id;
  ncErrorMacro(nc_def_var(ncId, "var2", NC_DOUBLE, 2, var2Dims, &var2Id));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "long_name", 4, "var2"));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "mesh", 16, "value not needed"));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "location", 4, "face"));
  ncErrorMacro(nc_put_att_text(ncId, var2Id, "coordinates", 16, "value not needed"));

  const int var3Dims[] = {nMesh2dFaceDimId, fullLevelsDimId};
  int var3Id;
  ncErrorMacro(nc_def_var(ncId, "var3", NC_DOUBLE, 2, var3Dims, &var3Id));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "long_name", 4, "var3"));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "mesh", 16, "value not needed"));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "location", 4, "face"));
  ncErrorMacro(nc_put_att_text(ncId, var3Id, "coordinates", 16, "value not needed"));

  int var4Id = -1;
  if (multiple_mesh)
  {
    const int var4Dims[] = {nMesh2dEdgeHalfLevelsEdgeDimId, halfLevelsDimId};
    ncErrorMacro(nc_def_var(ncId, "var4", NC_DOUBLE, 2, var4Dims, &var4Id));  }
  else
  {
    const int var4Dims[] = {nMesh2dEdgeDimId, halfLevelsDimId};
    ncErrorMacro(nc_def_var(ncId, "var4", NC_DOUBLE, 2, var4Dims, &var4Id));
  }
  ncErrorMacro(nc_put_att_text(ncId, var4Id, "long_name", 4, "var4"));
  ncErrorMacro(nc_put_att_text(ncId, var4Id, "units", 4, "none"));
  ncErrorMacro(nc_put_att_text(ncId, var4Id, "mesh", 16, "value not needed"));
  ncErrorMacro(nc_put_att_text(ncId, var4Id, "location", 4, "edge"));
  ncErrorMacro(nc_put_att_text(ncId, var4Id, "coordinates", 16, "value not needed"));

  // End of define mode
  ncErrorMacro(nc_enddef(ncId));

  const size_t start = 0;
  const size_t start2[] = {0,0};
  size_t count = 0;
  size_t count2[2];

  // Faces mesh
  const int Mesh2d_data[] = {-2147483647};
  count = 1;
  ncErrorMacro(nc_put_var1(ncId, Mesh2dId, &count, Mesh2d_data));

  const int Mesh2d_face_nodes_data[] =
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
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFaceNodesId, start2, count2,
                           Mesh2d_face_nodes_data));

  const int Mesh2d_edge_nodes_data[] =
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
  ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeNodesId, start2, count2,
                           Mesh2d_edge_nodes_data));

  const int Mesh2d_face_edges_data[] =
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
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFaceEdgesId, start2, count2,
                           Mesh2d_face_edges_data));

  const int Mesh2d_face_links_data[] =
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
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFaceLinksId, start2, count2,
                           Mesh2d_face_links_data));

  const double Mesh2d_node_x_data[] =
    {315, 345, 15, 315, 345, 15, 315, 345, 15, 45, 75, 105, 45, 75, 105,
     45, 75, 105, 135, 165, 195, 135, 165, 195, 135, 165, 195, 225, 255,
     285, 225, 255, 285, 225, 255, 285, 315, 225, 45, 135, 45, 75, 105,
     15, 45, 135, 345, 315, 225, 315, 285, 255, 135, 165, 195, 225};
  count = 56;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dNodeXId, &start, &count,
                           Mesh2d_node_x_data));

  const double Mesh2d_node_y_data[] =
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
  ncErrorMacro(nc_put_vara(ncId, Mesh2dNodeYId, &start, &count,
                           Mesh2d_node_y_data));

  // Note: edge coordinates are arithmetic means of edge-node coordinates
  const double Mesh2d_edge_x_data[] =
    {330, 315, 330, 180, 345, 180,  30,  15,  30, 315, 330, 345, 180,  15,
     30, 315, 330, 345, 180,  15,  30,  60,  45,  60,  90,  75,  90, 120,
     105, 120,  45,  60,  75,  90, 105, 120,  45,  60,  75,  90, 105, 120,
     150, 135, 150, 180, 165, 180, 210, 195, 210, 135, 150, 165, 180, 195,
     210, 135, 150, 165, 180, 195, 210, 240, 225, 240, 270, 255, 270, 300,
     285, 300, 225, 240, 255, 270, 285, 300, 225, 240, 255, 270, 285, 300,
     330, 270, 210,  30,  90, 150, 300, 240, 180, 180,  60, 120,  30,  90,
     150, 330, 270, 210,  60, 120, 180, 180, 300, 240};
  count = 108;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeXId, &start, &count,
                           Mesh2d_edge_x_data));

  const double Mesh2d_edge_y_data[] =
    {39.63570844, 22.9964864, 12.61970091, 44.0070272, 29.25892295,
     14.5108187, 39.63570844, 29.25892295, 12.61970091, 0.,
     -12.61970091, 0., -14.5108187, 0., -12.61970091,
     -22.9964864, -39.63570844, -29.25892295, -44.0070272, -29.25892295,
     -39.63570844, 39.63570844, 22.9964864, 12.61970091, 44.0070272,
     29.25892295, 14.5108187, 39.63570844, 29.25892295, 12.61970091,
     0., -12.61970091, 0., -14.5108187, 0.,
     -12.61970091, -22.9964864, -39.63570844, -29.25892295, -44.0070272,
     -29.25892295, -39.63570844, 39.63570844, 22.9964864, 12.61970091,
     44.0070272, 29.25892295, 14.5108187, 39.63570844, 29.25892295,
     12.61970091, 0., -12.61970091, 0., -14.5108187,
     0., -12.61970091, -22.9964864, -39.63570844, -29.25892295,
     -44.0070272, -29.25892295, -39.63570844, 39.63570844, 22.9964864,
     12.61970091, 44.0070272, 29.25892295, 14.5108187, 39.63570844,
     29.25892295, 12.61970091, 0., -12.61970091, 0.,
     -14.5108187, 0., -12.61970091, -22.9964864, -39.63570844,
     -29.25892295, -44.0070272, -29.25892295, -39.63570844, 56.62672811,
     69.24642902, 56.62672811, 56.62672811, 69.24642902, 56.62672811,
     56.62672811, 56.62672811, 69.24642902, 69.24642902, 56.62672811,
     56.62672811, -56.62672811, -69.24642902, -56.62672811, -56.62672811,
     -69.24642902, -56.62672811, -56.62672811, -56.62672811, -69.24642902,
     -69.24642902, -56.62672811, -56.62672811};
  count = 108;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeYId, &start, &count,
                           Mesh2d_edge_y_data));

  const double Mesh2d_face_x_data[] =
    {329.508305958902, 0, 30.4916940410979, 329.886509860388, 0, 30.1134901396116,
     329.508305958902, 2.9271427301566e-15, 30.4916940410979, 59.5083059589021,
     90, 120.491694041098, 59.8865098603884, 90, 120.113490139612, 59.5083059589021,
     90, 120.491694041098, 149.508305958902, 180, 210.491694041098, 149.886509860388,
     180, 210.113490139612, 149.508305958902, 180, 210.491694041098, 239.508305958902,
     270, 300.491694041098, 239.886509860388, 270, 300.113490139612, 239.508305958902,
     270, 300.491694041098, 315, 270, 225, 360, 296.565051177078, 180, 45, 90, 135,
     45, 90, 135, 7.5702706973617e-15, 63.434948822922, 180, 315, 270, 225};
  count = 54;
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFaceXId, &start, &count,
                           Mesh2d_face_x_data));

  const double Mesh2d_face_y_data[] =
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
  ncErrorMacro(nc_put_vara(ncId, Mesh2dFaceYId, &start, &count,
                           Mesh2d_face_y_data));

  // Edge-half-levels mesh - reuse full-levels data
  if (multiple_mesh)
  {
    count = 1;
    ncErrorMacro(nc_put_var1(ncId, Mesh2dEdgeHalfLevelsId, &count, Mesh2d_data));

    count2[0] = 108;
    count2[1] = 2;
    ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeHalfLevelsEdgeNodesId, start2, count2,
                             Mesh2d_edge_nodes_data));

    count = 56;
    ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeHalfLevelsNodeXId, &start, &count,
                             Mesh2d_node_x_data));
    ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeHalfLevelsNodeYId, &start, &count,
                             Mesh2d_node_y_data));

    count = 108;
    ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeHalfLevelsEdgeXId, &start, &count,
                             Mesh2d_edge_x_data));
    ncErrorMacro(nc_put_vara(ncId, Mesh2dEdgeHalfLevelsEdgeYId, &start, &count,
                             Mesh2d_edge_y_data));
  }

  // Vertical axes
  const float FullLevelsData[] = {0.0, 0.5, 1.0, 1.5};
  count = 4;
  ncErrorMacro(nc_put_vara(ncId, FullLevelsId, &start, &count, FullLevelsData));

  const float HalfLevelsData[] = {0.25, 0.75, 1.25};
  count = 3;
  ncErrorMacro(nc_put_vara(ncId, HalfLevelsId, &start, &count, HalfLevelsData));

  // Test variables

  std::vector<double> varData;

  // Half-level data
  count2[0] = levelsLen;
  count2[1] = nMesh2dFaceLen;
  varData.resize(count2[0]*count2[1]);

  // Store dimensions to help testing
  varData[0] = static_cast<double>(nMesh2dFaceLen);
  varData[1] = static_cast<double>(levelsLen);
  varData[2] = static_cast<double>(componentLen);

  // Fill the rest of the array with number sequence to test cell ordering
  for (size_t idx = 3; idx < varData.size(); idx++)
  {
    varData[idx] = static_cast<double>(idx);
  }
  ncErrorMacro(nc_put_vara_double(ncId, var1Id, start2, count2, varData.data()));

  // Half-level 2D data with multiple components
  count2[0] = nMesh2dFaceLen;
  count2[1] = componentLen;
  varData.resize(count2[0]*count2[1]);

  // Fill array with number sequence to test component ordering
  for (size_t idx = 0; idx < varData.size(); idx++)
  {
    varData[idx] = static_cast<double>(idx);
  }
  ncErrorMacro(nc_put_vara_double(ncId, var2Id, start2, count2, varData.data()));

  // Full-level 3D data with inverse dimension order
  count2[0] = nMesh2dFaceLen;
  count2[1] = levelsLen+1;
  varData.resize(count2[0]*count2[1]);

  // Fill each level with the same number to test level ordering
  size_t bufferIdx = 0;
  for (size_t iFace = 0; iFace < nMesh2dFaceLen; iFace++)
  {
    for (size_t iLevel = 0; iLevel < levelsLen+1; iLevel++)
    {
      varData[bufferIdx] = static_cast<double>(iLevel);
      bufferIdx++;
    }
  }
  ncErrorMacro(nc_put_vara_double(ncId, var3Id, start2, count2, varData.data()));

  // Edge-half-level 3D data with inverse dimension order
  count2[0] = nMesh2dEdgeLen;
  count2[1] = levelsLen;
  varData.resize(count2[0]*count2[1]);
  
  // Fill array with number sequence to test ordering
  for (size_t iLevel = 0; iLevel < levelsLen; iLevel++)
  {
    for (size_t iEdge = 0; iEdge < nMesh2dEdgeLen; iEdge++)
    {
      const size_t iCell = iLevel*nMesh2dEdgeLen+iEdge;
      varData[iEdge*levelsLen+iLevel] = static_cast<double>(iCell);
    }
  }
  
  // Store dimensions in the first 3 cells
  varData[0] = static_cast<double>(nMesh2dEdgeLen);
  varData[levelsLen] = static_cast<double>(levelsLen);
  varData[2*levelsLen] = static_cast<double>(componentLen);
  ncErrorMacro(nc_put_vara_double(ncId, var4Id, start2, count2, varData.data()));

  ncErrorMacro(nc_close(ncId));
}
