paraview_add_plugin(netCDFLFRicReader
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES vtkNetCDFLFRicReader
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")
