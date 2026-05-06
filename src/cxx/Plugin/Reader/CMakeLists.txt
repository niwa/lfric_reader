# ParaView requires at least C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Turn on DEBUG macro if CMAKE_BUILD_TYPE=Debug has been set
# This will enable debug output at runtime
add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

set(xmls netCDFLFRicReader.xml)
set(smsrcs vtkNetCDFLFRicReader.cxx)
set(srcs netCDFLFRicFile.cxx netCDFLFRicReaderUtils.cxx)
set(hdrs vtkNetCDFLFRicReader.h netCDFLFRicFile.h netCDFLFRicReaderUtils.h)

# Plugin build setup - tell CMake to build VTK module
vtk_module_add_module(vtkNetCDFLFRicReader
  CLASSES vtkNetCDFLFRicReader
  SOURCES ${smsrcs} ${srcs}
  PRIVATE_HEADERS ${hdrs})
vtk_module_definitions(vtkNetCDFLFRicReader
  PRIVATE HAVE_LIBNETCDF HAVE_NETCDF4)
paraview_add_server_manager_xmls(XMLS ${xmls})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
