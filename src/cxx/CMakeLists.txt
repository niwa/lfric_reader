cmake_minimum_required(VERSION 3.1)

# Source code requires at least C++ 11
set (CMAKE_CXX_STANDARD 11)

# Locate ParaView build and import CMake config, macros etc.
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

# Turn on DEBUG macro if CMAKE_BUILD_TYPE=Debug has been set
# This will enable debug output at runtime
add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

# Tell CMake to build a shared library plugin for ParaView
# Readers are only needed on the server side
ADD_PARAVIEW_PLUGIN(
    netCDFLFRicReader "1.0"
    SERVER_MANAGER_XML netCDFLFRicReader.xml
    SERVER_MANAGER_SOURCES vtkNetCDFLFRicReader.cxx 
    REQUIRED_ON_SERVER)

option(BUILD_TESTING "Build Testing" OFF)
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(test)
endif()
