cmake_minimum_required(VERSION 3.8)
project(sharedlibrary)

# Default to optimised build
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Adds option -DBUILD_TESTING
include(CTest)

# Locate ParaView build and import CMake config, macros etc.
FIND_PACKAGE(ParaView REQUIRED)

if (ParaView_VERSION GREATER_EQUAL 5.7)

  # New build system
  include(GNUInstallDirs)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  set(BUILD_SHARED_LIBS ON)
  paraview_plugin_scan(
    PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
    PROVIDES_PLUGINS plugins
    ENABLE_BY_DEFAULT ON)
  paraview_plugin_build(PLUGINS ${plugins})

else()

  # Old way, with directory structure based on new build system
  INCLUDE(${PARAVIEW_USE_FILE})
  add_subdirectory(Plugin/Reader)

endif()

option(BUILD_DOCUMENTATION "Build Doxygen documentation" OFF)

find_package(Doxygen)
if (DOXYGEN_FOUND)

    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../../Documentation/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    message("Doxygen build started")

    add_custom_target( doc_doxygen ALL
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM )

else (DOXYGEN_FOUND)
  message("Doxygen needs to be installed to generate the documentation")
endif (DOXYGEN_FOUND)
