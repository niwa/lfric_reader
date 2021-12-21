# Include Catch CMake integration module
include("${CMAKE_CURRENT_SOURCE_DIR}/Catch2/cmake/Catch.cmake")

# Test executable
set(Test_SRCS catch_main.cxx generate_testfile.cxx vtkNetCDFLFRicReaderUnitTest.cxx
    vtkNetCDFLFRicReaderLeakTest.cxx vtkNetCDFLFRicReaderRegressionTest.cxx
    netCDFLFRicFileUnitTest.cxx netCDFLFRicReaderUtilsUnitTest.cxx)

add_executable (netCDFLFRicReader_tests ${Test_SRCS})
# Include auto-generated header file for DLL symbol import/export when compiling with MSVC
target_include_directories(netCDFLFRicReader_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/..)

if (ParaView_VERSION GREATER_EQUAL 5.7)

  # New build system creates module libraries which cannot be linked,
  # need to build additional shared library for testing
  set(lib_SRCS ../vtkNetCDFLFRicReader.cxx ../netCDFLFRicFile.cxx ../netCDFLFRicReaderUtils.cxx)
  add_library(netCDFLFRicReaderShared ${lib_SRCS})
  # Include auto-generated header file for DLL symbol import/export when compiling with MSVC
  # and request symbol export
  target_include_directories(netCDFLFRicReaderShared PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/..)
  target_compile_definitions(netCDFLFRicReaderShared PRIVATE vtkNetCDFLFRicReader_EXPORTS)

  target_link_libraries(netCDFLFRicReaderShared PRIVATE
    VTK::netcdf
    VTK::CommonCore
    VTK::CommonExecutionModel)

  target_link_libraries(netCDFLFRicReader_tests PRIVATE
    netCDFLFRicReaderShared
    VTK::netcdf
    VTK::CommonCore
    VTK::CommonExecutionModel
    VTK::IOParallelXML)

else()

  target_link_libraries(netCDFLFRicReader_tests PRIVATE netCDFLFRicReader)

endif()

catch_discover_tests(netCDFLFRicReader_tests)
