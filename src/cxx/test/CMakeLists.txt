# Include Catch CMake integration module
include("${CMAKE_CURRENT_SOURCE_DIR}/Catch2/cmake/Catch.cmake")

set(Test_SRCS catch_main.cxx generate_testfile.cxx vtkNetCDFLFRicReaderUnitTest.cxx vtkNetCDFLFRicReaderLeakTest.cxx netCDFLFRicFileUnitTest.cxx)
add_executable (netCDFLFRicReader_tests ${Test_SRCS})
target_include_directories(netCDFLFRicReader_tests PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(netCDFLFRicReader_tests LINK_PRIVATE netCDFLFRicReader)
catch_discover_tests(netCDFLFRicReader_tests)
