cmake_minimum_required(VERSION 3.1)

# Source code requires at least C++ 11
set (CMAKE_CXX_STANDARD 11)

set(Test_SRCS catch_main.cxx netCDFLFRicReaderUnitTest.cxx netCDFLFRicReaderLeakTest.cxx)
add_executable (netCDFLFRicReader_tests ${Test_SRCS})
target_include_directories(netCDFLFRicReader_tests PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_link_libraries(netCDFLFRicReader_tests LINK_PRIVATE netCDFLFRicReader)

option(BUILD_TESTING "Build Testing" OFF)
if(BUILD_TESTING)
  include(CTest)
  add_test (NAME Tests COMMAND netCDFLFRicReader_tests)
  add_definitions(-DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR})
endif()
