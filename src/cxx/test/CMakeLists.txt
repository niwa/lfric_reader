cmake_minimum_required(VERSION 3.1)

# Source code requires at least C++ 11
set (CMAKE_CXX_STANDARD 11)

set(UnitTest_SRCS catch_main.cxx netCDFLFRicReaderUnitTest.cxx)
add_executable (netCDFLFRicReader_unittests ${UnitTest_SRCS})
target_include_directories(netCDFLFRicReader_unittests PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_link_libraries(netCDFLFRicReader_unittests LINK_PRIVATE netCDFLFRicReader)

option(BUILD_TESTING "Build Testing" OFF)
if(BUILD_TESTING)
  include(CTest)
  add_test (NAME UnitTests COMMAND netCDFLFRicReader_unittests)
  add_definitions(-DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR})
endif()
