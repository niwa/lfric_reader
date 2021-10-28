#define CATCH_CONFIG_RUNNER
#include "Catch2/catch.hpp"

#include "generate_testfile.h"

int main( int argc, char * argv[] )
{

  // Create netCDF files with test data,
  // valid without and with multiple meshes, and
  // invalid
  generate_testfile(false, true);
  generate_testfile(true, true);
  generate_testfile(false, false);

  const int result = Catch::Session().run( argc, argv );
  return result;
}
