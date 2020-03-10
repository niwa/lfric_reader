#define CATCH_CONFIG_RUNNER
#include "Catch2/catch.hpp"

#include "generate_testfile.h"

int main( int argc, char * argv[] )
{

  // Create netCDF files with test data, one that
  // is valid to read and one that is not
  generate_testfile(true);
  generate_testfile(false);

  const int result = Catch::Session().run( argc, argv );
  return result;
}
