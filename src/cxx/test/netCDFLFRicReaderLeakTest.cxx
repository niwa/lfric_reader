#include "Catch2/catch.hpp"
#include "vtkNetCDFLFRicReader.h"

// Two-level macro definition to stringize macro
#define MACRO2STRING(s) STRING(s)
#define STRING(s) #s

// ------------------------------------------------------------------------------------------

TEST_CASE( "Leak Test", "[leak]" )
{
  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  const std::string dataFilePath = MACRO2STRING(TEST_DATA_DIR);

  SECTION( "Load and Reload All Arrays Several Times" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();

    const int nTimes = 10;

    for (int irepeat = 0; irepeat < nTimes; irepeat++)
    {
      // Load all fields...
      for (int iarray = 0; iarray < reader->GetNumberOfCellArrays(); iarray++)
      {
        const char* arrayname = reader->GetCellArrayName(iarray);
        reader->SetCellArrayStatus(arrayname,1);
      }
      reader->Update();

      // ... and unload them...
      for (int iarray = 0; iarray < reader->GetNumberOfCellArrays(); iarray++)
      {
        const char* arrayname = reader->GetCellArrayName(iarray);
        reader->SetCellArrayStatus(arrayname,0);
      }
      reader->Update();
    }
  }

  reader->Delete();

}
