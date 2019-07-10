#include "Catch2/catch.hpp"
#include "vtkNetCDFLFRicReader.h"

// ------------------------------------------------------------------------------------------

TEST_CASE( "Leak Test", "[leak]" )
{
  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "Load and Reload All Arrays Several Times" )
  {
    reader->SetFileName("testdata_valid.nc");
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
