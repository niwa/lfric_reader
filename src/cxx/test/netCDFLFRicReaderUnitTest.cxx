#include "Catch2/catch.hpp"
#include "vtkNetCDFLFRicReader.h"

// Two-level macro definition to stringize macro
#define MACRO2STRING(s) STRING(s)
#define STRING(s) #s

// ------------------------------------------------------------------------------------------

TEST_CASE( "Basic Class Tests", "[basic]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "SetFileName/GetFileName Methods" )
  {
    reader->SetFileName("my_testname");
    const std::string result(reader->GetFileName());
    REQUIRE( result == "my_testname" );
  }

  SECTION( "SetUseCartCoords/GetUseCartCoords Methods" )
  {
    // This parameter defaults to 0, so start with 1
    reader->SetUseCartCoords(1);
    int result = reader->GetUseCartCoords();
    REQUIRE( result == 1 );

    reader->SetUseCartCoords(0);
    result = reader->GetUseCartCoords();
    REQUIRE( result == 0 );
  }

  SECTION( "SetVerticalScale/GetVerticalScale Methods" )
  {
    reader->SetVerticalScale(1.23456);
    const double result = reader->GetVerticalScale();
    REQUIRE( result == Approx(1.23456) );
  }

  SECTION( "SetVerticalBias/GetVerticalBias Methods" )
  {
    reader->SetVerticalBias(1.23456);
    const double result = reader->GetVerticalBias();
    REQUIRE( result == Approx(1.23456) );
  }

  SECTION( "SetVerticalBias/GetVerticalBias Methods" )
  {
    reader->SetVerticalBias(1.23456);
    const double result = reader->GetVerticalBias();
    REQUIRE( result == Approx(1.23456) );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "CanReadFile Method", "[vtk_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  const std::string dataFilePath = MACRO2STRING(TEST_DATA_DIR);

  SECTION( "Nullpointer as filename" )
  {
    const int result = reader->CanReadFile(nullptr);
    REQUIRE( result == 0 );
  }

  SECTION( "Nonexistent file" )
  {
    const int result = reader->CanReadFile("thisfiledoesnotexist.nc");
    REQUIRE( result == 0 );
  }

  SECTION( "Invalid file" )
  {
    const std::string dataFile = dataFilePath + "/testdata_invalid.nc";
    const int result = reader->CanReadFile(dataFile.c_str());
    REQUIRE( result == 0 );
  }

  SECTION( "Valid file" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    const int result = reader->CanReadFile(dataFile.c_str());
    REQUIRE( result == 1 );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "GetNumberOfCellArrays Method", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  const std::string dataFilePath = MACRO2STRING(TEST_DATA_DIR);

  SECTION( "No file loaded" )
  {
    const int result = reader->GetNumberOfCellArrays();
    REQUIRE( result == 0 );
  }

  SECTION( "File loaded" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    const int result = reader->GetNumberOfCellArrays();
    REQUIRE( result == 5 );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "GetCellArrayName Method", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  const std::string dataFilePath = MACRO2STRING(TEST_DATA_DIR);

  SECTION( "Invalid Index (negative)" )
  {
    std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    const char* result = reader->GetCellArrayName(-1);
    REQUIRE( result == nullptr );
  }

  SECTION( "Invalid Index (too large)" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    const char* result = reader->GetCellArrayName(100);
    REQUIRE( result == nullptr );
  }

  SECTION( "Valid Index" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    const std::string result(reader->GetCellArrayName(0));
    REQUIRE( result == "buoyancy" );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "SetCellArrayStatus/GetCellArrayStatus Methods", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  const std::string dataFilePath = MACRO2STRING(TEST_DATA_DIR);

  SECTION( "GetCellArrayStatus with Invalid Array Name" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    const int result = reader->GetCellArrayStatus("thisfielddoesnotexist");
    REQUIRE( result == 0 );
  }

  SECTION( "GetCellArrayStatus Default Status" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    for (int iarray = 0; iarray < reader->GetNumberOfCellArrays(); iarray++)
    {
      const char* arrayname = reader->GetCellArrayName(iarray);
      const int result = reader->GetCellArrayStatus(arrayname);
      REQUIRE( result == 0 );
    }
  }

  SECTION( "SetCellArrayStatus with Invalid Array Name" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();
    reader->SetCellArrayStatus("thisfielddoesnotexist",1);
    // Check if this has affected the status any of the existing arrays
    for (int iarray = 0; iarray < reader->GetNumberOfCellArrays(); iarray++)
    {
      const char* arrayname = reader->GetCellArrayName(iarray);
      const int result = reader->GetCellArrayStatus(arrayname);
      REQUIRE( result == 0 );
    }
  }

  SECTION( "SetCellArrayStatus with Valid Array Name" )
  {
    const std::string dataFile = dataFilePath + "/testdata_valid.nc";
    reader->SetFileName(dataFile.c_str());
    reader->Update();

    const int testarrayidx = 0;
    const char* testarrayname = reader->GetCellArrayName(testarrayidx);

    // Check that only status of the first array has been set
    reader->SetCellArrayStatus(testarrayname,1);
    for (int iarray = 0; iarray < reader->GetNumberOfCellArrays(); iarray++)
    {
      const char* arrayname = reader->GetCellArrayName(iarray);
      const int result = reader->GetCellArrayStatus(arrayname);
      if (iarray == testarrayidx)
      {
        REQUIRE( result == 1 );
      }
      else
      {
        REQUIRE( result == 0 );
      }
    }

    reader->SetCellArrayStatus(testarrayname,0);
    for (int iarray = 0; iarray < reader->GetNumberOfCellArrays(); iarray++)
    {
      const char* arrayname = reader->GetCellArrayName(iarray);
      const int result = reader->GetCellArrayStatus(arrayname);
      REQUIRE( result == 0 );
    }
  }

  reader->Delete();

}
