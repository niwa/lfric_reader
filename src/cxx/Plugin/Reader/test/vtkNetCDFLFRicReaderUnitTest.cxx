#include "Catch2/catch.hpp"
#include "vtkNetCDFLFRicReader.h"

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

  SECTION( "SetUseIndexAsVertCoord/GetUseIndexAsVertCoord Methods" )
  {
    // This parameter defaults to 0, so start with 1
    reader->SetUseIndexAsVertCoord(1);
    int result = reader->GetUseIndexAsVertCoord();
    REQUIRE( result == 1 );

    reader->SetUseIndexAsVertCoord(0);
    result = reader->GetUseIndexAsVertCoord();
    REQUIRE( result == 0 );
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

  SECTION( "SetOutputMode/GetOutputMode Methods" )
  {
    reader->SetOutputMode(1);
    int result = reader->GetOutputMode();
    REQUIRE( result == 1 );

    reader->SetOutputMode(0);
    result = reader->GetOutputMode();
    REQUIRE( result == 0 );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "CanReadFile Method", "[vtk_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "Nullpointer as filename" )
  {
    const int result = reader->CanReadFile(nullptr);
    REQUIRE( result == 0 );
  }

  SECTION( "Nonexistent File" )
  {
    const int result = reader->CanReadFile("filedoesnotexist");
    REQUIRE( result == 0 );
  }

  SECTION( "Invalid file" )
  {
    const int result = reader->CanReadFile("testdata_invalid.nc");
    REQUIRE( result == 0 );
  }

  SECTION( "Valid file" )
  {
    const int result = reader->CanReadFile("testdata_valid.nc");
    REQUIRE( result == 1 );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "GetNumberOfCellArrays Method", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "No file loaded" )
  {
    const int result = reader->GetNumberOfCellArrays();
    REQUIRE( result == 0 );
  }

  SECTION( "File loaded" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const int result = reader->GetNumberOfCellArrays();
    REQUIRE( result == 3 );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "GetCellArrayName Method", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "Invalid Index (negative)" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const char* result = reader->GetCellArrayName(-1);
    REQUIRE( result == nullptr );
  }

  SECTION( "Invalid Index (too large)" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const char* result = reader->GetCellArrayName(100);
    REQUIRE( result == nullptr );
  }

  SECTION( "Valid Index" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const std::string result(reader->GetCellArrayName(0));
    REQUIRE( result == "var1" );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "SetCellArrayStatus/GetCellArrayStatus Methods", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "GetCellArrayStatus with Invalid Array Name" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const int result = reader->GetCellArrayStatus("thisfielddoesnotexist");
    REQUIRE( result == 0 );
  }

  SECTION( "GetCellArrayStatus Default Status" )
  {
    reader->SetFileName("testdata_valid.nc");
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
    reader->SetFileName("testdata_valid.nc");
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
    reader->SetFileName("testdata_valid.nc");
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

// ------------------------------------------------------------------------------------------

TEST_CASE( "GetNumberOfPointArrays Method", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "No file loaded" )
  {
    const int result = reader->GetNumberOfPointArrays();
    REQUIRE( result == 0 );
  }

  SECTION( "File loaded" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const int result = reader->GetNumberOfPointArrays();
    REQUIRE( result == 1 );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "GetPointArrayName Method", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "Invalid Index (negative)" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const char* result = reader->GetPointArrayName(-1);
    REQUIRE( result == nullptr );
  }

  SECTION( "Invalid Index (too large)" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const char* result = reader->GetPointArrayName(100);
    REQUIRE( result == nullptr );
  }

  SECTION( "Valid Index" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const std::string result(reader->GetPointArrayName(0));
    REQUIRE( result == "var4" );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "SetPointArrayStatus/GetPointArrayStatus Methods", "[paraview_interface]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();

  SECTION( "GetPointArrayStatus with Invalid Array Name" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    const int result = reader->GetPointArrayStatus("thisfielddoesnotexist");
    REQUIRE( result == 0 );
  }

  SECTION( "GetPointArrayStatus Default Status" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    for (int iarray = 0; iarray < reader->GetNumberOfPointArrays(); iarray++)
    {
      const char* arrayname = reader->GetPointArrayName(iarray);
      const int result = reader->GetPointArrayStatus(arrayname);
      REQUIRE( result == 0 );
    }
  }

  SECTION( "SetPointArrayStatus with Invalid Array Name" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();
    reader->SetPointArrayStatus("thisfielddoesnotexist",1);
    // Check if this has affected the status any of the existing arrays
    for (int iarray = 0; iarray < reader->GetNumberOfPointArrays(); iarray++)
    {
      const char* arrayname = reader->GetPointArrayName(iarray);
      const int result = reader->GetPointArrayStatus(arrayname);
      REQUIRE( result == 0 );
    }
  }

  SECTION( "SetPointArrayStatus with Valid Array Name" )
  {
    reader->SetFileName("testdata_valid.nc");
    reader->Update();

    const int testarrayidx = 0;
    const char* testarrayname = reader->GetPointArrayName(testarrayidx);

    // Check that only status of the first array has been set
    reader->SetPointArrayStatus(testarrayname,1);
    for (int iarray = 0; iarray < reader->GetNumberOfPointArrays(); iarray++)
    {
      const char* arrayname = reader->GetPointArrayName(iarray);
      const int result = reader->GetPointArrayStatus(arrayname);
      if (iarray == testarrayidx)
      {
        REQUIRE( result == 1 );
      }
      else
      {
        REQUIRE( result == 0 );
      }
    }

    reader->SetPointArrayStatus(testarrayname,0);
    for (int iarray = 0; iarray < reader->GetNumberOfPointArrays(); iarray++)
    {
      const char* arrayname = reader->GetPointArrayName(iarray);
      const int result = reader->GetPointArrayStatus(arrayname);
      REQUIRE( result == 0 );
    }
  }

  reader->Delete();

}
