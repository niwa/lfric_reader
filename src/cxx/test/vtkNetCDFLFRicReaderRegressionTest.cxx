#include "Catch2/catch.hpp"
#include "vtkNetCDFLFRicReader.h"
#include "vtkCellTypes.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkXMLPUnstructuredGridWriter.h"

// ------------------------------------------------------------------------------------------

TEST_CASE( "UnstructuredGrid Properties", "[regression]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();
  reader->SetFileName("testdata_valid.nc");
  reader->Update();
  vtkUnstructuredGrid * grid = reader->GetOutput();

  SECTION( "GetNumberOfPoints Returns Correct Numbers" )
  {
    // Points are replicated for periodic grid
    REQUIRE( grid->GetNumberOfPoints() == 256 );
    reader->SetUseCartCoords(1);
    reader->Update();
    REQUIRE( grid->GetNumberOfPoints() == 224 );
  }

  SECTION( "GetNumberOfCells Returns Correct Number" )
  {
    // Test grid has 3x3x6 cells, must be the same
    // in lon-lat-rad and xyz
    REQUIRE( grid->GetNumberOfCells() == 162 );
    reader->SetUseCartCoords(1);
    reader->Update();
    REQUIRE( grid->GetNumberOfCells() == 162 );
  }

  SECTION( "Lon-Lat-Rad Grid Only Contains VTK_HEXAHEDRON Cells")
  {
    REQUIRE( grid->IsHomogeneous() == 1 );
    REQUIRE( grid->GetCellType(0) == VTK_HEXAHEDRON );
  }

  SECTION( "XYZ Grid Only Contains VTK_HEXAHEDRON Cells")
  {
    reader->SetUseCartCoords(1);
    reader->Update();
    REQUIRE( grid->IsHomogeneous() == 1 );
    REQUIRE( grid->GetCellType(0) == VTK_HEXAHEDRON );
  }

  SECTION ("Horizontal Lon-Lat-Rad Bounds Are Correct")
  {
    double gridBounds[6];
    grid->GetBounds(gridBounds);
    // Longitude bounds should add up to 360 degrees
    REQUIRE( (gridBounds[0] + gridBounds[1]) == Approx(360.0) );
    // Latitude bounds are symmetric
    REQUIRE( gridBounds[3] == Approx(-gridBounds[2]) );
    // Vertical bounds must not be symmetric
    REQUIRE( gridBounds[5] != Approx(-gridBounds[4]) );
  }

  SECTION ("XYZ Bounds Are Symmetric")
  {
    reader->SetUseCartCoords(1);
    reader->Update();
    double gridBounds[6];
    grid->GetBounds(gridBounds);
    REQUIRE( gridBounds[1] == Approx(-gridBounds[0]) );
    REQUIRE( gridBounds[3] == Approx(-gridBounds[2]) );
    REQUIRE( gridBounds[5] == Approx(-gridBounds[4]) );
  }

  SECTION ("Vertical Bounds Are Correct")
  {
    double gridBounds[6];
    grid->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(1.0) );
    REQUIRE( gridBounds[5] == Approx(2.5) );
  }

  SECTION ("Using Index As Vertical Coordinate Works")
  {
    double gridBounds[6];
    reader->SetUseIndexAsVertCoord(1);
    reader->Update();
    grid->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(1.0) );
    REQUIRE( gridBounds[5] == Approx(4.0) );
  }

  SECTION ("Setting Vertical Bias Works")
  {
    const double verticalBias = 3.7;
    reader->SetUseIndexAsVertCoord(1);
    reader->SetVerticalBias(verticalBias);
    reader->Update();
    double gridBounds[6];
    grid->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(verticalBias) );
    REQUIRE( gridBounds[5] == Approx(verticalBias+3.0) );
  }

  SECTION ("Setting Vertical Scale Works")
  {
    const double verticalScale = 9.17;
    reader->SetUseIndexAsVertCoord(1);
    reader->SetVerticalScale(verticalScale);
    reader->Update();
    double gridBounds[6];
    grid->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(verticalScale) );
    REQUIRE( gridBounds[5] == Approx(4.0*verticalScale) );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "Data Fields", "[regression]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();
  reader->SetFileName("testdata_valid.nc");
  reader->Update();
  reader->SetCellArrayStatus("var1",1);
  reader->SetCellArrayStatus("var2",1);
  reader->SetCellArrayStatus("var3",1);
  reader->Update();
  vtkUnstructuredGrid * grid = reader->GetOutput();

  SECTION( "Correct Number Of Cell Data Arrays Are Read" )
  {
    REQUIRE( grid->GetCellData()->GetNumberOfArrays() == 3 );
  }

  SECTION( "No Point Data Arrays Are Available" )
  {
    REQUIRE( grid->GetPointData()->GetNumberOfArrays() == 0 );
  }

  SECTION( "All Arrays Are Scalar" )
  {
    REQUIRE( grid->GetCellData()->GetArray(0)->GetNumberOfComponents() == 1 );
    REQUIRE( grid->GetCellData()->GetArray(1)->GetNumberOfComponents() == 1 );
    REQUIRE( grid->GetCellData()->GetArray(2)->GetNumberOfComponents() == 1 );
  }

  SECTION( "Cell Data Fields Contain Correct Data" )
  {
    double dataRange[2];
    grid->GetCellData()->GetArray(0)->GetRange(dataRange,0);
    REQUIRE( dataRange[0] == Approx(1.0) );
    REQUIRE( dataRange[1] == Approx(1.0) );
    grid->GetCellData()->GetArray(1)->GetRange(dataRange,0);
    REQUIRE( dataRange[0] == Approx(2.0) );
    REQUIRE( dataRange[1] == Approx(2.0) );
    grid->GetCellData()->GetArray(2)->GetRange(dataRange,0);
    REQUIRE( dataRange[0] == Approx(3.0) );
    REQUIRE( dataRange[1] == Approx(3.0) );
  }

  reader->Delete();

}

// ------------------------------------------------------------------------------------------

TEST_CASE( "Grid Partitioning", "[regression]" )
{

  vtkNetCDFLFRicReader * reader = vtkNetCDFLFRicReader::New();
  reader->SetFileName("testdata_valid.nc");
  reader->SetUseIndexAsVertCoord(1);
  reader->Update();

  // Use writer class to request pieces (partitions) from reader
  vtkXMLPUnstructuredGridWriter * writer = vtkXMLPUnstructuredGridWriter::New();
  writer->SetInputConnection(reader->GetOutputPort());
  writer->SetFileName("testoutput.pvtu");

  SECTION( "Single Piece Request Works" )
  {
    writer->SetNumberOfPieces(1);
    writer->SetStartPiece(0);
    writer->SetEndPiece(0);
    writer->Update();
    REQUIRE( reader->GetOutput()->GetNumberOfPieces() == 1 );
  }

  SECTION( "No Ghost Cells Are Produced By Default" )
  {
    writer->SetNumberOfPieces(1);
    writer->SetStartPiece(0);
    writer->SetEndPiece(0);
    writer->Update();
    REQUIRE( reader->GetOutput()->HasAnyGhostCells() == 0 );
  }

  SECTION( "No Ghost Points Are Produced By Default" )
  {
    writer->SetNumberOfPieces(1);
    writer->SetStartPiece(0);
    writer->SetEndPiece(0);
    writer->Update();
    REQUIRE( reader->GetOutput()->HasAnyGhostPoints() == 0 );
  }

  SECTION( "Bottom Piece Has Single Ghost Level" )
  {
    writer->SetNumberOfPieces(3);
    writer->SetStartPiece(0);
    writer->SetEndPiece(0);
    writer->SetGhostLevel(1);
    writer->Update();
    REQUIRE( reader->GetOutput()->GetGhostLevel() == 1 );
    double gridBounds[6];
    reader->GetOutput()->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(1.0) );
    REQUIRE( gridBounds[5] == Approx(3.0) );
  }

  SECTION( "Middle Piece Has Two Ghost Levels" )
  {
    writer->SetNumberOfPieces(3);
    writer->SetStartPiece(1);
    writer->SetEndPiece(1);
    writer->SetGhostLevel(1);
    writer->Update();
    REQUIRE( reader->GetOutput()->GetGhostLevel() == 1 );
    double gridBounds[6];
    reader->GetOutput()->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(1.0) );
    REQUIRE( gridBounds[5] == Approx(4.0) );
  }

  SECTION( "Top Piece Has Single Ghost Level" )
  {
    writer->SetNumberOfPieces(3);
    writer->SetStartPiece(2);
    writer->SetEndPiece(2);
    writer->SetGhostLevel(1);
    writer->Update();
    REQUIRE( reader->GetOutput()->GetGhostLevel() == 1 );
    double gridBounds[6];
    reader->GetOutput()->GetBounds(gridBounds);
    REQUIRE( gridBounds[4] == Approx(2.0) );
    REQUIRE( gridBounds[5] == Approx(4.0) );
  }

  writer->Delete();
  reader->Delete();

}
