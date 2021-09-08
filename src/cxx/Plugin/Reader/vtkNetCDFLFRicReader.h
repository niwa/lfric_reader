/*
 * @class   vtkNetCDFLFRicReader
 * @brief   Read netCDF files with LFRic output in ugrid format.
 *
 * Reads in a netCDF file with LFRic output in ugrid format and
 * produces an unstructured grid. The grid is unstructured in
 * horizontal dimension and extruded in the vertical.
*/

#ifndef vtkNetCDFLFRicReader_h
#define vtkNetCDFLFRicReader_h

#include "netCDFLFRicFile.h"

#include <vtkIONetCDFModule.h> // For export macro
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtk_netcdf.h>

#include <vector>
#include <map>
#include <string>

class VTKIONETCDF_EXPORT vtkNetCDFLFRicReader : public vtkUnstructuredGridAlgorithm
{

public:

  vtkTypeMacro(vtkNetCDFLFRicReader,vtkUnstructuredGridAlgorithm)
  static vtkNetCDFLFRicReader *New();
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Specify name of input data file
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Check file validity
  virtual int CanReadFile(const char* fileName);

  // ParaView interface for using index as vertical coordinate
  void SetUseIndexAsVertCoord(const int status);
  vtkGetMacro(UseIndexAsVertCoord, int);

  // ParaView interface for switching coordinates
  void SetUseCartCoords(const int status);
  vtkGetMacro(UseCartCoords, int);

  // ParaView interface for controlling radius coordinate
  void SetVerticalScale(const double value);
  vtkGetMacro(VerticalScale, double);
  void SetVerticalBias(const double value);
  vtkGetMacro(VerticalBias, double);

  // ParaView interface for choosing output mode
  void SetOutputMode(const int status);
  vtkGetMacro(OutputMode, int);

  // ParaView interface for selecting cell and point data fields
  int GetNumberOfCellArrays();
  const char* GetCellArrayName(const int index);
  int GetCellArrayStatus(const char* name);
  void SetCellArrayStatus(const char* name, const int status);

  int GetNumberOfPointArrays();
  const char* GetPointArrayName(const int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, const int status);

protected:

  vtkNetCDFLFRicReader();
  ~vtkNetCDFLFRicReader() override;

  // Return time steps from the input file and collect field (array) names
  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*) override;

  // Return grid and field data from the input file
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  // Build full VTK grid from UGRID description
  int CreateVTKGrid(netCDFLFRicFile& inputFile, vtkUnstructuredGrid *grid,
                    const size_t startLevel, const size_t numLevels,
                    const size_t numGhostsAbove, const size_t numGhostsBelow);

  // Points-only for W2 fields
  int CreateVTKPoints(netCDFLFRicFile& inputFile, vtkUnstructuredGrid *grid,
                      const size_t startLevel, const size_t numLevels,
                      const size_t numGhostsAbove, const size_t numGhostsBelow);

  // Read selected field data from netCDF file and add to the VTK grid
  int LoadFields(netCDFLFRicFile& inputFile, vtkUnstructuredGrid *grid,
                 const std::map<std::string, DataField> & fields,
                 const size_t timestep, const size_t startLevel,
                 const size_t numLevels);

private:

  vtkNetCDFLFRicReader(const vtkNetCDFLFRicReader&) = delete;
  void operator=(const vtkNetCDFLFRicReader&) = delete;

  char* FileName;
  int UseCartCoords, UseIndexAsVertCoord;
  double VerticalScale, VerticalBias;
  int OutputMode;
  std::map<std::string, DataField> CellFields;
  std::map<std::string, DataField> PointFields;
  std::vector<double> TimeSteps;

  UGRIDMeshDescription mesh2D;
  std::map<std::string, CFAxis> zAxes;
  CFAxis tAxis;

};
#endif
