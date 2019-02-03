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

#include <vtkIONetCDFModule.h> // For export macro
#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtk_netcdf.h>

#include <vector>

// Set automatically with PV plugin debug?
#define DEBUG 1

class VTKIONETCDF_EXPORT vtkNetCDFLFRicReader : public vtkUnstructuredGridAlgorithm {

public:

  vtkTypeMacro(vtkNetCDFLFRicReader,vtkUnstructuredGridAlgorithm);
  static vtkNetCDFLFRicReader *New();
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Specify name of input data file
  vtkSetStringMacro(FileName); // Defines SetFileName()
  vtkGetStringMacro(FileName); // Defines GetFileName()

  void SetUseCartCoords(bool SetCartCoords);
  vtkGetMacro(UseCartCoords, bool);

protected:

  vtkNetCDFLFRicReader();
  ~vtkNetCDFLFRicReader() override;

  // Data is defined on 4 different grids
  // nodal - full level face - half level face - half level edge
  enum mesh_types {nodal, full_level_face, half_level_face, half_level_edge};

  // Return time steps from the input file
  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*) override;

  // Return grid and field data from the input file
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  // Build VTK grid from UGRID description
  int CreateVTKGrid(const int ncid, vtkUnstructuredGrid *grid);

  // Read field data from netCDF file and add to the VTK grid
  int LoadFields(const int ncid, vtkUnstructuredGrid *grid, size_t timestep);

  // Utility functions
  size_t getNCDim(const int ncid, const char * dimname);

  std::vector<double> getNCVarDouble(const int ncid, const char * varname, const std::initializer_list<size_t> start, const std::initializer_list<size_t> count);
  std::vector<unsigned long long> getNCVarULongLong(const int ncid, const char * varname, const std::initializer_list<size_t> start, const std::initializer_list<size_t> count);
  void mirror_points(vtkSmartPointer<vtkUnstructuredGrid> grid);

  char *FileName;

  bool UseCartCoords;

private:

  vtkNetCDFLFRicReader(const vtkNetCDFLFRicReader&) = delete;
  void operator=(const vtkNetCDFLFRicReader&) = delete;

  std::vector<double> TimeSteps;
  size_t NumberOfLevels, NumberOfFaces2D, NumberOfEdges2D;

};
#endif
