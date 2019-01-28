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
#include <vtk_netcdf.h>

#define DEBUG 1

class VTKIONETCDF_EXPORT vtkNetCDFLFRicReader : public vtkUnstructuredGridAlgorithm {
public:
  vtkTypeMacro(vtkNetCDFLFRicReader,vtkUnstructuredGridAlgorithm);
  static vtkNetCDFLFRicReader *New();
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Specify name of input data file
  vtkSetStringMacro(FileName); // Defines SetFileName()
  vtkGetStringMacro(FileName); // Defines GetFileName()

protected:

  vtkNetCDFLFRicReader();
  ~vtkNetCDFLFRicReader() override;

  // Return time steps from the input file
  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*) override;

  // Return grid and field data from the input file
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  // Build VTK grid from UGRID description
  int CreateVTKGrid(const int ncid, vtkUnstructuredGrid *grid);

  // Read field data from netCDF file and add to the VTK grid
  int LoadFields(const int ncid, vtkUnstructuredGrid *grid);

  // Utility functions
  size_t getNCDim(const int ncid, const char * dimname);
  // FIXME: use rvalue here
  int getNCVar(const int ncid, const char * varname, const nc_type expect_vartype,
               const int expect_ndims, const size_t start [], const size_t count [], void * buffer);

  char *FileName;

private:

  vtkNetCDFLFRicReader(const vtkNetCDFLFRicReader&) = delete;
  void operator=(const vtkNetCDFLFRicReader&) = delete;

  double * TimeSteps;
  size_t NumberOfTimeSteps, NumberOfLevels, NumberOfFaces2D;

};
#endif
