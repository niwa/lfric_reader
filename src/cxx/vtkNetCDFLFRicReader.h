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
#include <map>
#include <string>

#define DEBUG 1

class VTKIONETCDF_EXPORT vtkNetCDFLFRicReader : public vtkUnstructuredGridAlgorithm
{

public:

  vtkTypeMacro(vtkNetCDFLFRicReader,vtkUnstructuredGridAlgorithm)
  static vtkNetCDFLFRicReader *New();
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Specify name of input data file
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // ParaView interface for switching coordinates
  void SetUseCartCoords(const int status);
  vtkGetMacro(UseCartCoords, int);

  // ParaView interface for selecting data fields
  int GetNumberOfCellArrays();
  const char* GetCellArrayName(const int index);
  int GetCellArrayStatus(const char* name);
  void SetCellArrayStatus(const char* name, const int status);

protected:

  vtkNetCDFLFRicReader();
  ~vtkNetCDFLFRicReader() override;

  // Data can be defined on these 4 different mesh types
  enum mesh_types {nodal, full_level_face, half_level_face, half_level_edge};

  // Return time steps from the input file and collect field (array) names
  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*) override;

  // Return grid and field data from the input file
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  // Build VTK grid from UGRID description
  int CreateVTKGrid(const int ncid, vtkUnstructuredGrid *grid);

  // Read selected field data from netCDF file and add to the VTK grid
  int LoadFields(const int ncid, vtkUnstructuredGrid *grid, const size_t timestep);

  // NetCDF utility functions
  size_t getNCDim(const int ncid, const char * dimname);
  std::vector<std::string> getNCVarNames(const int ncid);
  std::vector<double> getNCVarDouble(const int ncid, const char * varname,
                                     const std::initializer_list<size_t> start,
                                     const std::initializer_list<size_t> count);
  std::vector<unsigned long long> getNCVarULongLong(const int ncid, const char * varname,
                                                    const std::initializer_list<size_t> start,
                                                    const std::initializer_list<size_t> count);

  // Transforms periodic grid into non-periodic grid by replicating vertices (points)
  void mirror_points(vtkSmartPointer<vtkUnstructuredGrid> grid);

private:

  vtkNetCDFLFRicReader(const vtkNetCDFLFRicReader&) = delete;
  void operator=(const vtkNetCDFLFRicReader&) = delete;

  char *FileName;
  int UseCartCoords;
  std::map<std::string,bool> Fields;
  std::vector<double> TimeSteps;
  size_t NumberOfLevels, NumberOfFaces2D, NumberOfEdges2D;

};
#endif
