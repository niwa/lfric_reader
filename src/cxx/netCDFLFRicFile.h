/*
 * @class   netCDFLFRicFile
 * @brief   Utility class for netCDF file handling
 *
 * Reads metadata and data from a netCDF file using
 * convenient utility functions.
*/

#ifndef netCDFLFRicFile_h
#define netCDFLFRicFile_h
#include <vector>
#include <string>

class netCDFLFRicFile
{

public:

  netCDFLFRicFile(const char* fileName);
  ~netCDFLFRicFile();

  const char* GetFileName();

  size_t GetDimLen(const std::string& dimName);

  std::string GetAttText(const std::string& varName, const std::string& attName);

  bool VarHasDim(const std::string& varName, const std::string& dimName);

  std::vector<std::string> GetVarNames();

  std::vector<double> GetVarDouble(const std::string& varName,
                                   const std::initializer_list<size_t> start,
                                   const std::initializer_list<size_t> count);

  std::vector<unsigned long long> GetVarULongLong(const std::string& varName,
                                                  const std::initializer_list<size_t> start,
                                                  const std::initializer_list<size_t> count);

private:

  netCDFLFRicFile(const netCDFLFRicFile&) = delete;
  void operator=(const netCDFLFRicFile&) = delete;

  int ncId;
  bool fileOpen;
  std::string FileName;
};

#endif
