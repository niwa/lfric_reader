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

  bool IsFileOpen();

  const char* GetFileName();

  size_t GetDimLen(const std::string& dimName);

  std::string GetVarDimName(const std::string& varName, const size_t dim);

  std::string GetAttText(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetAttTextSplit(const std::string& varName,
                                           const std::string& attName);

  bool VarHasDim(const std::string& varName, const std::string& dimName);

  bool VarHasAtt(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetVarNames();

  std::vector<double> GetVarDouble(const std::string& varName,
                                   const std::initializer_list<size_t> start,
                                   const std::initializer_list<size_t> count);

  std::vector<long long> GetVarLongLong(const std::string& varName,
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
