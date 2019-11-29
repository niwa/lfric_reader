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

  bool HasDim(const std::string& dimName);

  size_t GetDimLen(const std::string& dimName);

  size_t GetVarNumDims(const std::string& varName);

  std::string GetVarDimName(const std::string& varName, const size_t dim);

  int GetAttInt(const std::string& varName, const std::string& attName);

  std::string GetAttText(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetAttTextSplit(const std::string& varName,
                                           const std::string& attName);

  bool HasVar(const std::string& varName);

  bool VarHasDim(const std::string& varName, const std::string& dimName);

  bool VarHasAtt(const std::string& varName, const std::string& attName);

  std::vector<std::string> GetVarNames();

  std::vector<double> GetVarDouble(const std::string& varName,
                                   const std::vector<size_t> start,
                                   const std::vector<size_t> count);

  std::vector<long long> GetVarLongLong(const std::string& varName,
                                        const std::vector<size_t> start,
                                        const std::vector<size_t> count);

private:

  netCDFLFRicFile(const netCDFLFRicFile&) = delete;
  void operator=(const netCDFLFRicFile&) = delete;

  int ncId;
  bool fileOpen;
  std::string FileName;
};

#endif
