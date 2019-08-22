#include "netCDFLFRicFile.h"
#include <vtk_netcdf.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <limits>

#define ncErrorMacro(e) if (e != NC_NOERR) {std::cerr << "NetCDF Error: " << nc_strerror(e) << "\n";}

//----------------------------------------------------------------------------

netCDFLFRicFile::netCDFLFRicFile(const char* fileName)
{
  this->FileName = fileName;
  const int err = nc_open(fileName, NC_NOWRITE, &this->ncId);
  ncErrorMacro(err);
  if (err != NC_NOERR)
  {
    this->fileOpen = false;
  }
  else
  {
    this->fileOpen = true;
  }
}

//----------------------------------------------------------------------------

netCDFLFRicFile::~netCDFLFRicFile()
{
  if (this->fileOpen)
  {
    ncErrorMacro(nc_close(this->ncId));
    this->fileOpen = false;
  }
  this->FileName.clear();
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::IsFileOpen()
{
  return this->fileOpen;
}

//----------------------------------------------------------------------------

const char* netCDFLFRicFile::GetFileName()
{
  return this->FileName.c_str();
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::HasDim(const std::string& dimName)
{
  int dimId;
  const int result = nc_inq_dimid(this->ncId, dimName.c_str(), &dimId);
  if (result == NC_NOERR)
  {
    return true;
  }
  else if (result == NC_EBADDIM)
  {
    return false;
  }
  else
  {
    ncErrorMacro(result);
  }
}

//----------------------------------------------------------------------------

size_t netCDFLFRicFile::GetDimLen(const std::string& dimName)
{
  int dimId;
  size_t dimLen;
  ncErrorMacro(nc_inq_dimid(this->ncId, dimName.c_str(), &dimId));
  ncErrorMacro(nc_inq_dimlen(this->ncId, dimId, &dimLen));
  return dimLen;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetVarDimName(const std::string& varName,
                                           const size_t dim)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));

  std::string dimNameStr;
  dimNameStr.clear();

  // Return empty string if out of range
  if (dim >= 0 and dim < numDims)
  {
    std::vector<int> dimIds;
    dimIds.resize(numDims);
    ncErrorMacro(nc_inq_vardimid(this->ncId, varId, dimIds.data()));

    char *dimName = new char[NC_MAX_NAME+1];
    ncErrorMacro(nc_inq_dimname(this->ncId, dimIds[dim], dimName));

    // NetCDF docs say that dimension names are guaranteed to be null-terminated
    dimNameStr = dimName;

    delete[] dimName;
  }

  return dimNameStr;
}

//----------------------------------------------------------------------------

std::string netCDFLFRicFile::GetAttText(const std::string& varName, const std::string& attName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  size_t attTextLen;
  ncErrorMacro(nc_inq_attlen(this->ncId, varId, attName.c_str(), &attTextLen));

  char *attText = new char[attTextLen+1];
  ncErrorMacro(nc_get_att_text(this->ncId, varId, attName.c_str(), attText));

  // We cannot assume that attText is always null-terminated, so copy only attTextLen bytes
  std::string attTextStr(attText, attTextLen);
  delete[] attText;

  return attTextStr;
}

//----------------------------------------------------------------------------

std::vector<std::string> netCDFLFRicFile::GetAttTextSplit(const std::string& varName,
                                                          const std::string& attName)
{
  std::string attText = this->GetAttText(varName, attName);
  std::istringstream attStream(attText);
  std::vector<std::string> textSplit(std::istream_iterator<std::string>{attStream},
                                     std::istream_iterator<std::string>());
  return textSplit;
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::VarHasDim(const std::string& varName, const std::string& dimName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));

  std::vector<int> dimIds;
  dimIds.resize(numDims);
  ncErrorMacro(nc_inq_vardimid(this->ncId, varId, dimIds.data()));

  int dimIdOfInterest;
  ncErrorMacro(nc_inq_dimid(this->ncId, dimName.c_str(), &dimIdOfInterest));

  bool dimFound = false;
  for (int const dimId : dimIds)
  {
    if (dimId == dimIdOfInterest) dimFound = true;
  }

  return dimFound;
}

//----------------------------------------------------------------------------

bool netCDFLFRicFile::VarHasAtt(const std::string& varName, const std::string& attName)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  int attId;
  const int result = nc_inq_attid(this->ncId, varId, attName.c_str(), &attId);
  if (result == NC_NOERR)
  {
    return true;
  }
  else if (result == NC_ENOTATT)
  {
    return false;
  }
  else
  {
    ncErrorMacro(result);
  }
}

//----------------------------------------------------------------------------

std::vector<std::string> netCDFLFRicFile::GetVarNames()
{
  std::vector<std::string> varNames;
  int numVars;
  ncErrorMacro(nc_inq_nvars(this->ncId, &numVars));
  varNames.resize(numVars);

  // NetCDF variables can use up to NC_MAX_NAME bytes
  char *varName = new char[NC_MAX_NAME+1];
  for (int iVar = 0; iVar < numVars; iVar++)
  {
    ncErrorMacro(nc_inq_varname(this->ncId, iVar, varName));
    // NetCDF docs say that variable names are guaranteed to be null-terminated
    varNames[iVar] = varName;
  }
  delete[] varName;
  return varNames;
}

//----------------------------------------------------------------------------

std::vector<double> netCDFLFRicFile::GetVarDouble(
                    const std::string& varName,
                    const std::vector<size_t> start,
                    const std::vector<size_t> count)
{
  // Find variable by name
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  // Compute total number of elements to read
  size_t size = 1;
  for (size_t const n : count)
  {
    size *= n;
  }

  std::vector<double> varData;
  varData.resize(size);

  // Check if number of dimensions matches netCDF variable
  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));
  if (numDims != start.size() or numDims != count.size())
  {
    std::cerr << "netCDFLFRicFile::GetVarDouble: number of dimensions does not match netCDF variable.\n";
    std::fill(varData.begin(), varData.end(), std::numeric_limits<double>::quiet_NaN());
  }
  else
  {
    // NetCDF will automatically convert non-double numeric data into double
    // This function will also check index ranges automatically
    ncErrorMacro(nc_get_vara_double(this->ncId, varId, start.data(),
                                    count.data(), varData.data()));
  }

  return varData;
}

//----------------------------------------------------------------------------

std::vector<long long> netCDFLFRicFile::GetVarLongLong(
                       const std::string& varName,
                       const std::vector<size_t> start,
                       const std::vector<size_t> count)
{
  int varId;
  ncErrorMacro(nc_inq_varid(this->ncId, varName.c_str(), &varId));

  size_t size = 1;
  for (size_t const n : count)
  {
    size *= n;
  }

  std::vector<long long> varData;
  varData.resize(size);

  // Check if number of dimensions matches netCDF variable
  int numDims;
  ncErrorMacro(nc_inq_varndims(this->ncId, varId, &numDims));
  if (numDims != start.size() or numDims != count.size())
  {
    std::cerr << "netCDFLFRicFile::GetVarLongLong: number of dimensions does not match netCDF variable.\n";
    std::fill(varData.begin(), varData.end(), 0);
  }
  else
  {
    // netCDF will automatically convert non-double numeric data into long long
    ncErrorMacro(nc_get_vara_longlong(this->ncId, varId, start.data(),
                                      count.data(), varData.data()));
  }

  return varData;
}
