/* hmmvp: Software to form and apply Hierarchical Matrices
 *   Version 1.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 *   https://pangea.stanford.edu/research/CDFM/software
 * hmmvp is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
*/

#ifndef INCLUDE_UTIL_DEBUG
#define INCLUDE_UTIL_DEBUG

#include <sstream>
#include "Matrix.hpp"

namespace util {

template<typename T> void WriteMatrix(const Matrix<T>& A, FILE* fid)
  throw(FileException)
{
  int64_t sz[2];
  for (size_t i = 1; i <= 2; ++i) sz[i-1] = A.Size(i);
  if (fwrite(sz, sizeof(int64_t), 2, fid) < 2)
    throw FileException("WriteMatrix: fwrite failed");
  int64_t n = sz[0]*sz[1];
  if (fwrite(A.GetPtr(), sizeof(T), n, fid) < (size_t) n)
    throw FileException("WriteMatrix: fwrite failed");
}

template<typename T> void WriteMatrix(const Matrix<T>& A, const std::string& fn)
  throw(FileException)
{
  FILE* fid;
  if (!(fid = fopen(fn.c_str(), "wb")))
    throw FileException(std::string("WriteMatrix: fopen failed on ") + fn);
  try { WriteMatrix(A, fid); }
  catch (FileException) {
    fclose(fid);
    throw FileException(std::string("WriteMatrix: fwrite failed on ") + fn); }
  fclose(fid);
}

template<typename T> void ReadMatrix(Matrix<T>& A, FILE* fid)
  throw(FileException,  OutOfMemoryException)
{
  int64_t sz[2];
  if (fread(sz, sizeof(int64_t), 2, fid) < 2)
    throw FileException("ReadMatrix: fread failed");
  A.Resize(sz[0], sz[1]);
  int64_t n = sz[0]*sz[1];
  int64_t k;
  if ((k = fread(A.GetPtr(), sizeof(T), n, fid)) < n)
    throw FileException("ReadMatrix: fread failed");
}
  
template<typename T> void ReadMatrix(Matrix<T>& A, const std::string& fn)
  throw(FileException,  OutOfMemoryException)
{
  FILE* fid;
  if (!(fid = fopen(fn.c_str(), "rb")))
    throw FileException(std::string("ReadMatrix: fopen failed on ") + fn);
  try { ReadMatrix(A, fid); }
  catch (FileException) {
    fclose(fid);
    throw FileException(std::string("ReadMatrix: fread failed on ") + fn); }
  fclose(fid);
}

template<typename T> void WriteMatrices
(const std::vector<Matrix<T> >& As, const std::string& fn)
  throw(FileException)
{
  FILE* fid;
  if (!(fid = fopen(fn.c_str(), "wb")))
    throw FileException(std::string("WriteMatrices: fopen failed on ") + fn);
  int64_t n = As.size();
  if (fwrite(&n, sizeof(int64_t), 1, fid) < 1) {
    fclose(fid);
    throw FileException(std::string("WriteMatrices: fwrite failed on ") + fn); }
  for (int64_t i = 0; i < n; ++i) {
    try { WriteMatrix(As[i], fid); }
    catch (FileException) {
      fclose(fid);
      throw FileException(std::string("WriteMatrices: fwrite failed on ") + fn);
    }}
  fclose(fid);
}

template<typename T> void ReadMatrices
(std::vector<Matrix<T> >& As, const std::string& fn)
  throw(FileException,  OutOfMemoryException)
{
  FILE* fid;
  if (!(fid = fopen(fn.c_str(), "rb")))
    throw FileException(std::string("ReadMatrices: fopen failed on ") + fn);
  int64_t n;
  if (fread(&n, sizeof(int64_t), 1, fid) < 1) {
    fclose(fid);
    throw FileException(std::string("ReadMatrices: fread failed on ") + fn); }
  As.resize(n);
  for (int64_t i = 0; i < n; ++i) {
    try { ReadMatrix(As[i], fid); }
    catch (FileException) {
      fclose(fid);
      throw FileException(std::string("ReadMatrices: fread failed on ") + fn);
    }}
  fclose(fid);
}

}

#endif
