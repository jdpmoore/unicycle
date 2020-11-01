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

#ifndef INCLUDE_UTIL_MATRIX
#define INCLUDE_UTIL_MATRIX

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "util/include/Exception.hpp"
#include "util/include/IO.hpp"

namespace util {
template<typename T, size_t idx_start=1>
class Matrix {
private:
  T* _data;
  size_t _m, _n;
  bool _am_rsbl; // Am responsible to delete _data.
  // NOT YET: T _stat_data[static_data_size];

public:
  Matrix () : _data(0), _m(0), _n(0), _am_rsbl(false) { Resize(0, 0); }
  explicit Matrix (size_t m) : _data(0), _m(0), _n(0), _am_rsbl(false) {
    Resize(m, 1); }
  Matrix (size_t m, size_t n) : _data(0), _m(0), _n(0), _am_rsbl(false) {
    Resize(m, n); }
  Matrix (size_t m, size_t n, T* data)
    : _data(0), _m(0), _n(0), _am_rsbl(false) { SetPtr(m, n, data); }
  Matrix (size_t m, T* data)
    : _data(0), _m(0), _n(0), _am_rsbl(false) { SetPtr(m, data); }
  ~Matrix() { if (_am_rsbl && _data) free(_data); }

  void Serialize(FILE* fid) const throw (FileException) {
    util::write(&_m, 1, fid);
    util::write(&_n, 1, fid);
    util::write(_data, _m*_n, fid);
  }
  void Deserialize(FILE* fid) throw (FileException) {
    size_t m, n;
    util::read(&m, 1, fid);
    util::read(&n, 1, fid);
    Resize(m, n);
    util::read(_data, _m*_n, fid);
  }

  Matrix (const Matrix<T>& A) throw (OutOfMemoryException)
    : _m(0), _n(0), _am_rsbl(false) {
    Resize(A._m, A._n);
    memcpy(_data, A._data, _m*_n*sizeof(T));
  }
  template<typename T1> Matrix (const Matrix<T1>& A)
    throw (OutOfMemoryException) : _m(0), _n(0), _am_rsbl(false) {
    Resize(A.Size(1), A.Size(2));
    size_t sz = Size();
    for (size_t i = 0; i < sz; ++i) _data[i] = (T)A(i+1);
  }

  Matrix<T>& operator= (const Matrix<T>& A) throw (OutOfMemoryException) {
    if (this != &A) {
      Resize(A._m, A._n);
      memcpy(_data, A._data, _m*_n*sizeof(T)); }
    return *this;
  }
  template<typename T1> Matrix<T>& operator= (const Matrix<T1>& A)
    throw (OutOfMemoryException) {
    if (this != (const Matrix<T>*) &A) {
      Resize(A.Size(1), A.Size(2));
      size_t sz = Size();
      for (size_t i = 0; i < sz; ++i) _data[i] = (T)A(i+1); }
    return *this;
  }

  // Get size of row (dim=1),  col (dim=2),  or row*col (dim=0)
  size_t Size (size_t dim = 0) const {
    switch (dim) {
    case 1: return _m;
    case 2: return _n;
    default: return _m*_n; }
  }

  Matrix<T>& Resize (size_t m) throw (OutOfMemoryException) {
    return Resize(m, 1); }
  Matrix<T>& Reshape (size_t m) throw (OutOfMemoryException) {
    return Reshape(m, 1); }

  // Resize matrix. Memory now handled by Matrix.
  Matrix<T>& Resize (size_t m, size_t n) throw (OutOfMemoryException) {
    if (_m*_n != m*n) {
      if (_am_rsbl) {
        if (_data) free(_data);
        if (m == 0 || n == 0) _am_rsbl = false;
        else _data = (T*)malloc(m*n*sizeof(T));
        if (!_data) throw OutOfMemoryException();
      } else if (m > 0 && n > 0) {
        _data = (T*)malloc(m*n*sizeof(T));
        if (!_data) throw OutOfMemoryException();
        _am_rsbl = true;
      }
    }
    _m = m; _n = n;
    if (_m == 0 || _n == 0) _m = _n = 0;
    return *this;
  }
  
  // Same as Resize, but copies old data over. The extra is set to zero.
  Matrix<T>& Reshape (size_t m, size_t n) throw (OutOfMemoryException) {
    if (_m*_n != m*n) {
      if (m == 0 || n == 0) {
        if (_data) free(_data);
        _am_rsbl = false;
      } else {
        if (_am_rsbl) {
          _data = (T*)realloc(_data, m*n*sizeof(T));
          if (!_data) throw OutOfMemoryException();
        } else if (m > 0 && n > 0) {
          T* data = (T*)malloc(m*n*sizeof(T));
          if (!data) throw OutOfMemoryException();
          _am_rsbl = true;
          if (_m > 0 && _n > 0) memcpy(data, _data, _m*_n*sizeof(T));
          _data = data; }
        // Zero the remaining memory
        if (m*n > _m*_n) memset(_data + _m*_n, 0, (m*n - _m*_n)*sizeof(T));
      }
    }
    _m = m; _n = n;
    if (_m == 0 || _n == 0) _m = _n = 0;
    return *this;
  }
  
  Matrix<T>& SetPtr (T* data) { return SetPtr(_m, _n, data); }
  Matrix<T>& SetPtr (size_t m, T* data) { return SetPtr(m, 1, data); }
  Matrix<T>& SetPtr (size_t m, size_t n, T* data) {
    if (_am_rsbl) {
      if (_data) free(_data);
      _am_rsbl = false;
    }
    _data = data;
    _m = m; _n = n;
    return *this;
  }

  T* GetPtr () { return _data; }
  const T* GetPtr () const { return _data; }

  T& operator() (size_t i) {
    if (!(i >= idx_start && i < _m*_n + idx_start)) {
      int* bogus = 0; *bogus = 0;
    }
    assert(i >= idx_start);
    assert(i < _m*_n + idx_start);
    return _data[i-idx_start];
  }
  const T& operator() (size_t i) const {
    assert(i >= idx_start && i < _m*_n + idx_start);
    return _data[i-idx_start];
  }
  T& operator() (size_t i, size_t j) {
    assert(i >= idx_start && j >= idx_start
           && i < _m + idx_start && j < _n + idx_start);
    return _data[(j - idx_start)*_m + i - idx_start];
  }
  const T& operator() (size_t i, size_t j) const {
    assert(i >= idx_start && j >= idx_start
           && i < _m + idx_start && j < _n + idx_start);
    return _data[(j - idx_start)*_m + i - idx_start];
  }

  Matrix<T>& Bzero () { memset(_data, 0, _m*_n*sizeof(T)); return *this; }

  Matrix<T>& Set (T v) {
    for (size_t i = 0; i < _m*_n; ++i) _data[i] = v;
    return *this;
  }

  Matrix<T>& Zero() { return Set((T)0); }

  Matrix<T>& operator*= (T s) {
    for (size_t i = 0; i < _m*_n; ++i) _data[i] *= s;
    return *this;
  }

  Matrix<T>& operator+= (T s) {
    for (size_t i = 0; i < _m*_n; ++i) _data[i] += s;
    return *this;
  }
};

template<typename T1, typename T2>
Matrix<T2>& Transpose (const Matrix<T1>& A, Matrix<T2>& At) {
  At.Resize(A.Size(2), A.Size(1));
  for (size_t i = 1; i <= A.Size(1); ++i)
    for (size_t j = 1; j <= A.Size(2); ++j)
      At(j,i) = A(i,j);
  return At;
}

template<typename T>
void mvp (const T* A, size_t m, size_t n, const T* x, T* y) {
  size_t os = 0;
  memset(y, 0, m*sizeof(T));
  for (size_t j = 0; j < n; j++) {
    double xj = x[j];
    for (size_t i = 0; i < m; i++)
      y[i] += A[os + i]*xj;
    os += m;
  }
}

template<typename T>
Matrix<T>& mvp (const Matrix<T>& A, const Matrix<T>& x, Matrix<T>& y) {
  size_t m = A.Size(1), n = A.Size(2);
  y.Resize(m);
  mvp(A.GetPtr(), m, n, x.GetPtr(), y.GetPtr());
  return y;
}

template<typename T>
Matrix<T>& mmp (const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C) {
  size_t m = A.Size(1), n = A.Size(2), Bm = B.Size(1), Bn = B.Size(2);
  size_t Bos = 0, Cos = 0;
  C.Resize(m, Bn);
  for (size_t k = 0; k < Bn; k++) {
    mvp(A.GetPtr(), m, n, B.GetPtr() + Bos, C.GetPtr() + Cos);
    Bos += Bm;
    Cos += m;
  }
  return C;
}

template<typename T>
std::ostream& operator<< (std::ostream& os, const Matrix<T>& A) {
  size_t m,n;
  m = A.Size(1);
  n = A.Size(2);
  os << "[";
  if (n == 1)
    for (size_t r = 1; r <= m; ++r) os << A(r) << " ";
  else
    for (size_t r = 1; r <= m; ++r) {
      if (r>1) os << "\n";
      for (size_t c = 1; c < n; ++c) os << A(r,c) << " ";
      os << A(r,n) << ";";
    }
  os << "]";
  return os;
}
}

#endif
