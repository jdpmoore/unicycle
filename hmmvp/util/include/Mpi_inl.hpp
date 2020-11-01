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

#ifdef UTIL_MPI
#else
#endif

#include <string.h>
#include <iostream>
#include "util/include/CodeAnalysis.hpp"

namespace util {
namespace mpi {

template<typename T> MPI_Datatype GetType () {
#ifdef UTIL_MPI
  T v = (T) 1.01;
  if (v == 1) {
    switch (sizeof(T)) {
    case sizeof(char):  return MPI_CHAR;
    case sizeof(short): return MPI_SHORT;
    case sizeof(int):   return MPI_INT;
    case sizeof(long):  return MPI_LONG;
    default: ;
    }
  } else {
    switch (sizeof(T)) {
    case sizeof(float):  return MPI_FLOAT;
    case sizeof(double): return MPI_DOUBLE;
    default: ;
    }
  }
  cerr << "Pid " << mpi::Pid() << " failed in mpi::GetType(); " <<
    "v = " << (double) v << ", sizeof(T) = " << sizeof(T) << endl;
  return 0;
#else
  return 1;
#endif
}

template<typename T>
int Reduce (const T* sendbuf, T* rcvbuf, int count, MPI_Op op, int root) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Reduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op, root,
                    MPI_COMM_WORLD);
#else
  if (rcvbuf) memcpy(rcvbuf, sendbuf, count*sizeof(T));
  return 0;
#endif
}

template<typename T>
int Allreduce (const T* sendbuf, T* rcvbuf, int count, MPI_Op op) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Allreduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op,
                       MPI_COMM_WORLD);
#else
  if (rcvbuf) memcpy(rcvbuf, sendbuf, count*sizeof(T));
  return 0;
#endif
}

template<typename T>
int Bcast (T* buf, int count, int root) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Bcast(buf, count, dt, root, MPI_COMM_WORLD);
#else
  return 0;
#endif
}

template<typename T>
int Gather (const T* sendbuf, int sendcnt, T* rcvbuf, int rcvcnt, int root) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Gather(const_cast<T*>(sendbuf), sendcnt, dt, rcvbuf,
                    rcvcnt, dt, root, MPI_COMM_WORLD);
#else
  if (rcvbuf) memcpy(rcvbuf, sendbuf, sendcnt*sizeof(T));
  return 0;
#endif      
}

template<typename T>
int Allgather (const T* sendbuf, int sendcnt, T* rcvbuf, int rcvcnt) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Allgather(const_cast<T*>(sendbuf), sendcnt, dt, rcvbuf,
                       rcvcnt, dt, MPI_COMM_WORLD);
#else
  if (rcvbuf) memcpy(rcvbuf, sendbuf, sendcnt*sizeof(T));
  return 0;
#endif
}

template<typename T>
int Scatter (const T* sendbuf, int sendcnt, T* rcvbuf, int rcvcnt, int root) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Scatter(const_cast<T*>(sendbuf), sendcnt, dt, rcvbuf,
                     rcvcnt, dt, root, MPI_COMM_WORLD);
#else
  if (rcvbuf) memcpy(rcvbuf, sendbuf, sendcnt*sizeof(T));
  return 0;
#endif
}

template<typename T>
int Isend (const T* buf, int count, int dest, int tag, MPI_Request* ireq) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? ireq : &ureq;
  int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag,
                      MPI_COMM_WORLD, req);
  if (!ireq) MPI_Request_free(req);
  return ret;
#else
  return 0;
#endif
}

template<typename T>
int Send (const T* buf, int count, int dest, int tag) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  return MPI_Send(const_cast<T*>(buf), count, dt, dest, tag,
                  MPI_COMM_WORLD);
#else
  return 0;
#endif
}

template<typename T>
int Irecv (T* buf, int count, int src, int tag, MPI_Request* ireq) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? ireq : &ureq;
  int ret = MPI_Irecv(buf, count, dt, src, tag, MPI_COMM_WORLD, req);
  if (!ireq) MPI_Request_free(req);
  return ret;
#else
  return 0;
#endif
}

template<typename T>
int Recv (T* buf, int count, int src, int tag) {
#ifdef UTIL_MPI
  MPI_Datatype dt = GetType<T>();
  MPI_Status stat;
  return MPI_Recv(buf, count, dt, src, tag, MPI_COMM_WORLD, &stat);
#else
  return 0;
#endif
}

template<typename T>
int ArraySegmenter::Gather (const T* sndbuf, T* rcvbuf, int root, int factor)
  const
{
  int cnt = factor*GetN();
#ifdef UTIL_MPI
  if (_all_equal_size)
    return mpi::Gather(const_cast<T*>(sndbuf), cnt, rcvbuf, cnt, root);
  else {
    int nproc = _bds.size() - 1;
    int tag = 101;
    int ret;
    if (mpi::AmRoot()) {
      _reqs.resize(nproc - 1);
      _stats.resize(nproc - 1);
      for (int i = 0, os = 0, ir = 0; i < nproc; i++) {
        int n = factor * (_bds[i+1] - _bds[i]);
        if (i == root) {
          // No need for MPI.
          memcpy(rcvbuf + os, sndbuf, n * sizeof(T));
        } else {
          // Receive async'ly from everyone else.
          if ((ret = mpi::Irecv(rcvbuf + os, n, i, tag, &_reqs[ir])))
            return ret;
          ir++;
        }
        os += n;
      }
      // Wait for all my async receives to finish.
      if ((ret = mpi::Waitall(nproc - 1, &_reqs[0], &_stats[0])))
        return ret;
      // Is this necessary? I don't think so, but I want to make sure I don't
      // get leaks. [Later:] Actually, it seems to result in MPI_REQUEST_ERR.
      //for (int i = 0; i < nproc - 1; i++) mpi::Request_free(&_reqs[i]);
    } else {
      // Send sync'ly to root.
      int n = factor * this->GetN();
      if ((ret = mpi::Send(sndbuf, n, root, tag))) return ret;
    }
    return 0;
  }
#else
  if (rcvbuf) memcpy(rcvbuf, sndbuf, cnt*sizeof(T));
  return 0;
#endif
}

template<typename T>
int ArraySegmenter::Allgather (const T* sndbuf, T* rcvbuf, int factor) const {
  int cnt = factor*GetN();
#ifdef UTIL_MPI
  if (_all_equal_size)
    return mpi::Allgather(sndbuf, cnt, rcvbuf, cnt);
  else {
    int ret;
    int tag = 102;
    int pid = mpi::Pid(), nproc = mpi::GetNproc();
    int nreqs = 2*(nproc - 1);
    _reqs.resize(nreqs);
    _stats.resize(nreqs);
    int ir = 0;
    // Send async'ly to everyone else.
    int n = factor * (_bds[pid + 1] - _bds[pid]);
    for (int i = 0; i < nproc; i++)
      if (i != pid) {
        if ((ret = mpi::Isend(sndbuf, n, i, tag, &_reqs[ir]))) return ret;
        ir++;
      }
    // Receive async'ly from everyone.
    for (int i = 0, os = 0; i < nproc; i++) {
      int n = factor * (_bds[i+1] - _bds[i]);
      if (i == pid) {
        memcpy(rcvbuf + os, sndbuf, n * sizeof(T));
      } else {
        if ((ret = mpi::Irecv(rcvbuf + os, n, i, tag, &_reqs[ir])))
          return ret;
        ir++;
      }
      os += n;
    }
    // Wait on async sends and receives.
    if ((ret = mpi::Waitall(nreqs, &_reqs[0], &_stats[0]))) return ret;
    for (int i = 0; i < nreqs; i++) mpi::Request_free(&_reqs[i]);
    return 0;
  }
#else
  if (rcvbuf) memcpy(rcvbuf, sndbuf, cnt*sizeof(T));
  return 0;
#endif
}

template<typename T>
int ArraySegmenter::Scatter (const T* sndbuf, T* rcvbuf, int root,
                             int factor) const {
  int cnt = factor*GetN();
#ifdef UTIL_MPI
  if (_all_equal_size)
    return mpi::Scatter(sndbuf, cnt, rcvbuf, cnt, root);
  else {
    int ret;
    int tag = 103;
    int nproc = mpi::GetNproc();
    if (mpi::AmRoot()) {
      _reqs.resize(nproc - 1);
      _stats.resize(nproc - 1);
      for (int i = 0, os = 0, ir = 0; i < nproc; i++) {
        int n = factor * (_bds[i+1] - _bds[i]);
        if (i == root) {
          // No need for MPI.
          memcpy(rcvbuf, sndbuf + os, n * sizeof(T));
        } else {
          // Send async'ly to everyone else.
          if ((ret = mpi::Isend(sndbuf + os, n, i, tag, &_reqs[ir])))
            return ret;
          ir++;
        }
        os += n;
      }
      // Wait for all my async sends to finish.
      if ((ret = mpi::Waitall(nproc - 1, &_reqs[0], &_stats[0])))
        return ret;
      //todo I still can't figure out when it's ok/necessary to call this when
      // using Isend.
      //for (int i = 0; i < nproc - 1; i++) mpi::Request_free(&_reqs[i]);
    } else {
      // Receive sync'ly from root.
      int n = factor * this->GetN();
      if ((ret = mpi::Recv(rcvbuf, n, root, tag))) return ret;
    }
    return 0;
  }
#else
  if (rcvbuf) memcpy(rcvbuf, sndbuf, cnt*sizeof(T));
  return 0;
#endif
}

inline ByteBufferWriter::ByteBufferWriter(size_t init_sz)
{ if (init_sz > 0) _buf.reserve(init_sz); }

inline void ByteBufferWriter::Reset() { _buf.clear(); }

template<typename T>
inline void ByteBufferWriter::Write(const T* data, size_t n)
{
  size_t nbytes = n * sizeof(T);
  _buf.insert(_buf.end(), (const char*) data, (const char*) data + nbytes);
}

template<typename T>
inline void ByteBufferWriter::WriteScalar(T data) { Write(&data, 1); }

template<typename T>
inline T* ByteBufferWriter::CopyTo (size_t n) {
  size_t nbytes = n * sizeof(T);
  size_t i = _buf.size();
  _buf.insert(_buf.end(), nbytes, 0);
  return (T*) &_buf[i];
}

template<typename T>
inline void ByteBufferWriter::Resize (size_t n) { _buf.resize(n); }

inline size_t ByteBufferWriter::Size () const { return _buf.size(); }

inline const char* ByteBufferWriter::GetPtr () const { return &_buf[0]; }

inline void ByteBufferWriter::Bcast () {
  size_t sz = _buf.size();
  mpi::Bcast(&sz, 1);
  mpi::Bcast(&_buf[0], sz);
}

inline void ByteBufferWriter::Isend (int dest) {
  int tag = 0;
  size_t sz = _buf.size();
  mpi::Isend(&sz, 1, dest, tag);
  if (sz > 0) mpi::Isend(&_buf[0], sz, dest, tag);
}

inline void ByteBufferWriter::Send (int dest) {
  int tag = 0;
  size_t sz = _buf.size();
  mpi::Send(&sz, 1, dest, tag);
  if (sz > 0) mpi::Send(&_buf[0], sz, dest, tag);
}

inline ByteBufferReader::ByteBufferReader (const char* buf)
  : _buf(buf), _i(0) {}

inline const char* ByteBufferReader::CurrentPosition () const
{ return _buf + _i; }

template<typename T>
inline void ByteBufferReader::Read (T* data, size_t n) {
  size_t nbytes = n * sizeof(T);
  memcpy((char*) data, _buf + _i, nbytes);
  _i += nbytes;
}

template<typename T>
inline T ByteBufferReader::ViewScalar (size_t pos) const {
  T s;
  memcpy(&s, _buf + pos, sizeof(T));
  return s;
}

template<typename T>
inline T ByteBufferReader::ViewScalarAndAdvance () {
  T s = ViewScalar<T>(_i);
  _i += sizeof(T);
  return s;
}

inline ByteBufferReaderMpi::ByteBufferReaderMpi()
  : ByteBufferReader(NULL)
{}

inline void ByteBufferReaderMpi::Bcast () {
  size_t sz;
  mpi::Bcast(&sz, 1);
  _data.resize(sz);
  _buf = &_data[0];
  _i = 0;
  if (sz > 0) mpi::Bcast(&_data[0], sz);
}

inline void ByteBufferReaderMpi::Recv (int src) {
  int tag = 0;
  size_t sz;
  mpi::Recv(&sz, 1, src, tag);
  _data.resize(sz);
  _buf = &_data[0];
  _i = 0;
  if (sz > 0) mpi::Recv(&_data[0], sz, src, tag);
}

inline size_t ByteBufferReaderMpi::Size() const { return _data.size(); }

}}
