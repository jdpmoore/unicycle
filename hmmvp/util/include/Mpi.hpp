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

#ifndef INCLUDE_UTIL_MPI
#define INCLUDE_UTIL_MPI

// Simplified set of MPI routines. If the symbol UTIL_MPI is not defined, then
// the program can be compiled as a serial, non-MPI-dependent program.

#include <vector>

#ifdef UTIL_MPI
#  include <mpi.h>
#endif
#include "util/include/OpenMP.hpp"

//todo Since I'm using Mpi.?pp across multiple projects, I might take away this
// outer namespace. I could still protect the name "mpi" by including a wrapper
// header to this file.
namespace util {
  using namespace std;

#ifndef UTIL_MPI
    typedef int MPI_Datatype;
    typedef int MPI_Op;
    typedef int MPI_Request;
    typedef int MPI_Status;
# define MPI_SUM 1
# define MPI_MAX 1
# define MPI_LAND 1
# define MPI_SUCCESS 0
#endif

  namespace mpi {

    void Init(int argc, char** argv);
    void Finalize();
    int  GetNproc();
    int  Pid();
    void Barrier();
    int  Root();
    bool AmRoot();

    // For convenience, provide an AmRoot-like functions that handles both MPI
    // and OpenMP.
    bool AmTPid0();

    // Report whether the root process says we're ok. Input is significant only
    // for root.
    bool IsTrue(bool p = true);
    // I need to make this a macro, I think, because it's the only way to get
    // the short-circuit && to work properly. If I make it a function, then "b"
    // is evaluated before mpi_IsTrue would be called. I don't want that.
    #define mpi_IsTrue(b)                                               \
      (( mpi::AmRoot() && mpi::IsTrue(b)) ||                            \
       (!mpi::AmRoot() && mpi::IsTrue( )))

    // True only if all procs say so.
    bool AllOk(bool p = true);
    
    template<typename T>
    int Reduce(const T* sendbuf, T* rcvbuf, int count, MPI_Op op, int root = 0);
    template<typename T>
    int Allreduce(const T* sendbuf, T* rcvbuf, int count, MPI_Op op);
    template<typename T>
    int Bcast(T* buf, int count, int root = 0);
    template<typename T>
    int Gather(const T* sendbuf, int sendcnt, T* rcvbuf, int rcvcnt,
               int root = 0);
    template<typename T>
    int Allgather(const T* sendbuf, int sendcnt, T* rcvbuf, int rcvcnt);
    template<typename T>
    int Scatter(const T* sendbuf, int sendcnt, T* rcvbuf, int rcvcnt,
                int root = 0);
    template<typename T>
    // Don't pass a request if it isn't needed for test/wait.
    int Isend(const T* buf, int count, int dest, int tag,
              MPI_Request* req = NULL);
    template<typename T>
    int Send(const T* buf, int count, int dest, int tag = 0);
    template<typename T>
    int Recv(T* buf, int count, int src, int tag = 0);
    template<typename T>
    int Irecv(T* buf, int count, int src, int tag, MPI_Request* req = NULL);
    int Waitall(int count, MPI_Request* reqs, MPI_Status* stats);
    int Waitany(int count, MPI_Request* reqs, int* index, MPI_Status* stats);
    int Test(MPI_Request* req, int* flag, MPI_Status* stat);
    int Request_free(MPI_Request* req);

    // -------------------------------------------------------------------------
    // The client implements the ParforManager and ParforWorker interfaces.
    class ParforManager {
    public:
      virtual ~ParforManager() {}
      // *Non-blocking* send whatever is needed to worker dest. If the worker
      // just needs to know the job index job_idx, then do nothing and return 0;
      // otherwise, return the output of the MPI routine on failure.
      //   dest and src have ranks 1:GetNProc(). *Blocking* send here or receive
      // in the next method breaks Parfor.
      virtual int Isend(int job_idx, int dest) = 0;
      // *Non-blocking* receive data from worker src for job index job_idx. If
      // the manger doesn't need to receive anything, then return 0.
      virtual int Irecv(int job_idx, int src) = 0;
      // Optional routine that is called when job_idx is completed.
      virtual void IsDone(int worker, int job_idx) {}

    private:
      static const int _tag = 101;
    };

    class ParforWorker {
    public:
      virtual ~ParforWorker() {}
      // Fulfill the job with index job_idx. Receive any data sent in
      // ParformManager::Isend from rank root. Send anything you want for for
      // receipt in ParformManager::Irecv to root. Return the result of the MPI
      // send, or 0.
      virtual int Work(int job_idx, int root) = 0;
    };

    // Then the client calls this function. pw is NULL on the root node. pm is
    // NULL on all the others. njobs need not be valid if pm is NULL.
    int Parfor(ParforManager* pm, ParforWorker* pw, int njobs);

    // -------------------------------------------------------------------------
    // Segment an array into equal pieces. Even if a computational domain has
    // greater than one dimension, one organizes the array so that each chunk is
    // contiguous in memory. Hence a 1D model here is all that is needed.
    class ArraySegmenter {
    public:
      ArraySegmenter ();

      // Segment with roughly n per segment.
      void ApportionN(int n);
      // Segment so that mpi::Pid() has n.
      void ApportionToMe(int n);

      int GetNproc () const { return _bds.size() - 1; }

      // The caller is responsible for indices bounds[0]:bounds[1]-1. Hence a
      // typical for loop appears as
      //     as.GetIndexBounds(bds, 2);
      //     for (int i = bds[0]; i < bds[1]; i++)
      void GetIndexBounds(int bounds[2], int factor = 1) const;

      // This is the offset from index 0.
      int GetOffset () const;
      // Size of the block associated with the calling process.
      int GetN () const;
      int GetNtot () const { return _bds.back(); }
      // True if every block is the same size.
      bool AllEqualSize () const { return _all_equal_size; }
      bool AmRoot () const { return _am_root; }

      // Like their MPI_* counterparts, except they compensate for the uneven
      // sizes of the blocks.
      template<typename T>
      int Gather (const T* sndbuf, T* rcvbuf, int root, int factor = 1) const;
      template<typename T>
      int Allgather (const T* sndbuf, T* rcvbuf, int factor = 1) const;
      template<typename T>
      int Scatter (const T* sndbuf, T* rcvbuf, int root, int factor = 1) const;

    private:
      vector<int> _bds;
      bool _all_equal_size;
      bool _am_root;
      // Keep some buffers at the ready to minimize allocs.
      mutable vector<MPI_Request> _reqs;
      mutable vector<MPI_Status> _stats;
    };

    // -------------------------------------------------------------------------
    // Write arbitary data to a byte buffer. This is useful to serialize an
    // object into a single MPI message. HOWEVER, my treatment here is terrible
    // because I reduce everything to a stream of bytes, which can break if
    // there are computers having different endianness on the network. I should
    // migrate to Boost.MPI.

    class ByteBufferWriter {
    public:
      ByteBufferWriter(size_t init_sz = 1);

      // Reset the buffer. The buffer's capacity remains as large as it was
      // before the call.
      void Reset();

      // Write data to the buffer. Advance buffer pointer by n*sizeof(T).
      template<typename T> void Write(const T* data, size_t n);
      template<typename T> void WriteScalar(T data);
      
      // Get the current size of the buffer (the part containing data).
      size_t Size() const;

      // Get the address of the beginning of the buffer.
      const char* GetPtr() const;

      // A bit awkward, but useful. Copy to this location and advance the buffer
      // by n.
      template<typename T> T* CopyTo(size_t n);
      // But then allow erasure of part. These two are helpful when reading from
      // a file using fread and then sending.
      template<typename T> void Resize(size_t n);

      void Bcast();
      void Isend(int dest);
      void Send(int dest);

    private:
      vector<char> _buf;
    };

    // Read typed data from the buffer. This is useful to unpack data from an
    // MPI message to build an object.
    class ByteBufferReader {
    public:
      ByteBufferReader(const char* attach_buf);

      // Copy the next n*sizeof(T) bytes into data. data must be
      // preallocated. No bounds checking.
      template<typename T> void Read(T* data, size_t n);

      const char* CurrentPosition() const;

      // View a scalar in the buffer, interpreted as type T, from any position
      // pos, where pos refers to the position in bytes. No bounds
      // checking. Alignment is not an issue, though, because memory is copied
      // to a temp.
      template<typename T> T ViewScalar(size_t pos) const;

      // Mix of View and Read: view the data, and advance the buffer pointer by
      // sizeof(T).
      template<typename T> T ViewScalarAndAdvance();

    protected:
      const char* _buf;
      size_t _i;
    };

    class ByteBufferReaderMpi : public ByteBufferReader {
    public:
      ByteBufferReaderMpi();

      // Every call to these two routines resets the buffer and buffer pointer.
      void Bcast();
      void Recv(int src);

      // Size in bytes.
      size_t Size() const;

    protected:
      vector<char> _data;
    };

  }
}

#include "Mpi_inl.hpp"

#endif
