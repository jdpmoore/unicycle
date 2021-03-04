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

#include <stdio.h>
#include "util/include/Mpi.hpp"

namespace util {
namespace mpi {

#ifdef UTIL_MPI
static void MPI_Handler_Crash (MPI_Comm*, int*, ...) {
  // Make valgrind dump out information.
  fprintf(stderr, "Forcing a seg fault. Buckle up.\n");
  double* a = NULL;
  *a = 3.14;
}
#endif

void Init (int argc, char** argv) {
#ifdef UTIL_MPI
  MPI_Init(&argc, &argv);
  MPI_Errhandler eh;
  MPI_Errhandler_create(MPI_Handler_Crash, &eh);
  MPI_Errhandler_set(MPI_COMM_WORLD, eh);
  MPI_Errhandler_free(&eh);
#else
#endif
}

void Finalize () {
#ifdef UTIL_MPI
  MPI_Finalize();
#else
#endif
}

int GetNproc () {
#ifdef UTIL_MPI
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  return nproc;
#else
  return 1;
#endif
}

int Pid () {
#ifdef UTIL_MPI
  static int pid = -1;
  if (pid < 0) MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  return pid;
#else
  return 0;
#endif
}

void Barrier () {
#ifdef UTIL_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#else
#endif
}

int Root () {
  return 0;
}

bool AmRoot () {
#ifdef UTIL_MPI
  return Pid() == Root();
#else
  return true;
#endif
}

bool AmTPid0() { return Pid() == Root() && omp_get_thread_num() == 0; }

bool IsTrue (bool p) {
#ifdef UTIL_MPI
  char msg;
  if (AmRoot()) msg = (char) p;
  Bcast(&msg, 1);
  return (bool) msg;
#else
  return p;
#endif
}

bool AllOk (bool iok) {
#ifdef UTIL_MPI
  int ok = iok, msg;
  Allreduce<int>(&ok, &msg, 1, MPI_LAND);
  return (bool) msg;
#else
  return iok;
#endif
}

int Waitall (int count, MPI_Request* reqs, MPI_Status* stats) {
#ifdef UTIL_MPI
  return MPI_Waitall(count, reqs, stats);
#else
  return 0;
#endif
}

int Waitany (int count, MPI_Request* reqs, int* index, MPI_Status* stats) {
#ifdef UTIL_MPI
  return MPI_Waitany(count, reqs, index, stats);
#else
  *index = 0;
  return 0;
#endif
}

int Test (MPI_Request* req, int* flag, MPI_Status* stat) {
#ifdef UTIL_MPI
  return MPI_Test(req, flag, stat);
#else
  *flag = 1;
  return 0;
#endif      
}

int Request_free (MPI_Request* req) {
#ifdef UTIL_MPI
  return MPI_Request_free(req);
#else
  return 0;
#endif      
}

int Parfor (ParforManager* pm, ParforWorker* pw, int njobs) {
  const int tag = 100;
  if (mpi::GetNproc() < 2) return -1;

  if (mpi::AmRoot()) {
    // I am the manager.
    int nwrkr = std::min(mpi::GetNproc() - 1, njobs);
    vector<MPI_Status> stats(nwrkr);
    vector<MPI_Request> reqs(nwrkr);
    vector<int> acks(nwrkr, 0), ibuf(nwrkr);
    // Initial distribution of jobs. This primes the various bookkeeping
    // data.
    int ijob = 0, n = nwrkr;
    for (ijob = 0; ijob < n; ijob++) {
      ibuf[ijob] = ijob;
      mpi::Isend(&ibuf[ijob], 1, ijob + 1, tag);
      pm->Isend(ijob, ijob + 1);
      pm->Irecv(ijob, ijob + 1);
      mpi::Irecv(&acks[ijob], 1, ijob + 1, tag, &reqs[ijob]);
    }
    // Now the rest of the jobs. In most cases (njobs >> nwrkr), most of the
    // time is spent in this loop.
    for (;;) {
      // Even though at the end there is likely < n valid reqs, I think
      // Waitany is fine.
      int idx;
      mpi::Waitany(nwrkr, &reqs[0], &idx, &stats[0]);
      // At least one worker has reported back. Poll them all. (If we just
      // deal with worker idx and the workers are faster than the comm,
      // we'll idle all but worker idx.)
      for (int pid = 1; pid <= nwrkr; pid++) {
        // ack == -1 indicates we're no longer using this req because we're
        // wrapping up the parfor.
        if (acks[pid - 1] == -1) continue;
        // Test for receipt of the acknowledgement.
        int flag;
        mpi::Test(&reqs[pid - 1], &flag, &stats[pid - 1]);
        if (flag) {
          // Test takes care of this free:
          //mpi::Request_free(&reqs[pid - 1]);
          n--; // One fewer outstanding job.
          pm->IsDone(pid, ibuf[pid - 1]); // Tell the mgr this job is done.
          if (ijob < njobs) {
            // Tell worker pid what job number it has.
            ibuf[pid - 1] = ijob;
            mpi::Isend(&ibuf[pid - 1], 1, pid, tag);
            // Do application-dependent comm.
            pm->Isend(ijob, pid);
            pm->Irecv(ijob, pid);
            // We'll wait on the returned acknowledgement with value >= 0.
            mpi::Irecv(&acks[pid - 1], 1, pid, tag, &reqs[pid - 1]);
            ijob++; // This will complete another job, but until then ...
            n++;    // ... we have another outstanding job.
          } else if (n == 0) {
            // We've done all the jobs. Break from the pid loop and then
            // break from for (;;).
            break;
          } else {
            // We're nearly done with the jobs. Idle worker pid.
            acks[pid - 1] = -1;
          }
        }
      }
      if (n == 0) break;
    }
    // Tell all the workers it's quittin' time.
    int fin = -1;
    for (int i = 0, np = mpi::GetNproc() - 1; i < np; i++)
      mpi::Send(&fin, 1, i + 1, tag);

  } else {

    // Workers wait for orders and then carry them out. Blocking comm makes
    // sense here, as each step depends on completing the previous.
    for (;;) {
      // Receive a job number.
      int job;
      mpi::Recv(&job, 1, 0, tag);
      // Oh, hey, it's quittin' time!
      if (job < 0) break;
      // Do some work.
      pw->Work(job, 0);
      // Tell the manager I'm done with this job.
      mpi::Send(&job, 1, 0, tag);
    }

  }
  return 0;
}

ArraySegmenter::ArraySegmenter ()
  : _am_root(mpi::AmRoot())
{}

void ArraySegmenter::ApportionN (int n) {
  int nthreads = mpi::GetNproc();

  _bds.resize(nthreads + 1);
  int k = n / nthreads, extra = n - k*nthreads;
  _all_equal_size = extra == 0;
  for (int i = 0, b = 0; i < nthreads; i++) {
    _bds[i] = b;
    b += k;
    if (extra > 0) {
      b++;
      extra--;
    }
  }
  _bds[nthreads] = n;
}

void ArraySegmenter::ApportionToMe (int n) {
  vector<int> ns(mpi::GetNproc());
  mpi::Allgather(&n, 1, &ns[0], ns.size());
  _bds.resize(ns.size() + 1);
  int b = 0;
  for (size_t i = 0; i < ns.size(); i++) {
    _bds[i] = b;
    b += ns[i];
  }
  _bds[ns.size()] = b;

  if (AmRoot()) {
    printf("ns: ");
    for (size_t i = 0; i < ns.size(); i++) printf("%d ", ns[i]);
    printf("\n_bds: ");
    for (size_t i = 0; i < _bds.size(); i++) printf("%d ", _bds[i]);
    printf("\n");
  }
}

void ArraySegmenter::GetIndexBounds (int bounds[2], int factor) const {
  int tid = mpi::Pid();
  bounds[0] = _bds[tid] * factor;
  bounds[1] = _bds[tid + 1] * factor;
}

int ArraySegmenter::GetOffset () const {
  int tid = mpi::Pid();
  return _bds[tid];
}

int ArraySegmenter::GetN () const {
  int tid = mpi::Pid();
  return _bds[tid + 1] - _bds[tid];
}

}}
