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

/* hmmvpbuild is an example driver to create an H-matrix. It can be run as a
   serial code, as a parallel code using OpenMP, or as a parallel code using
   MPI.

   Run ./hmmvpbuild on the command line without arguments to see documentation.

   In your particular application, you will likely want to use this file as an
   example of how to use hmmvp, but will customize many parts to your
   setting. In particular, your application likely handles problem setup
   differently than is done here.

   The objective of this driver is to demonstrate how to integrate a new Green's
   function. Look for the string DEV in this file to see how to add a new
   Green's function.

   Build as follows. For parallelization by MPI:
       make build mode=p opt=-O
   For serial or OpenMP with no MPI dependencies:
       make build mode=s opt=-O
   Edit the top lines of Makefile to profile pointers to LAPACK and BLAS and
   compilers.
*/

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "util/include/Mpi.hpp"
#include "util/include/KeyValueFile.hpp"
#include "util/include/Util.hpp"
#include "util/include/CodeAnalysis.hpp"
#include "hmmvp/include/Hd.hpp"
#include "hmmvp/include/Compress.hpp"
#include "Help.hpp"
using namespace std;
using namespace hmmvp;
using namespace util;

/* DEV Implement Green's functions here. To add a new Green's function, follow
   these steps:
     
   1. Inherit from ImplGreensFn or one of its descendants.
   2. Add the new string identifier to NewGreensFn().

   You needn't worry about MPI or OpenMP parallelization (with one caveat for
   the latter). Implement your class as though this is a serial code. For
   examples, see the classes that follow.
     If you are using OpenMP parallelization, your Green's function needs to be
   thread-safe. At the interface level, 'Call' is const as a first step to
   thread safety, but at the very least, mutable variables can break thread
   safety. More subtly, Fortran routines often are not thread safe, most
   obviously in the form of COMMON blocks. Use !$OMP THREADPRIVATE annotations
   to protect COMMON blocks. */

class ImplGreensFn : public GreensFn {
public:
  virtual ~ImplGreensFn () {}

  /* DEV Use the key-value file to initialize the Green's function. If
     something in the input is amiss, throw an exception stating the
     problem. */
  virtual void Init(const KeyValueFile* kvf) throw (Exception) = 0;

  /* DEV Return an Hd corresponding to the model discretization. First
     determine the element centers. Then return the result of
     hmmvp::NewHd(). Only the root process calls this method. */
  virtual Hd* ComputeHd(double eta) = 0;

  /* DEV Implement these if you want to compute Green's functions differently
     when estimating ||B||_F than when forming the H-matrix. For example, a
     lower-quality approximation to a source-reciever GF may be possible when
     estimating ||B||_F. */
  virtual void AmEstimatingBfro () {}
  virtual void AmFormingHmat () {}

  /* DEV Compute B(rs,cs). Indexing starts at 1. B is preallocated.  Return true
     if all is well; false if there is an error in computing the Green's
     function and you want compression to stop. You almost certainly want to
     ignore 'cbi'. */
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const = 0;

  /* DEV Optionally do extra tasks, such as computing the boundary
     conditions. The 'Serial' version is called only by the root process. You
     can include OpenMP/pthread statements in the 'Serial' version. The 'Mpi'
     version is called by all processes. */
  virtual void DoExtraTasksSerial () throw (Exception) {}
  virtual void DoExtraTasksMpi () throw (Exception) {}
};

/* DEV For neatness, you might put your code into a file whose name
   corresponds to the string identifier and include it here. */
#include "GreensFnInverseR.cpp"
#include "GreensFnOkada92.cpp"
#include "GreensFnTgf.cpp"

ImplGreensFn* NewGreensFn (const string& id, const KeyValueFile* kvf)
  throw (Exception)
{
  ImplGreensFn* gf = NULL;
  /* DEV Add your GF's string identifier here. */
  if (id == "inverse-r") {
    gf = new InverseRGreensFn();
  } else if (id == "okada92") {
    gf = new gfb::GreensFnOkada92();
  } else if (id == "gimbutas+12") {
    gf = new tgf::GreensFn();
  } else {
    throw Exception("No such Green's function string identifier.");
  }
  if (!gf) throw Exception("gf is NULL.");
  try {
    gf->Init(kvf);
    return gf;
  } catch (const Exception& e) {
    delete gf;
    throw e;
  }
}

/* DEV You should not have to modify anything below here to run hmmvpbuild.
   However, you might use the following code as an example of how to write your
   own application-dependent driver.
     I use a util::KeyValueFile to communicate problem setup data to
   hmmvpbuild. A KeyValueFile is independent of libhmmvp and so just one
   optional problem setup approach. An established simulator quite likely
   already has its own mechanism for representing problem data at the level of
   communication between modules. */

namespace {
int Finalize (int code) {
  mpi::Finalize();
  return code;
}

int Error (const string& s, KeyValueFile* kvf) {
  fprintf(stderr, "No key %s.\n", s.c_str());
  return Finalize(-1);
}

inline string RemoveExtensionIfAny (const string& fn, const string& ext) {
  size_t p;
  string fno;
  if ((p = fn.rfind("." + ext)) != string::npos)
    fno = fn.substr(0, p);
  else fno = fn;
  return fno;
}

struct Inputs {
  bool allow_overwrite, do_extra_tasks_only, write_Bfro_only;
  string base_fn, write_hmat_filename, use_hmat_filename;
  string write_hd_filename, use_hd_filename;
  string greens_fn_name;
  string err_method;
  double eta, tol, Bfro;
  size_t nthreads;
};

int ProcessKvf (KeyValueFile* kvf, Inputs& in) {
  double d;
  const string* s;
  bool am_root = mpi::AmRoot();

  in.allow_overwrite = false;
  if (kvf->GetDouble("allow_overwrite", d))
    in.allow_overwrite = (bool) d;

  in.do_extra_tasks_only = false;
  if (kvf->GetDouble("do_extra_tasks_only", d))
    in.do_extra_tasks_only = (bool) d;

  in.write_Bfro_only = false;
  if (kvf->GetDouble("write_Bfro_only", d))
    in.write_Bfro_only = (bool) d;

  in.write_hmat_filename = "";
  if (!kvf->GetString("write_hmat_filename", s))
    return Error("write_hmat_filename", kvf);
  in.base_fn = RemoveExtensionIfAny(*s, "hm");
  in.write_hmat_filename = in.base_fn + ".hm";
  if (!in.allow_overwrite &&
      access(in.write_hmat_filename.c_str(), F_OK) == 0) {
    if (am_root)
      fprintf(stderr, "%s exists; won't overwrite.\n",
              in.write_hmat_filename.c_str());
    return Finalize(-1);
  }

  in.use_hd_filename = "";
  in.write_hd_filename = "";
  if (kvf->GetString("use_hd_filename", s))
    in.use_hd_filename = RemoveExtensionIfAny(*s, "hd") + ".hd";
  if (kvf->GetString("write_hd_filename", s))
    in.write_hd_filename = RemoveExtensionIfAny(*s, "hd") + ".hd";
  if (!in.use_hd_filename.empty() && !in.write_hd_filename.empty()) {
    if (am_root)
      fprintf(stderr, "Must provide none or only one of use_hd_filename or "
              "write_hd_filename.\n");
    return Finalize(-1);
  }
  if (!in.allow_overwrite &&
      access(in.write_hd_filename.c_str(), F_OK) == 0) {
    if (am_root)
      fprintf(stderr, "%s exists; won't overwrite.\n",
              in.write_hd_filename.c_str());
    return Finalize(-1);
  }

  if (!kvf->GetString("err_method", s)) in.err_method = "mrem-fro";
  else {
    in.err_method = *s;
    if (!(in.err_method == "mrem-fro" || in.err_method == "brem-fro" ||
          in.err_method == "mrem-abs")) {
      if (am_root)
        fprintf(stderr, "err_method must be one of: mrem-fro, brem-fro, "
                "mrem-abs.\n");
      return Finalize(-1);
    }
  }

  in.tol = 1e-6;
  kvf->GetDouble("tol", in.tol);
  if (in.tol <= 0) {
    if (am_root) fprintf(stderr, "tol must be > 0.\n");
    return Finalize(-1);
  }

  in.eta = 3;
  kvf->GetDouble("eta", in.eta);
  if (in.eta <= 0) {
    if (am_root) fprintf(stderr, "eta must be > 0.\n");
    return Finalize(-1);
  }

  if (!kvf->GetString("greens_fn", s)) return Error("greens_fn", kvf);
  in.greens_fn_name = *s;

  in.Bfro = -1.0;
  if (kvf->GetDouble("Bfro", in.Bfro) && in.Bfro <= 0.0) {
    if (am_root) fprintf(stderr, "Bfro must be > 0.\n");
    return Finalize(-1);
  }

  in.use_hmat_filename = "";
  if (kvf->GetString("use_hmat_filename", s))
    in.use_hmat_filename = RemoveExtensionIfAny(*s, "hm") + ".hm";
  if (in.use_hmat_filename == in.write_hmat_filename) {
    if (am_root) fprintf(stderr, "use_hmat_filename == write_hmat_filename.\n");
    return Finalize(-1);
  }

  in.nthreads = 1;
  if (mpi::GetNproc() == 1 && kvf->GetDouble("nthreads", d) && d >= 1)
    in.nthreads = (size_t) d;

  return 0;
}

double SetBfro (Compressor* c, double in_Bfro) {
  double b;
  if (in_Bfro > 0)
    c->SetBfroEstimate(b = in_Bfro);
  else if (c->HaveOldHmat())
    c->SetBfroEstimate(b = c->GetOldHmatBfro());
  else
    c->SetBfroEstimate(b = c->EstimateBfro());
  return b;
}

void
WriteBfroToFile (const string& fn, const double Bfro) throw (FileException) {
  FILE* fid = fopen(fn.c_str(), "wa");
  if (!fid) throw FileException(fn + string(" can't be written."));
  fprintf(fid, "%1.15e\n", Bfro);
  fclose(fid);
}

}

int main (int argc, char** argv) {
  mpi::Init(argc, argv);
  bool am_root = mpi::AmRoot();

  if (argc == 1 || (argc == 2 && string(argv[1]) == "help")) {
    if (am_root) {
      printf(_hmmvpbuild_PrintHelp_text_);
      printf(_hmmvp_Header_text_);
    }
    return Finalize(-1);
  } else if (argc == 3 && string(argv[1]) == "help") {
    if (am_root) {
      string cmd = argv[2];
      if (cmd == "compress") {
        printf(_hmmvpbuild_compress_PrintHelp_text_);
        printf(_hmmvp_Header_text_);
      } else {
        printf("%s is not a valid command.\n\n", cmd.c_str());
        printf(_hmmvpbuild_PrintHelp_text_);
        printf(_hmmvp_Header_text_);
      }
    }
    return Finalize(-1);
  }

  // Handle key-value file.
  Inputs in;
  KeyValueFile* kvf;
  { kvf = NewKeyValueFile();
    const string kvf_fn = argv[1];
    if (!kvf->Read(kvf_fn) && !kvf->Read(kvf_fn + ".kvf")) {
      if (am_root)
        fprintf(stderr, "Failed to read the key-value file %s.\n",
                kvf_fn.c_str());
      DeleteKeyValueFile(kvf);
      return Finalize(-1);
    }
    int ret;
    if ((ret = ProcessKvf(kvf, in))) {
      DeleteKeyValueFile(kvf);
      return ret;
    }
  }

  // Create the Green's function.
  if (am_root) printf("... Initializing the Green's function.\n");
  ImplGreensFn* gf = NULL;
  try {
    gf = NewGreensFn(in.greens_fn_name, kvf);
  } catch (const Exception& e) {
    if (am_root)
      fprintf(stderr, "From greens_fn_name = %s comes the exception: %s\n",
              in.greens_fn_name.c_str(), e.GetMsg().c_str());
  }
  DeleteKeyValueFile(kvf);
  if (!gf) return Finalize(-1);
  
  if (am_root) printf("... Doing extra tasks.\n");
  if (am_root) {
    try {
      gf->DoExtraTasksSerial();
      mpi_IsTrue(true);
    } catch (const Exception& e) {
      mpi_IsTrue(false);
      fprintf(stderr, "DoExtraTasksSerial gives the exception: %s\n",
              e.GetMsg().c_str());
      delete gf;
      return Finalize(-1);
    }
  } else if (!mpi_IsTrue()) {
    delete gf;
    return Finalize(-1);
  }
  try {
    gf->DoExtraTasksMpi();
    mpi::Barrier();
  } catch (const Exception& e) {
    if (am_root)
      fprintf(stderr, "DoExtraTasksMpi gives the exception: %s\n",
              e.GetMsg().c_str());
    delete gf;
    return Finalize(-1);
  }
  if (in.do_extra_tasks_only) {
    if (am_root) printf("... Exiting after doing extra tasks only.\n");
    delete gf;
    return Finalize(0);
  }
  
  // Spatial decomposition.
  hmmvp::Hd* hd = NULL;
  if (am_root) {
    if (!in.use_hd_filename.empty()) {
      try {
        hd = NewHd(in.use_hd_filename);
      } catch (const FileException& e) {
        fprintf(stderr, "Could not read %s.\n", in.use_hd_filename.c_str());
        delete gf;
        return Finalize(-1);
      } catch (const Exception& e) {
        fprintf(stderr, "NewHd, given %s, gives the exception: %s\n",
                in.use_hd_filename.c_str(), e.GetMsg().c_str());
        delete gf;
        return Finalize(-1);
      }
    } else {
      printf("... Computing spatial decomposition.\n");
      hd = gf->ComputeHd(in.eta);
      if (!in.write_hd_filename.empty())
        try {
          WriteHd(hd, in.write_hd_filename);
        } catch (const FileException& e) {
          fprintf(stderr, "Could not write %s.\n",
                  in.write_hd_filename.c_str());
          delete gf;
          return Finalize(-1);
        }
    }
  }

  // Set up Compressor.
  Compressor* c = NULL;
  if (am_root) printf("... Setting up compressor.\n");
  try {
    c = hmmvp::NewCompressor(hd, gf);
    c->SetOutputLevel(1);
    if (in.err_method == "mrem-fro")
      c->SetTolMethod(Compressor::tm_mrem_fro);
    else if (in.err_method == "brem-fro")
      c->SetTolMethod(Compressor::tm_brem_fro);
    else if (in.err_method == "mrem-abs")
      c->SetTolMethod(Compressor::tm_mrem_abs);
    else {
      // Protected by check in ProcessKvf.
      assert(false);
    }
    c->SetTol(in.tol);
    if (mpi::GetNproc() == 1) c->SetOmpNthreads(in.nthreads);
    if (!in.use_hmat_filename.empty())
      try {
        c->UseHmatFile(in.use_hmat_filename);
      } catch (const Exception& e) {
        if (am_root)
          fprintf(stderr, "Warning: Ignoring old H-matrix file %s; received "
                  "the exception: %s\n", in.use_hmat_filename.c_str(),
                  e.GetMsg().c_str());
      }
    if (c->GetTolMethod() == Compressor::tm_mrem_fro) {
      gf->AmEstimatingBfro();
      in.Bfro = SetBfro(c, in.Bfro);
      gf->AmFormingHmat();
    }
  } catch (const Exception& fe) {
    if (am_root) {
      fprintf(stderr, "Exception creating Compressor: %s\n",
              fe.GetMsg().c_str());
      DeleteHd(hd);
    }
    return Finalize(-1);
  }

  if (am_root) DeleteHd(hd);

  if (in.write_Bfro_only) {
    if (am_root) {
      WriteBfroToFile(in.base_fn + ".Bfroest", in.Bfro);
      printf("... Exiting without compressing at user request "
             "(write_Bfro_only 1).\n");
    }
  } else {
    // Compress.
    if (am_root) printf("... Compressing.\n");
    try {
      c->CompressToFile(in.write_hmat_filename);
    } catch (const UserReqException& e) {
      if (am_root) fprintf(stderr, "Green's function requested stop.\n");
    } catch (const Exception& e) {
      fprintf(stderr, "Compressor threw: %s\n", e.GetMsg().c_str());
    }
  }

  DeleteCompressor(c);
  delete gf;
  return Finalize(0);
}
