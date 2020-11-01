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

#include <ctype.h>
#include "mex.h"
#include "hmmvp/include/Hmat.hpp"
using namespace std;

typedef util::UInt UInt;

// Sits in memory between mex calls.
struct Data {
  hmmvp::Hmat* hm;
  UInt nthreads;
  bool reorganize_memory;

  Data () : hm(0), nthreads(0), reorganize_memory(false) {}
  void Clear () {
    if (hm) DeleteHmat(hm);
    hm = 0;
    nthreads = 0;
    reorganize_memory = false;
  }
  void NewHmat (const string& filename, UInt ncol, UInt max_nthreads) {
    Clear();
    hm = hmmvp::NewHmat(filename, ncol, max_nthreads);
    nthreads = max_nthreads;
  }
  UInt SetThreads (size_t new_nthreads) {
    if (!hm) return 0;
    const UInt old_nthreads = nthreads;
    nthreads = hm->SetThreads(new_nthreads);
    if (reorganize_memory && nthreads != old_nthreads)
      hm->ReorganizeMemory();
    return nthreads;
  }
  bool ReorganizeMemory () {
    if (hm) {
      hm->ReorganizeMemory();
      reorganize_memory = true;
    }
  }
};

vector<Data> ghms;

void Cleanup (UInt id) {
  if (id + 1 > ghms.size()) mexErrMsgTxt("Cleanup: Invalid id.");
  ghms[id].Clear();
  if (id == ghms.size() - 1) ghms.pop_back();
}

mxArray* Mvp (hmmvp::Hmat& h, const double* x, UInt m, UInt ncol,
              const vector<hmmvp::Blint>* rs = 0,
              const vector<hmmvp::Blint>* cs = 0) {
  if (h.GetN() != m || ncol == 0)
    mexErrMsgTxt("mvp: x is the wrong size.");

  mxArray* my = mxCreateDoubleMatrix(h.GetM(), ncol, mxREAL);
  try {
    h.Mvp(x, mxGetPr(my), ncol, rs, cs);
  } catch (hmmvp::Exception& e) {
    mexErrMsgTxt(e.GetMsg().c_str());
  }
  return my;
}

mxArray* MvpT (hmmvp::Hmat& h, const double* x, UInt m, UInt ncol,
               const vector<hmmvp::Blint>* rs = 0,
               const vector<hmmvp::Blint>* cs = 0) {
  if (h.GetM() != m || ncol == 0)
    mexErrMsgTxt("tmvp: x is the wrong size.");
  if (rs || cs)
    mexErrMsgTxt("MvpT: rs and cs functionality is not implemented; sorry.");

  mxArray* my = mxCreateDoubleMatrix(h.GetN(), ncol, mxREAL);
  try {
#if 0
    if (rs) {
      if (cs) h.MvpT(x, mxGetPr(my), ncol, *rs, *cs);
      else    h.MvpT(x, mxGetPr(my), ncol, *rs);
    } else
#endif
    h.MvpT(x, mxGetPr(my), ncol);
  } catch (hmmvp::Exception& e) {
    mexErrMsgTxt(e.GetMsg().c_str());
  }
  return my;
}

int Init (const char *filename, UInt nthreads, UInt ncol, UInt& nnz) {
  if (nthreads < 0) mexErrMsgTxt("Init: Invalid nthreads.");
  if (ncol < 1) mexErrMsgTxt("Init: Invalid ncol.");

  // Get a free id
  int id = -1;
  for (UInt i = 0; i < ghms.size(); i++)
    if (!ghms[i].hm) {
      id = i;
      break;
    }
  if (id == -1) {
    ghms.push_back(Data());
    id = ghms.size() - 1;
  }

  const char* exts[] = {"", ".hm", ".hmat"};
  for (size_t i = 0, N = sizeof(exts)/sizeof(exts[0]); i < N; i++) {
    try {
      ghms[id].NewHmat(string(filename) + exts[i], ncol, nthreads);
      break;
    } catch (hmmvp::FileException e) {
      if (i == N - 1) {
        // Every extension failed.
        Cleanup(id);
        mexErrMsgTxt(e.GetMsg().c_str());
      }
    }
  }
  nnz = ghms[id].hm->GetNnz();

  return id;
}

void mexExitFunction () {
  for (UInt i = 0; i < ghms.size(); i++) Cleanup(i);
}

inline void GetIdxs (const mxArray* ma, vector<hmmvp::Blint>& rs, UInt maxr) {
  UInt nrs = mxGetNumberOfElements(ma);
  const double* drs = mxGetPr(ma);
  rs.resize(nrs);
  for (UInt i = 0; i < nrs; i++) {
    rs[i] = (hmmvp::Blint)drs[i] - 1;
    if (rs[i] < 0 || rs[i] > maxr)
      mexErrMsgTxt("Index out of bounds.");
  }
}

inline hmmvp::Hmat& GetHmat (const mxArray* ma, UInt* rid = NULL) {
  if (mxGetNumberOfElements(ma) != 1)
    mexErrMsgTxt("numel(id) should be 1.");
  UInt id = (UInt) mxGetScalar(ma);
  if (ghms.empty() || id > ghms.size() - 1 || !ghms[id].hm)
    mexErrMsgTxt("mvp: Invalid id.");
  if (rid) *rid = id;
  return *ghms[id].hm;
}

// Implements
//   [id nnz] = hmmvp('init',filename,[nthreads],[ncol])
//   hmmvp('cleanup',id)
//   y = hmmvp('mvp',id,x,[rs],[cs])
//     rs,cs are optional subsets of rows, cols. indexing starts at 1. Just rs
//     can be specified. See Hmat.hpp for more details.
//   y = hmmvp('tmvp',id,x,[rs],[cs])
//     Transpose MVP.
//   nthreads = hmmvp('threads',id,nthreads);
//   hmmvp('stsave',id);
//   hmmvp('strelease',id);
//   hmmvp('e2bon',id);
//   hmmvp('e2boff',id);
//   hmmvp('pon',id);
//   hmmvp('poff',id);
//   nf2 = hmmvp('fronorm2',id)
//   n1 = hmmvp('norm1',id)
//   Brc = hmmvp('extract',id,rs,cs).
//   [I J S] = hmmvp('bblocks',id,[cutoff]);
//   m/n = hmmvp('getm/n',id);
//   hmmvp('reorganizememory',id);
void mexFunction (int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  if (nrhs == 0) mexErrMsgTxt("Type 'help hmmvp' for help.");

  // Clean up command.
  int strlen = mxGetNumberOfElements(prhs[0]) + 1;
  vector<char> vfn(strlen);
  mxGetString(prhs[0], &vfn[0], strlen);
  for (int i = 0; i < strlen - 1; i++) vfn[i] = tolower(vfn[i]);
  string fn(&vfn[0]);

  // init
  if (fn[0] == 'i' || fn == "init") {
    mexAtExit(&mexExitFunction);
    if (nrhs < 2 || nrhs > 4 || nlhs < 1 || !mxIsChar(prhs[1]))
      mexErrMsgTxt("[id nnz] = hmmvp('init',filename,[nthreads],[ncol])");
    UInt n = mxGetNumberOfElements(prhs[1]) + 1;
    char* filename = new char[n];
    mxGetString(prhs[1], filename, n);
    UInt nthreads = 1;
    if (nrhs > 2) {
      if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
	mexErrMsgTxt("nthreads should be a positive integer");
      nthreads = (int)mxGetScalar(prhs[2]);
    }
    UInt ncol = 1;
    if (nrhs > 3) {
      if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
	mexErrMsgTxt("ncol should be a positive integer");
      ncol = (int)mxGetScalar(prhs[3]);
    }
    UInt nnz;
    UInt id = Init(filename, nthreads, ncol, nnz);
    delete[] filename;
    plhs[0] = mxCreateDoubleScalar((double)id);
    if (nlhs > 1) plhs[1] = mxCreateDoubleScalar((double)nnz);
  }

  // cleanup
  else if (fn[0] == 'c' || fn == "cleanup") {
    if (nrhs != 2 || !mxIsDouble(prhs[1])) mexErrMsgTxt("hmmvp('cleanup',id)");
    UInt id = (int)mxGetScalar(prhs[1]);
    Cleanup(id);
  }

  // mvp
  else if (fn[0] == 'm' || fn == "mvp") {
    if (nrhs < 3 || nrhs > 5 || nlhs != 1 || !mxIsDouble(prhs[1]) ||
	!mxIsDouble(prhs[2]) || (nrhs >= 4 && !mxIsDouble(prhs[3])) ||
	(nrhs == 5 && !mxIsDouble(prhs[4])))
      mexErrMsgTxt("y = hmmvp('mvp',id,x,[rs],[cs])");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    const mxArray* mx = prhs[2];
    vector<hmmvp::Blint> *prs = NULL, *pcs = NULL;
    vector<hmmvp::Blint> rs, cs;
    if (nrhs >= 4) {
      GetIdxs(prhs[3], rs, h.GetM() - 1);
      if (!rs.empty()) prs = &rs;
      if (nrhs == 5) {
        GetIdxs(prhs[4], cs, h.GetN() - 1);
        if (!cs.empty()) pcs = &cs;
      }
    }
    plhs[0] = Mvp(h, mxGetPr(mx), mxGetM(mx), mxGetN(mx), prs, pcs);
  }

  // tranpose mvp (not supported yet)
#if 0
  else if (fn == "tmvp") {
    if (nrhs < 3 || nrhs > 5 || nlhs != 1 || !mxIsDouble(prhs[1]) ||
	!mxIsDouble(prhs[2]) || (nrhs >= 4 && !mxIsDouble(prhs[3])) ||
	(nrhs == 5 && !mxIsDouble(prhs[4])))
      mexErrMsgTxt("y = hmmvp('tmvp',id,x,[rs],[cs])");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    const mxArray* mx = prhs[2];
    vector<hmmvp::Blint> *prs = NULL, *pcs = NULL;
    vector<hmmvp::Blint> rs, cs;
    if (nrhs >= 4) {
      GetIdxs(prhs[3], rs, h.GetM() - 1);
      if (!rs.empty()) prs = &rs;
      if (nrhs == 5) {
        GetIdxs(prhs[4], cs, h.GetN() - 1);
        if (!cs.empty()) pcs = &cs;
      }
    }
    plhs[0] = MvpT(h, mxGetPr(mx), mxGetM(mx), mxGetN(mx), prs, pcs);
  }
#endif

  // Frobenius norm
  else if (fn == "fronorm2") {
    if (nlhs != 1 || nrhs > 2 || !mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("nf2 = hmmvp('fronorm2',id)");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    plhs[0] = mxCreateDoubleScalar(h.NormFrobenius2());
  }

  // 1-norm
  else if (fn == "norm1") {
    if (nlhs != 1 || nrhs > 2 || !mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("nf2 = hmmvp('norm1',id)");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    plhs[0] = mxCreateDoubleScalar(h.NormOne());
  }

  // State
  else if (fn.substr(0, 2) == "st") {
    if (nrhs < 2 || mxGetNumberOfElements(prhs[0]) < 3)
      mexErrMsgTxt("mvp: Invalid function");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    switch (fn[2]) {
    case 's':
      h.SaveState();
      break;
    case 'r':
      h.ReleaseState();
      break;
    }
  }

  // Permutations
  else if (fn.substr(0, 2) == "po") {
    if (nrhs < 2 || mxGetNumberOfElements(prhs[0]) < 3)
      mexErrMsgTxt("mvp: Invalid function");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    switch (fn[2]) {
    case 'n': case 'N':
      h.TurnOnPermute();
      break;
    case 'f': case 'F':
      h.TurnOffPermute();
      break;
    }
  }

  // Element-to-block map
  else if (fn.substr(0, 4) == "e2bo") {
    if (nrhs < 2 || mxGetNumberOfElements(prhs[0]) < 5)
      mexErrMsgTxt("mvp: Invalid function");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    switch (fn[4]) {
    case 'n': case 'N':
      h.TurnOnElementToBlockMap();
      break;
    case 'f': case 'F':
      h.TurnOffElementToBlockMap();
      break;
    }    
  }

  // Extract
  else if (fn == "extract") {
    if (nlhs != 1 || nrhs != 4 || !mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("Brc = hmmvp('extract',id,rs,cs)");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    vector<hmmvp::Blint> rs, cs;
    GetIdxs(prhs[2], rs, h.GetM() - 1);
    GetIdxs(prhs[3], cs, h.GetN() - 1);
    plhs[0] = mxCreateDoubleMatrix(rs.size(), cs.size(), mxREAL);
    h.Extract(rs, cs, mxGetPr(plhs[0]));
  }

  // [I J S] triple of full-block blocks
  else if (fn == "bblocks") {
    if (nlhs != 3 || nrhs < 2 || !mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("[I J S] = hmmvp('bblocks',id,[cutoff])");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    double cutoff = -1;
    if (nrhs > 2) cutoff = (double) mxGetScalar(prhs[2]);
    vector<hmmvp::Blint> I, J;
    vector<double> S;
    h.FullBlocksIJS(I, J, S, cutoff);
    size_t n = I.size();
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double* p = mxGetPr(plhs[0]);
    for (size_t i = 0; i < n; i++) p[i] = (double) I[i];
    I.clear();
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    p = mxGetPr(plhs[1]);
    for (size_t i = 0; i < n; i++) p[i] = (double) J[i];
    J.clear();
    plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL);
    p = mxGetPr(plhs[2]);
    for (size_t i = 0; i < n; i++) p[i] = (double) S[i];
    S.clear();
  }

  // Number of threads
  else if (fn == "threads") {
    if (nlhs != 1 || nrhs != 3 || mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("nthreads = hmmvp('threads',id,nthreads)");
    UInt id;
    hmmvp::Hmat& h = GetHmat(prhs[1], &id);
    int nthreads = (int) mxGetScalar(prhs[2]);
    if (nthreads < 1) mexErrMsgTxt("nthreads must be > 0");
    nthreads = ghms[id].SetThreads(nthreads);
    plhs[0] = mxCreateDoubleScalar((double) nthreads);
  }

  else if (fn == "getm" || fn == "getn") {
    if (nlhs != 1 || nrhs > 2 || !mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("m/n = hmmvp('getm/n',id)");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    if (fn[3] == 'm')
      plhs[0] = mxCreateDoubleScalar(h.GetM());
    else
      plhs[0] = mxCreateDoubleScalar(h.GetN());
  }

  else if (fn == "nnz") {
    if (nlhs != 1 || nrhs != 4 || !mxIsDouble(prhs[1]) ||
	mxGetNumberOfElements(prhs[1]) != 1)
      mexErrMsgTxt("nnz = hmmvp('nnz',id,rs,cs)");
    hmmvp::Hmat& h = GetHmat(prhs[1]);
    vector<hmmvp::Blint> rs, cs;
    GetIdxs(prhs[2], rs, h.GetM() - 1);
    GetIdxs(prhs[3], cs, h.GetN() - 1);
    hmmvp::Blint nnz = h.GetNnz(rs, cs);
    plhs[0] = mxCreateDoubleScalar((double) nnz);
  }

  else if (fn.substr(0, 5) == "reorg") {
    if (nrhs < 2) mexErrMsgTxt("hmmvp('reorganizememory',id)");
    UInt id;
    GetHmat(prhs[1], &id);
    ghms[id].ReorganizeMemory();
  }
  
  else {
    mexErrMsgTxt((string("Invalid function: ") + fn).c_str());
  }
}
