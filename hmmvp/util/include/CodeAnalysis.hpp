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

#ifndef INCLUDE_UTIL_CODEANALYSIS
#define INCLUDE_UTIL_CODEANALYSIS

#include <vector>

#ifdef ANALYZE_CODE
#  include <sys/time.h>
#else
// Ugly, but I can't think of a better way to make timeval disappear.
#  define timeval int
#endif

namespace util {
  using namespace std;

#ifdef ANALYZE_CODE
# define caprint(...) printf(__VA_ARGS__)
#else
# define caprint(...)
#endif

#define errpr(...) fprintf(stderr, __VA_ARGS__);

#ifndef NDEBUG
#define assertpr(b, fmt, ...)                            \
  if (!(b)) {                                            \
    fprintf(stderr, "%s (%d): (" #b "): " fmt "\n",      \
            __FILE__, __LINE__, __VA_ARGS__);            \
    exit(-1); }
#else
#define assertpr(...)
#endif

  class Timer {
  public:
    static void GetTime (timeval& t) {
#ifdef ANALYZE_CODE
      gettimeofday(&t, 0);
#endif
    }      

    static double Et (const timeval& t1, const timeval& t2) {
#ifdef ANALYZE_CODE
      static const double us = 1.0e6;
      return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
#else
      return 0.0;
#endif
    }

    timeval Time () {
      timeval t;
#ifdef ANALYZE_CODE
      gettimeofday(&t, 0);
#else
      t = 0;
#endif
      return t;
    }

    void Reset (size_t idx) {
#ifdef ANALYZE_CODE
      if (idx >= _tv.size()) _tv.resize(idx + 1);
      _tv[idx].tot_et = 0.0;
#endif
    }

    void Tic (size_t idx) {
#ifdef ANALYZE_CODE
      if (idx >= _tv.size()) _tv.resize(idx + 1);
      gettimeofday(&_tv[idx].tv, 0);
#endif
    }

    double Toc (size_t idx) {
#ifdef ANALYZE_CODE
      if (idx >= _tv.size()) return 0.0;
      timeval t;
      gettimeofday(&t, 0);
      double et = Et(_tv[idx].tv, t);
      _tv[idx].tot_et += et;
      return et;
#else
      return 0.0;
#endif
    }

    double TotEt (size_t idx) {
#ifdef ANALYZE_CODE
      if (idx >= _tv.size()) Reset(idx);
      return _tv[idx].tot_et;
#else
      return 0.0;
#endif
    }

  private:
    struct Slot {
      timeval tv;
      double tot_et;
      Slot() : tot_et(0.0) {}
    };

    vector<Slot> _tv;
  };

  namespace Ca {
    util::Timer* GetTimer();
  }
}

/* Some handy debug code:
#include <sys/resource.h>
  static void ReportStuff()
  {
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    errpr("%1.1f %1.1f\n",
          (r.ru_utime.tv_sec*1.0e6 + r.ru_utime.tv_usec)*1.0e-6,
          r.ru_maxrss*1.0e-3);
  }
*/

#ifndef ANALYZE_CODE
#  undef timeval
#endif

#endif
