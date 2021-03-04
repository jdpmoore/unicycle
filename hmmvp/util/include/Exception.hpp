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

#ifndef INCLUDE_UTIL_EXCEPTION
#define INCLUDE_UTIL_EXCEPTION

#include <string>

namespace util {

  class Exception {
  public:
    Exception(const std::string& msg = "e") { _msg = msg; }
    virtual ~Exception() {}
    const std::string& GetMsg() const { return _msg; }
  private:
    std::string _msg;
  };

  class FileException : public Exception {
  public:
    FileException(const std::string& msg = "fe") : Exception(msg) {}
  };

  class OutOfMemoryException : public Exception {
  public:
    OutOfMemoryException(const std::string& msg = "oom") : Exception(msg) {}
  };

  class UserReqException : public Exception {
  public:
    UserReqException(const std::string& msg = "ur") : Exception(msg) {}
  };

}

#endif
