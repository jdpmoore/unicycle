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

#ifndef INCLUDE_UTIL_KEYVALUEFILE_PRI
#define INCLUDE_UTIL_KEYVALUEFILE_PRI

#include <map>
#include <string>
#include "util/include/Defs.hpp"
#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"

namespace util {
  using namespace std;

  typedef Matrix<double> Matd;

  class KeyValueFile {
  public:
    enum ValueType { vt_none, vt_string, vt_Matd };

    KeyValueFile() {};
    virtual ~KeyValueFile();

    void AddString(const string& key, const string& str);
    void AddMatd(const string& key, const Matd& m);
    void AddDouble(const string& key, double d);

    ValueType GetType(const string& key) const;
    bool GetString(const string& key, const string*& s) const;
    bool GetMatd(const string& key, const Matd*& m) const;
    bool GetDouble(const string& key, double& d) const;
    
    bool Write(const string& filename) const;
    bool Write(ofstream& os) const;

    bool Read(const string& filename);
    bool Read(ifstream& is);

  private:
    map<string, string*> _s_string;
    map<string, Matd*> _s_Matd;
    map<string, ValueType> _s_vt;

    void AddString(const string& key, string* str);
    void AddMatd(const string& key, Matd* m);
    void AddType(const string& key, ValueType vt);

    KeyValueFile(const KeyValueFile&);
    KeyValueFile& operator=(const KeyValueFile&);
  };

};

#endif
