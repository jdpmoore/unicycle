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

#ifndef INCLUDE_UTIL_KEYVALUEFILE
#define INCLUDE_UTIL_KEYVALUEFILE

#include "util/include/Exception.hpp"
#include "util/include/Matrix.hpp"

// A simple means to create files containing key-value pairs for input and
// output encapsulation.

namespace util {
using namespace std;

typedef Matrix<double> Matd;

class KeyValueFile {
private:
  virtual ~KeyValueFile();

public:
  enum ValueType { vt_none, vt_string, vt_Matd };

  // Add values of the types: string, matrix, and scalar.

  void AddString(const string& key, const string& str);
  void AddMatd(const string& key, const Matd& m);
  // Convenience routine to add a 1x1 matrix.
  void AddDouble(const string& key, double d);

  // Determine the type of the value associated with the key.
  ValueType GetType(const string& key) const;

  // Get the values.
  bool GetString(const string& key, const string*& s) const;
  bool GetMatd(const string& key, const Matd*& m) const;
  // Convenience routine to get a 1x1 matrix. Fails if the matrix isn't 1x1.
  bool GetDouble(const string& key, double& d) const;
    
  // Write the key-value data to a file having a name ...
  bool Write(const string& filename) const;
  // ... or to a file output stream.
  bool Write(ofstream& os) const;

  // Read the contents of a key-value file.
  bool Read(const string& filename);
  bool Read(ifstream& is);

private:
  KeyValueFile();
  KeyValueFile(const KeyValueFile&);
  KeyValueFile& operator=(const KeyValueFile&);
};

// Create a KeyValueFile object.
KeyValueFile* NewKeyValueFile();
// Destroy it.
void DeleteKeyValueFile(KeyValueFile* kvf);
}

#endif
