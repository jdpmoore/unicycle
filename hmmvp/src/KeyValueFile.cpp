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

#include <assert.h>
#include <fstream>
#include <vector>
#include "util/include/Defs.hpp"
#include "KeyValueFile_pri.hpp"
using namespace std;

namespace util {

typedef int64 kvf_int;

KeyValueFile* NewKeyValueFile () { return new KeyValueFile(); }

void DeleteKeyValueFile (KeyValueFile* kvf) { delete kvf; }

template<typename T>
static void FreeMap (map<string, T>& m) {
  for (typename map<string, T>::iterator it = m.begin(), end = m.end();
       it != end; ++it)
    delete it->second;
}

template<typename T>
static inline void EraseIfNeeded (map<string, T>& m, const string& key) {
  typename map<string, T>::iterator it = m.find(key);
  if (it != m.end()) {
    delete it->second;
    m.erase(it);
  }
}

template<typename T>
static inline void Insert (map<string, T*>& m, const string& key, T* v) {
  m.insert(pair<string, T*>(key, v));
}

template<typename T> static inline void
CopyAndInsert (map<string, T*>& m, const string& key, const T& v) {
  T* vc = new T(v);
  Insert(m, key, vc);
}

KeyValueFile::~KeyValueFile () {
  FreeMap(_s_string);
  FreeMap(_s_Matd);
}

void KeyValueFile::AddType (const string& key, ValueType vt) {
  map<string, ValueType>::iterator it = _s_vt.find(key);
  if (it != _s_vt.end()) _s_vt.erase(it);
  _s_vt.insert(pair<string, ValueType>(key, vt));
}

void KeyValueFile::AddString (const string& key, const string& str) {
  AddString(key, new string(str));
}

void KeyValueFile::AddString (const string& key, string* str) {
  EraseIfNeeded(_s_string, key);
  Insert(_s_string, key, str);
  AddType(key, vt_string);
}

void KeyValueFile::AddMatd (const string& key, const Matd& m) {
  AddMatd(key, new Matd(m));
}

void KeyValueFile::AddDouble (const string& key, double d) {
  Matd* s = new Matd(1, 1);
  (*s)(1,1) = d;
  AddMatd(key, s);
}

void KeyValueFile::AddMatd (const string& key, Matd* m) {
  EraseIfNeeded(_s_Matd, key);
  Insert(_s_Matd, key, m);
  AddType(key, vt_Matd);
}

KeyValueFile::ValueType KeyValueFile::GetType (const string& key) const {
  map<string, ValueType>::const_iterator it = _s_vt.find(key);
  if (it == _s_vt.end()) return vt_none;
  return it->second;
}

template<typename T> static inline bool
GetValue (const map<string, T*>& m, const string& key, const T*& v) {
  typename map<string, T*>::const_iterator it = m.find(key);
  if (it == m.end()) {
    v = NULL;
    return false;
  }
  v = it->second;
  return true;
}

bool KeyValueFile::GetString (const string& key, const string*& s) const {
  return GetValue(_s_string, key, s);
}

bool KeyValueFile::GetMatd (const string& key, const Matd*& v) const {
  return GetValue(_s_Matd, key, v);
}

bool KeyValueFile::GetDouble (const string& key, double& d) const {
  const Matd* m;
  if (!GetValue(_s_Matd, key, m) || m->Size() > 1) return false;
  d = (*m)(1);
  return true;
}

static bool WriteCode (ofstream& os, const char* code) {
  return !os.write(code, 2).bad();
}

static bool ReadCode (ifstream& is, char* code) {
  return !is.read(code, 2).bad();
}

static bool WriteString (ofstream& os, const string& s) {
  kvf_int n = s.length();
  return !(os.write((char*) &n, sizeof(kvf_int)).bad() ||
           os.write(s.c_str(), n).bad());
}

static bool ReadStringAndTestEof (ifstream& is, string& s, vector<char>& work) {
  kvf_int n;
  if (is.read((char*) &n, sizeof(kvf_int)).bad()) return false;
  if (is.eof()) return false;
  work.resize(n + 1);
  if (is.read(&work[0], n).bad()) return false;
  s.resize(0);
  s.insert(0, &work[0], n);
  return true;
}

static bool WriteInts (ofstream& os, const vector<kvf_int>& vi) {
  kvf_int n = vi.size();
  return !(os.write((char*) &n, sizeof(kvf_int)).bad() ||
           os.write((char*) &vi[0], n*sizeof(kvf_int)).bad());
}

static bool ReadInts (ifstream&is, vector<kvf_int>& vi) {
  kvf_int n;
  if (is.read((char*) &n, sizeof(kvf_int)).bad()) return false;
  vi.resize(n);
  return (!is.read((char*) &vi[0], n*sizeof(kvf_int)).bad());
}

template<typename T>
static bool WriteArrayData (ofstream& os, const Matrix<T>& m) {
  kvf_int n = m.Size();
  return !os.write((char*) m.GetPtr(), n*sizeof(T)).bad();
}

template<typename T>
static bool ReadArrayData (ifstream& is, Matrix<T>& m) {
  return !is.read((char*) m.GetPtr(), m.Size()*sizeof(T)).bad();
}

bool KeyValueFile::Write (const string& filename) const {
  ofstream os(filename.c_str(), ofstream::binary);
  if (!os.is_open()) return false;
  bool ret = Write(os);
  os.close();
  return ret;
}
  
bool KeyValueFile::Write (ofstream& os) const {
  for (map<string, string*>::const_iterator it = _s_string.begin(),
         end = _s_string.end(); it != end; ++it) {
    if (!WriteString(os, it->first) ||
        !WriteCode(os, "st") ||
        !WriteString(os, *it->second))
      return false;
  }
    
  vector<kvf_int> sz(2);
  for (map<string, Matd*>::const_iterator it = _s_Matd.begin(),
         end = _s_Matd.end(); it != end; ++it) {
    sz[0] = it->second->Size(1);
    sz[1] = it->second->Size(2);
    if (!WriteString(os, it->first) ||
        !WriteCode(os, "da") ||
        !WriteInts(os, sz) ||
        !WriteArrayData(os, *it->second))
      return false;
  }

  return true;
}
  
bool KeyValueFile::Read (const string& filename) {
  ifstream is(filename.c_str(), ofstream::binary);
  if (!is.is_open()) return false;
  bool ret = Read(is);
  is.close();
  return ret;
}

static inline bool IsCode (const char c1[2], const char c2[2]) {
  return c1[0] == c2[0] && c1[1] == c2[1];
}
  
bool KeyValueFile::Read (ifstream& is) {
  vector<char> work;
  while (true) {
    string key;
    if (!ReadStringAndTestEof(is, key, work)) return is.eof();
      
    char code[2];
    if (!ReadCode(is, code)) return false;

    if (IsCode(code, "st")) {
      string s;
      if (!ReadStringAndTestEof(is, s, work)) return false;
      AddString(key, s);
    } else if (IsCode(code, "da")) {
      vector<kvf_int> sz;
      if (!ReadInts(is, sz)) return false;
      Matd* m;
      switch (sz.size()) {
      case 1:
        m = new Matd(sz[0]);
        break;
      case 2:
        m = new Matd(sz[0], sz[1]);
        break;
      default:
        int numel = 1;
        for (int i = 0, n = sz.size(); i < n; i++) numel *= sz[i];
        m = new Matd(numel);
      }
      ReadArrayData(is, *m);
      AddMatd(key, m);
    } else {
      return false;
    }
  }
}

};
