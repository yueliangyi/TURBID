/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_FILE_INI_H_
#define YLY_PS_LVBASE_FILE_INI_H_

#include <list>
#include <map>
#include "base.h"
#include "exprtk.hpp"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline string TrimTwoEnds(const string& target, const string& delimiter = " \f\n\r\t\v") {
  string result = target.substr(target.find_first_not_of(delimiter));
  result = result.substr(0,result.find_last_not_of(delimiter)+1);
  return result;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline string TrimComment(const string& target, const string& delimiter = ";#") {
  int pleft = target.find_first_of(delimiter);
  if (pleft == (int)string::npos) { return target; }
  return target.substr(0,pleft);
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline string ConvertCase(const string& target, const int issensitive) {
  string result = target;
  if (issensitive > 0) {
    std::transform(result.begin(),result.end(),result.begin(),::toupper);
  } else if (issensitive < 0) {
    std::transform(result.begin(),result.end(),result.begin(),::tolower);
  }
  return result;
}
  
  
  
class IniFile : public File {
  
public:
  
  IniFile(const mpi_t& mpi,
          const dcp_t& dcp,
          const string filename = IO_DEFAULT_FILE_NAME,
          const string directory = IO_DEFAULT_DIRECTORY,
          const int issensitive = 0);
  
  // Get different types of key values.
  string GetValueString (const string& sectionname, const string& keyname);
  int    GetValueInt    (const string& sectionname, const string& keyname);
  bool   GetValueBool   (const string& sectionname, const string& keyname);
  double GetValueDouble (const string& sectionname, const string& keyname);
  ivec   GetVectorInt   (const string& sectionname, const string& keyname);
  dvec   GetVectorDouble(const string& sectionname, const string& keyname);
  
  
  int  isSensitive(void) const { return issensitive_; }
  void SetSensitive(const int issensitive) { issensitive_ = issensitive; }
  
  bool Load(const string subfolder = "", const ioctrl_t level = IOControlLevelOne);
  
protected:
  
  struct ConfigurationSection {
    string name;
    std::map<string,string> keypair;
  }; typedef ConfigurationSection csec_t;
  
  // Synchronize information among all processors.
  bool Synchronize(void);

  bool AddSection(const string& sectionname);
  bool AddKey(const string& sectionname, const string& keyname, const string& keyvalue);
  int  FindSection(const string& sectionname) const;
  int  FindKey(const int& sectionid, const string& keyname) const;
  
  // Disable the public save operation.
  bool Save(const string subfolder, const ioctrl_t level) { return true; }
  
  std::list<csec_t> section_;
  int issensitive_;
  
  exprtk::symbol_table<double> symboltable_;
  exprtk::expression<double> expression_;
  exprtk::parser<double> parser_;
  
  string name_ = "IniFile::";

private:
  
};
  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_FILE_INI_H_
