/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <iomanip>
#include "file_ini.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
//string trim_right_copy(const string& str, const string& dlter = INI_DELIMITER) {
//  return str.substr(0,str.find_last_not_of(dlter)+1);
//}
//string trim_left_copy(const string& str, const string& dlter = INI_DELIMITER) {
//  return str.substr(str.find_first_not_of(dlter));
//}
//string trim_copy(const string& str, const string& dlter = INI_DELIMITER) {
//  return trim_left_copy(trim_right_copy(str,dlter),dlter);
//}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
IniFile::IniFile(const mpi_t& mpi, const dcp_t& dcp, const string filename,
                 const string directory, const int issensitive) :
  File(mpi,dcp,filename,directory),issensitive_(issensitive) {
  PS_DEBUG_TRACE_ENTER(name_+"IniFile(mpi_t,dcp_t,string,string,bool)")
    symboltable_.add_constants();
    expression_.register_symbol_table(symboltable_);
  PS_DEBUG_TRACE_LEAVE(name_+"IniFile(mpi_t,dcp_t,string,string,bool)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool IniFile::Load(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Load(string,ioctrl_t)")
  
  int failurelocal = 0;
  int failureglobal = 0;
  
  if (mpi_.isMain()) {
    filestream_.open(directory_+filename_,fstream::in);
    if (!filestream_.is_open()) { ++failurelocal; }
  }
  
  // Check whether open the configuration file successfully.
  MPI_Allreduce(&failurelocal,&failureglobal,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if (failureglobal != 0) {
    name_ += "Load(string,ioctrl_t) > Can not Open the Configuration File!";
    throw std::runtime_error(name_);
  }
  
  section_.clear();
  string line;
  string sectionname;
  string keyname;
  string keyvalue;
  while(mpi_.isMain() && getline(filestream_,line)) {
    
    int pleft = 0;
    int pright = 0;
    
    // Convert string case if required.
    line = ConvertCase(line,issensitive_);
    
    // Get rid of uninteresting lines.
    pleft = line.find_first_of(";#[=");
    if (pleft == (int)string::npos) { continue; }
    
    switch (line[pleft]) {
        
      case '[':
        pright = line.find_last_of("]");
        if (pright==(int)string::npos || pright<pleft) { break; }
        sectionname = TrimTwoEnds(line.substr(pleft+1,pright-pleft-1));
        AddSection(sectionname);
        break;
        
      case '=':
        keyname = TrimTwoEnds(line.substr(0,pleft));
        keyvalue = TrimComment(line.substr(pleft+1));
        if (keyvalue.find_first_not_of(' ') == string::npos) { keyvalue = ""; }
        else { keyvalue = TrimTwoEnds(keyvalue); }
        AddKey(sectionname,keyname,keyvalue);
        break;
        
      case ';':
      case '#':
      default:
        // Ignore all comments and other unknown symbols.
        break;
        
    }
    
  }
  
  if (mpi_.isMain()) { filestream_.close(); }
  MPI_Barrier(MPI_COMM_WORLD);
  Synchronize();
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Load(string,ioctrl_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool IniFile::Synchronize(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Synchronize(void)")
  
  int sectionnum = (int)section_.size();
  mpi_.BcastInt(sectionnum);
  if (!mpi_.isMain()) { section_.resize(sectionnum); }
  
  for (auto it = section_.begin(); it != section_.end(); ++it) {
    
    // Fill in section name.
    string sectionname = it->name;
    mpi_.BcastString(sectionname);
    if (!mpi_.isMain()) { it->name = sectionname; }
    
    // Fill in key pair.
    int keynum = (int)it->keypair.size();
    mpi_.BcastInt(keynum);
    if (keynum == 0) { continue; }
    for (int index = 0; index < keynum; ++index) {
      auto targetkey = it->keypair.begin();
      string keyname;
      string keyvalue;
      if (mpi_.isMain()) {
        advance(targetkey,index);
        keyname = targetkey->first;
        keyvalue = targetkey->second;
      }
      mpi_.BcastString(keyname);
      mpi_.BcastString(keyvalue);
      if (mpi_.isMain()) { continue; }
      it->keypair.insert(std::pair<string,string>(keyname,keyvalue));
    }
    
  }
  
  
  
  
//  for (auto it = section_.cbegin(); it != section_.cend(); ++it) {
//    // Section title.
//    _message_
//    << std::endl
//    << std::endl
//    << '['
//    << it->name
//    << ']'
//    << std::endl;
//    // Key and corresponding value.
//    for (auto it_key = it->keypair.cbegin(); it_key != it->keypair.cend(); ++it_key) {
//      _message_
//      << std::left
//      << std::setw(25)
//      << it_key->first
//      << " = "
//      << it_key->second
//      << std::endl;
//    }
//  }
  
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Synchronize(void)")
  return true;
}
  
    
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
string IniFile::GetValueString(const string& sectionname, const string& keyname) {
  PS_DEBUG_TRACE_ENTER(name_+"GetValueString(string,string)")
  
  int sectionid = FindSection(sectionname);
  int keyid = FindKey(sectionid,keyname);
  
  if (sectionid < 0) {
    name_ += "GetValueString(string,string) > Can not Find Section: "+sectionname;
    throw std::runtime_error(name_);
  }
  if (keyid < 0) {
    name_ += "GetValueString(string,string) > Can not Find Key: ";
    name_ += keyname+" in Section "+sectionname;
    throw std::runtime_error(name_);
  }
  
  auto targetsection = section_.cbegin();
  advance(targetsection,sectionid);
  auto targetkey = targetsection->keypair.begin();
  advance(targetkey,keyid);
  
  PS_DEBUG_TRACE_LEAVE(name_+"GetValueString(string,string)")
  return targetkey->second;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int IniFile::GetValueInt(const string& sectionname, const string& keyname) {
  PS_DEBUG_TRACE_ENTER(name_+"GetValueInt(string,string)")
  PS_DEBUG_TRACE_LEAVE(name_+"GetValueInt(string,string)")
  return std::stoi(GetValueString(sectionname,keyname));
}
    
    
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool IniFile::GetValueBool(const string& sectionname, const string& keyname) {
  PS_DEBUG_TRACE_ENTER(name_+"GetValueBool(string,string)")
  PS_DEBUG_TRACE_LEAVE(name_+"GetValueBool(string,string)")
  return GetValueInt(sectionname,keyname);
}
    
    
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double IniFile::GetValueDouble(const string& sectionname, const string& keyname) {
  PS_DEBUG_TRACE_ENTER(name_+"GetValueDouble(string,string)")
  string result = GetValueString(sectionname,keyname);
  parser_.compile(result,expression_);
  PS_DEBUG_TRACE_LEAVE(name_+"GetValueDouble(string,string)")
  return expression_.value();
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
ivec IniFile::GetVectorInt(const string& sectionname, const string& keyname) {
  PS_DEBUG_TRACE_ENTER(name_+"GetVectorInt(string,string)")
  
  std::istringstream stream(GetValueString(sectionname,keyname));
  std::istream_iterator<std::string> begin(stream),end;
  std::vector<std::string> tokens(begin,end);
  
  int size = (int)tokens.size();
  ivec result(size);
  for (int id = 0; id < size; ++id) {
    parser_.compile(tokens[id],expression_);
    result(id) = expression_.value();
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"GetVectorInt(string,string)")
  return result;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
dvec IniFile::GetVectorDouble(const string& sectionname, const string& keyname) {
  PS_DEBUG_TRACE_ENTER(name_+"GetVectorDouble(string,string)")
  
  std::istringstream stream(GetValueString(sectionname,keyname));
  std::istream_iterator<std::string> begin(stream),end;
  std::vector<std::string> tokens(begin,end);
  
  int size = (int)tokens.size();
  dvec result(size);
  for (int id = 0; id < size; ++id) {
    parser_.compile(tokens[id],expression_);
    result(id) = expression_.value();
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"GetVectorDouble(string,string)")
  return result;
}
  
    
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool IniFile::AddSection(const string& sectionname) {
  PS_DEBUG_TRACE_ENTER(name_+"AddSection(string)")
  csec_t section;
  section.name = sectionname;
  section_.push_back(section);
  PS_DEBUG_TRACE_LEAVE(name_+"AddSection(string)")
  return true;
}
    
    
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool IniFile::AddKey(const string& sectionname, const string& keyname, const string& keyvalue) {
  PS_DEBUG_TRACE_ENTER(name_+"TrimTwoEnds(string,string)")
  
  // Insert a new section if necessary.
  int sectionid = FindSection(sectionname);
  if (sectionid < 0) {
    AddSection(sectionname);
    sectionid = FindSection(sectionname);
  }
  
  auto targetsection = section_.begin();
  advance(targetsection,sectionid);
  
  // Insert a new key if necessary.
  // Otherwise, replace the corresponding key value.
  int keyid = FindKey(sectionid,keyname);
  if (keyid < 0) {
    targetsection->keypair.insert(std::pair<string,string>(keyname,keyvalue));
  } else {
    auto targetkey = targetsection->keypair.begin();
    advance(targetkey,keyid);
    targetkey->second = keyvalue;
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"TrimTwoEnds(string,string)")
  return true;
}
    
    
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int IniFile::FindSection(const string& sectionname) const {
  PS_DEBUG_TRACE_ENTER(name_+"FindSection(string)")
  int sectionid = -1;
  for (auto it = section_.cbegin(); it != section_.cend(); it++) {
    if (it->name == sectionname) {
      sectionid = std::distance(section_.cbegin(),it);
    }
  }
  PS_DEBUG_TRACE_LEAVE(name_+"FindSection(string)")
  return sectionid;
}
    

  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int IniFile::FindKey(const int& sectionid, const string& keyname) const {
  PS_DEBUG_TRACE_ENTER(name_+"FindKey(string,string)")
  int keyid = -1;
  if (section_.size()==0 || sectionid>=(int)section_.size()) { return keyid; }
  auto targetsection = section_.cbegin();
  advance(targetsection,sectionid);
  auto it = targetsection->keypair.find(keyname);
  if(it != targetsection->keypair.end()) {
    keyid = std::distance(targetsection->keypair.cbegin(),it);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"FindKey(string,string)")
  return keyid;
}
  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi
