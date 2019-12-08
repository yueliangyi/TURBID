/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_FILE_LOG_H_
#define YLY_PS_LVBASE_FILE_LOG_H_

#include "base.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
class LogFile : public File {
  
public:
  
  
  LogFile(const mpi_t& mpi,
          const dcp_t& dcp,
          const string filename = IO_DEFAULT_FILE_NAME,
          const string directory = IO_DEFAULT_DIRECTORY,
          const bool isprint = true,
          const bool iswrite = true);
  ~LogFile(void);
  
  bool isPrint(void) { return isprint_; }
  
  const ostringstream& GetMessage(void) const { return message_; }
  ostringstream& GetMessage(void) { return message_; }
  
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  
  bool Configure(const bool isprint = true, const bool iswrite = true);
  
  
protected:
  
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne) { return true; }
  
  
  ostringstream message_;
  size_t count_ = 1;
  
  bool isprint_ = true;
  bool iswrite_ = true;
  
  
  string name_ = "LogFile::";

private:
  
};
  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_FILE_LOG_H_
