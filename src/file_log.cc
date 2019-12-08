/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "file_log.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
LogFile::LogFile(const mpi_t& mpi,
                 const dcp_t& dcp,
                 const string filename,
                 const string directory,
                 const bool isprint,
                 const bool iswrite) :
  File(mpi,dcp,filename,directory),isprint_(isprint),iswrite_(iswrite) {
  PS_DEBUG_TRACE_ENTER(name_+"LogFile(mpi_t,dcp_t,string,string,bool,bool)")
    
    
  PS_DEBUG_TRACE_LEAVE(name_+"LogFile(mpi_t,dcp_t,string,string,bool,bool)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
LogFile::~LogFile(void) {
  
  if (mpi_.isMain()) {
    
    if (count_ != 0) {
      if (isprint_) { std::cout << message_.str(); }
      if (iswrite_) { filestream_ << message_.str(); filestream_.close(); }
      message_.str("");
      message_.clear();
      count_ = 0;
    }
  }
  
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool LogFile::Save(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Save(string,ioctrl_t)")
  
  
  if (mpi_.isMain()) {
    ++count_;
    if (count_ > IO_FILE_LOG_WRITE_INTERVAL) {
      if (isprint_) { std::cout << message_.str(); }
      if (iswrite_) { filestream_ << message_.str(); }
      message_.str("");
      message_.clear();
      count_ = 1;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Save(string,ioctrl_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool LogFile::Configure(const bool isprint, const bool iswrite) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(bool,bool)")
  
  
  isprint_ = isprint;
  iswrite_ = iswrite;
  
  
  if (mpi_.isMain() && iswrite_) {
    
    filestream_.open(directory_+filename_,fstream::out|fstream::app);
    
    if (!filestream_.is_open()) {
      /* ok, proceed with output */
    }
    
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(bool,bool)")
  return true;
}
  
  
  

  
  
  
  
  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi
