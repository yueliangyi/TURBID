/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "stepper_time_rk3.h"
#include <iomanip>
#include <sys/stat.h>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
RK3TimeStepper::RK3TimeStepper(const mpi_t& mpi, const dcp_t& dcp, const Base::cfg_t* cfg) :
  TimeStepper(mpi,dcp,cfg) {
  PS_DEBUG_TRACE_ENTER("RK3TimeStepper::RK3TimeStepper(mpi_t,dcp_t,cfg_t)")
    
  // Coefficients for terms at different stages.
  // Refer corresponding report for details.
  // Column: the first, second nonlinear and diffusion term;
  // Row   : different stages.
  coefficients_.set_size(3,3);
  coefficients_(0,0) = 0.0/1.0; coefficients_(0,1) = -5.0/9.0 ; coefficients_(0,2) = -153.0/128.0;
  coefficients_(1,0) = 1.0/3.0; coefficients_(1,1) = 15.0/16.0; coefficients_(1,2) = 8.0   /15.0 ;
  coefficients_(2,0) = 1.0/6.0; coefficients_(2,1) = 5.0 /24.0; coefficients_(2,2) = 1.0   /8.0  ;
    
  PS_DEBUG_TRACE_LEAVE("RK3TimeStepper::RK3TimeStepper(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool RK3TimeStepper::Save(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Save(int,ioctrl_t,string)")
  
  // Check whether the target folder exists or not.
  // If not exist (the usual case), try to create one.
  struct stat status = {0};
  const mode_t mode = S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH;
  const string path = directory_+subfolder;
  if (mpi_.isMain() &&
      stat(path.c_str(),&status)==-1 &&
      mkdir(path.c_str(),mode)==-1) {
    name_ += "Save(int,ioctrl_t,string) > Failed to Create Target Folder!";
    throw std::runtime_error(name_);
  }
  
  if (mpi_.isMain()) {
    filestream_.open(path+filename_,fstream::out|fstream::binary);
  }
  
  
  
  
  
  
  if (mpi_.isMain()) {
    filestream_.close();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Save(int,ioctrl_t,string)")
  return true;
}
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool RK3TimeStepper::Load(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Load(int,ioctrl_t,string)")
  
  if (mpi_.isMain()) {
    filestream_.open(directory_+subfolder+filename_,fstream::in|fstream::binary);
  }
  
  
  
  
  
  if (mpi_.isMain()) {
    filestream_.close();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Load(int,ioctrl_t,string)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void RK3TimeStepper::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("RK3TimeStepper::Print(ostringstream)")
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "RK3TimeStepper > ";
  TimeStepper::Print(stream);
  stream << std::endl;
  PS_DEBUG_TRACE_LEAVE("RK3TimeStepper::Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
