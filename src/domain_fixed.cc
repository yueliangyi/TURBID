/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "domain_fixed.h"
#include <iomanip>
#include <sys/stat.h>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FixedDomain::FixedDomain(const mpi_t& mpi, const dcp_t& dcp, const Base::cfg_t* cfg) :
  Domain(mpi,dcp,cfg) {
  PS_DEBUG_TRACE_ENTER("FixedDomain::FixedDomain(mpi_t,dcp_t,cfg_t)")
  PS_DEBUG_TRACE_LEAVE("FixedDomain::FixedDomain(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FixedDomain::Save(const string subfolder, const ioctrl_t level) {
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
  
  
  
  // save top and bottom profiles first
  
  
  
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
bool FixedDomain::Load(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Load(int,ioctrl_t,string)")
  
  if (mpi_.isMain()) {
    filestream_.open(directory_+subfolder+filename_,fstream::in|fstream::binary);
  }
  
  
  // Only need to load top and bottom profiles
  // Others will be generated later
//  if (mpi_.isMain()) {
//    
//    profiletop_    = arma::zeros(size(0)+1,size(1)+1);
//    profilebottom_ = arma::zeros(size(0)+1,size(1)+1);
//    
//  }
  
  
  
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
void FixedDomain::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("FixedDomain::Print(ostringstream)")
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FixedDomain > ";
  Domain::Print(stream);
  stream << std::endl;
  PS_DEBUG_TRACE_LEAVE("FixedDomain::Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
