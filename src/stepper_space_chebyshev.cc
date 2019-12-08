/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "stepper_space_chebyshev.h"
#include <iomanip>
#include <sys/stat.h>
#include <math.h>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
ChebyshevSpaceStepper::ChebyshevSpaceStepper(const mpi_t& mpi,
                                             const dcp_t& dcp,
                                             Domain* domain,
                                             const Base::cfg_t* cfg) :
  SpaceStepper(mpi,dcp,domain,cfg) {
  PS_DEBUG_TRACE_ENTER("ChebyshevSpaceStepper::ChebyshevSpaceStepper(mpi_t,dcp_t,Domain,cfg_t)")
  PS_DEBUG_TRACE_LEAVE("ChebyshevSpaceStepper::ChebyshevSpaceStepper(mpi_t,dcp_t,Domain,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool ChebyshevSpaceStepper::Save(const string subfolder, const ioctrl_t level) {
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
bool ChebyshevSpaceStepper::Load(const string subfolder, const ioctrl_t level) {
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
bool ChebyshevSpaceStepper::Configure(void) {
  PS_DEBUG_TRACE_ENTER("ChebyshevSpaceStepper::Configure(void)")

  //
  SpaceStepper::Configure();
  
  
  const ivec& size = dcp_.GetSize();
  
  cmdz1_.set_size(size(2),size(2));
  cmdz2_.set_size(size(2),size(2));
  
  const dmat& gridz = domain_->GetGridZ();
  
  
  // Construct members.
  // For details, refer to the corresponding report.
  
  // Generate the coefficient matrix for the first z-derivative computation.
  // Keep using the formula cos()-cos() in denominator other than sin()*sin()
  // to reduce the round off error.
  for (int iz = 1; iz < size(2); ++iz) {
    for (int ix = 0; ix < iz; ++ix) {
      cmdz1_(iz,ix) = pow(-1,iz+ix)/(gridz(iz)-gridz(ix));
    }
  }
  // Fix the edges of matrix and fill in the upper triangle using antisymmetry.
  cmdz1_.col(0        ) *= 0.5;
  cmdz1_.row(size(2)-1) *= 2.0;
  for (int iz = 0; iz < size(2)-1; ++iz) {
    for (int ix = iz+1; ix < size(2); ++ix) {
      cmdz1_(iz,ix) = -cmdz1_(size(2)-iz-1,size(2)-ix-1);
    }
  }
  
  
  // Compute the diagonal terms from the other terms on the same line.
  // Reducing the round off error by following this way.
  for (int iz = 0; iz < size(2); ++iz) {
    cmdz1_(iz,iz) = 0;
    for (int ix = 0; ix < size(2); ++ix) {
      if (ix == iz) { continue; }
      cmdz1_(iz,iz) = cmdz1_(iz,iz)-cmdz1_(iz,ix);
    }
  }
  
  
  
  
  // adjust coefficient matrix for grid clustering
  cfg_t* pcfg = (cfg_t*)cfg_;
  dvec gridzcluster;
  switch (pcfg->mappingfunction) {
      
    case kTwoEndTangent:
    {
      const double ccoef = pcfg->mappingcoefficient(0);
      const dvec tantmp = arma::tan(arma::datum::pi/2*gridz);
      gridzcluster = 2/arma::datum::pi*arma::atan(ccoef*tantmp);
      const dvec dclusterdx3 = 1/ccoef*(1+arma::pow(ccoef*tantmp,2))/(1+arma::pow(tantmp,2));
      cmdz1_ = arma::diagmat(dclusterdx3)*cmdz1_;
    }
      break;
      
    case kTwoEndSin:
    {
      const double ccoef = pcfg->mappingcoefficient(0);
      gridzcluster = arma::asin(ccoef*gridz)/asin(ccoef);
      const dvec dclusterdx3 = 1/ccoef*arma::cos(arma::asin(ccoef*gridz))*asin(ccoef);
      cmdz1_ = arma::diagmat(dclusterdx3)*cmdz1_;
    }
      break;
      
    case kOneEndTangent:
    {
      const double ccoef = pcfg->mappingcoefficient(0);
      const dvec tantmp = arma::tan(arma::datum::pi/4*(gridz-1));
      gridzcluster = 4/arma::datum::pi*arma::atan(ccoef*tantmp)+1;
      const dvec dclusterdx3 = 1/ccoef*(1+arma::pow(ccoef*tantmp,2))/(1+arma::pow(tantmp,2));
      cmdz1_ = arma::diagmat(dclusterdx3)*cmdz1_;
    }
      break;
      
    default:
      break;
  }
  
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  if (pcfg->mappingfunction != kNone) {

    domain_->GetGridZCluster() = gridzcluster;
    dcube& gridlenz = domain_->GetGridLenZ();

    for (int id = 0; id < dcpdimension(10*2+8,rank); ++id) {
      if (dcpdimension(10*2+2,rank) == 1 && id == 0) {
        gridlenz.slice(id).fill((gridzcluster(1)-gridzcluster(0))/2);
        continue;
      }
      if (dcpdimension(10*2+5,rank) == size(2) && id == dcpdimension(10*2+8,rank)-1) {
        gridlenz.slice(id).fill((gridzcluster(size(2)-1)-gridzcluster(size(2)-2))/2);
        continue;
      }
      gridlenz.slice(id).fill(std::min(gridzcluster(dcpdimension(10*2+2,rank)-1+id+1)-gridzcluster(dcpdimension(10*2+2,rank)-1+id),
                                       gridzcluster(dcpdimension(10*2+2,rank)-1+id)-gridzcluster(dcpdimension(10*2+2,rank)-1+id-1)));
    }

  }
  
  
  

  
  
  
  
  
  // Generate the coefficient matrix for the second z-derivative computation.
  // The matrix comes from the first z-derivative differentiation matrix.
  cmdz2_ = cmdz1_*cmdz1_;
  for (int iz = 0; iz < size(2); ++iz) {
    cmdz2_(iz,iz) = 0;
    for (int ix = 0; ix < size(2); ++ix) {
      if (ix == iz) { continue; }
      cmdz2_(iz,iz) = cmdz2_(iz,iz)-cmdz2_(iz,ix);
    } // end for ix
  } // end for iz
  
  
  
//  if (mpi_.isMain()) {
//    dmat cmdz1 = arma::real(cmdz1_);
//    dmat cmdz2 = arma::real(cmdz2_);
//    cmdz1.save("cmdz1",arma::arma_ascii);
//    cmdz2.save("cmdz2",arma::arma_ascii);
//  }
  
  
  
  
//  tmz1_.set_size(size(2),size(2));
//  tmz2_.set_size(size(2),size(2));
  
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE("ChebyshevSpaceStepper::Configure(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void ChebyshevSpaceStepper::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("ChebyshevSpaceStepper::Print(ostringstream)")
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "ChebyshevSpaceStepper > ";
  SpaceStepper::Print(stream);
//  stream << std::endl;
  PS_DEBUG_TRACE_LEAVE("ChebyshevSpaceStepper::Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
