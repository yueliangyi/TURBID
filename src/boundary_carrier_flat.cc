/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <sys/stat.h>
#include <iomanip>
#include "boundary_carrier_flat.h"
#include "variable_carrier_flat.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
// Function : Constructor
//------------------------------------------------------------------------------
// Parameter: IN  - mpi       > mpi object
//            IN  - dcp      > domain decomposition
//------------------------------------------------------------------------------
FlatCarrierBoundary::FlatCarrierBoundary(const mpi_t& mpi,
                                         const dcp_t& dcp,
                                         Domain      * domain,
                                         SpaceStepper* sstepper,
                                         TimeStepper * tstepper,
                                         Variable    * variable,
                                         const int index) :
  Boundary(mpi,dcp,domain,sstepper,tstepper,variable,index) {
  PS_DEBUG_TRACE_ENTER(name_+"FlatCarrierBoundary(mpi_t,dcp_t)")
    
    const int rank = mpi_.GetRank();
    const imat& dimension = dcp_.GetDimension();
    const int sizexyz = dimension(10*3+7,rank)*dimension(10*3+8,rank)*3;
    
    dp1dxyzm1top    = arma::zeros<arma::Row<cx_double>>(sizexyz);
    dp2dxyzm1top    = arma::zeros<arma::Row<cx_double>>(sizexyz);
    dp1dxyzm1bottom = arma::zeros<arma::Row<cx_double>>(sizexyz);
    dp2dxyzm1bottom = arma::zeros<arma::Row<cx_double>>(sizexyz);
    
    dp1dxyzm1top.fill(cx_double(0,0));
    dp2dxyzm1top.fill(cx_double(0,0));
    dp1dxyzm1bottom.fill(cx_double(0,0));
    dp2dxyzm1bottom.fill(cx_double(0,0));
    
    dp1dxyzm2top = dp1dxyzm1top;
    dp2dxyzm2top = dp2dxyzm1top;
    dp1dxyzm2bottom = dp1dxyzm1bottom;
    dp2dxyzm2bottom = dp2dxyzm1bottom;
    
  PS_DEBUG_TRACE_LEAVE(name_+"FlatCarrierBoundary(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierBoundary::Save(const string subfolder, const ioctrl_t level) {
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
bool FlatCarrierBoundary::Load(const string subfolder, const ioctrl_t level) {
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
// Function : The update step for the application boundary condition.
//            Usually, this function is called in the correction step,
//            collecting information on boundaries mainly. Still, for the
//            carrier phase, pressure gradient in z-direction is modified here.
//------------------------------------------------------------------------------
// Parameter: IN  - domain    > domain object
//            IN  - tstepper  > time stepper
//            IN  - sstepper  > space stepper
//            IN  - variable  > variable object
//            IN  - stage     > number of stage
//------------------------------------------------------------------------------
bool FlatCarrierBoundary::Update(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Update(Domain,TimeStepper,SpaceStepper,Variable,int)")
  
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  const int index = pvariable->GetIndex();
  const size_t itop = GetDcp().GetSize(2)-1;
  const size_t ibottom = 0;
  
  
  if (stage == 2) {
    olddt_ = tstepper_->GetDt();
  }
  
  
  // 1. Update variable (m2) from (m1).
  dp1dxyzm2top = dp1dxyzm1top;
  dp2dxyzm2top = dp2dxyzm1top;
  dp1dxyzm2bottom = dp1dxyzm1bottom;
  dp2dxyzm2bottom = dp2dxyzm1bottom;
  
  
  // 3. Always set dpdz=0 on top and bottom boundaries.
  // boundary condition for pressre is moved to the solving process.
//  // following part is added for better performance.
//  pvariable->dpdxyz[index]->GetCMatCplx(2).row(ibottom).fill(cx_double(0,0));
//  pvariable->dpdxyz[index]->GetCMatCplx(2).row(itop   ).fill(cx_double(0,0));
  
  
  // 2. Collect dpdx and dpdy on top and bottom boundaries.
  dp1dxyzm1bottom = pvariable->dpdxyz[index]->GetCMatCplxJoin().row(ibottom);
  dp1dxyzm1top    = pvariable->dpdxyz[index]->GetCMatCplxJoin().row(itop);
  if (isopentop_) {
    dp2dxyzm1bottom = pvariable->dudxyztmp[index]->GetCMatCplxJoin().row(ibottom);
    dp2dxyzm1top    = pvariable->dudxyztmp[index]->GetCMatCplxJoin().row(itop);
  }
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Update(Domain,TimeStepper,SpaceStepper,Variable,int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : The apply step for the application boundary condition.
//            Usually, this function is called in the prediction step,
//            specifying boundary values for intermediate velocity.
//------------------------------------------------------------------------------
// Parameter: IN  - domain    > domain object
//            IN  - tstepper  > time stepper
//            IN  - sstepper  > space stepper
//            IN  - variable  > variable object
//            IN  - stage     > number of stage
//------------------------------------------------------------------------------
bool FlatCarrierBoundary::Apply(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Apply(Domain,TimeStepper,SpaceStepper,Variable,int)")
  
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  const double dt   = tstepper_->GetDt();
  const double c3m0 = tstepper_->GetCoefficients(2,stage      );
  const double c3m1 = tstepper_->GetCoefficients(2,(stage+2)%3);
  const double c3m2 = tstepper_->GetCoefficients(2,(stage+1)%3);
  const size_t itop = GetDcp().GetSize(2)-1;
  const size_t ibottom = 0;
  const int rank = mpi_.GetRank();
  const imat& dimension = dcp_.GetDimension();
  const int sizexyz = dimension(10*3+7,rank)*dimension(10*3+8,rank);
  
  double xm0 = 0;
  double xm1 = 0;
  double xm2 = 0;
  
  switch (stage) {
      
    case 0:
      xm0 = 2*(0+c3m0)*dt;
      xm1 = 0;
      xm2 = 2*(0-c3m1)*olddt_;
      break;
      
    case 1:
      xm0 = 2*(c3m0+c3m1+0)*dt;
      xm1 = 2*(0   +c3m1+0)*dt;
      xm2 = 0;
      break;
      
    case 2:
      xm0 = 2*(c3m0+c3m1+c3m2)*dt;
      xm1 = 2*(0   +c3m1+c3m2)*dt;
      xm2 = 2*(0+  +0   +c3m2)*dt;
      break;
      
    default:
      break;
      
  }
  
  const double intepscale = (xm0-xm1)/(xm1-xm2);
  
  pvariable->rhsuvw->GetCMatCplxJoin().row(itop) =
  2*dt*c3m0*((1+intepscale)*dp1dxyzm1top   +(0-intepscale)*dp1dxyzm2top);
  pvariable->rhsuvw->GetCMatCplxJoin().row(ibottom) =
  2*dt*c3m0*((1+intepscale)*dp1dxyzm1bottom+(0-intepscale)*dp1dxyzm2bottom);
  
  if (isopentop_) {
    pvariable->rhsuvw->GetCMatCplx(0).row(itop) =
    2*dt*c3m0*((1+intepscale)*dp2dxyzm1top.subvec(0,sizexyz-1)+
               (0-intepscale)*dp2dxyzm2top.subvec(0,sizexyz-1));
    pvariable->rhsuvw->GetCMatCplx(1).row(itop) =
    2*dt*c3m0*((1+intepscale)*dp2dxyzm1top.subvec(sizexyz,2*sizexyz-1)+
               (0-intepscale)*dp2dxyzm2top.subvec(sizexyz,2*sizexyz-1));
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"Apply(Domain,TimeStepper,SpaceStepper,Variable,int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Configuration
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
bool FlatCarrierBoundary::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  Boundary::Configure();
  
  dp1dxyzm2top = dp1dxyzm1top;
  dp2dxyzm2top = dp2dxyzm1top;
  dp1dxyzm2bottom = dp1dxyzm1bottom;
  dp2dxyzm2bottom = dp2dxyzm1bottom;
  
  olddt_ = tstepper_->GetDt();
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Print current status
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream as container
//------------------------------------------------------------------------------
void FlatCarrierBoundary::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FlatCarrierBoundary > ";
  Boundary::Print(stream);
  stream << std::endl;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
