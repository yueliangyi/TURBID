/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <sys/stat.h>
#include <iomanip>
#include "boundary_scalar_flat.h"
#include "variable_scalar_flat.h"

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
FlatScalarBoundary::FlatScalarBoundary(const mpi_t& mpi,
                                         const dcp_t& dcp,
                                         Domain      * domain,
                                         SpaceStepper* sstepper,
                                         TimeStepper * tstepper,
                                         Variable    * variable,
                                         const int index) :
  Boundary(mpi,dcp,domain,sstepper,tstepper,variable,index) {
  PS_DEBUG_TRACE_ENTER(name_+"FlatScalarBoundary(mpi_t,dcp_t)")
    
//    const int rank = mpi_.GetRank();
//    const imat& dimension = dcp_.GetDimension();
//    const int sizexyz = dimension(10*3+7,rank)*dimension(10*3+8,rank)*3;
//    
//    dpdxyzm1top    = arma::zeros<cmat>(1,sizexyz);
//    dpdxyzm1bottom = arma::zeros<cmat>(1,sizexyz);
//    dpdxyzm1top.fill(cx_double(0,0));
//    dpdxyzm1bottom.fill(cx_double(0,0));
//    
//    dpdxyzm2top    = dpdxyzm1top;
//    dpdxyzm2bottom = dpdxyzm1bottom;

  PS_DEBUG_TRACE_LEAVE(name_+"FlatScalarBoundary(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarBoundary::Save(const string subfolder, const ioctrl_t level) {
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
bool FlatScalarBoundary::Load(const string subfolder, const ioctrl_t level) {
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
bool FlatScalarBoundary::Update(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Update(Domain,TimeStepper,SpaceStepper,Variable,int)")
  
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
bool FlatScalarBoundary::Apply(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Apply(Domain,TimeStepper,SpaceStepper,Variable,int)")
  
  FlatScalarVariable* pvariable = (FlatScalarVariable*)variable_;
  const size_t itop = GetDcp().GetSize(2)-1;
  const size_t ibottom = 0;
  
  pvariable->rhsphi->GetCMatCplxJoin().row(ibottom).fill(cx_double(0,0)*BOUNDARY_COEF_SCALE);
  pvariable->rhsphi->GetCMatCplxJoin().row(itop   ).fill(cx_double(0,0)*BOUNDARY_COEF_SCALE);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Apply(Domain,TimeStepper,SpaceStepper,Variable,int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Configuration
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
bool FlatScalarBoundary::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  Boundary::Configure();
  
//  dpdxyzm2top    = dpdxyzm1top;
//  dpdxyzm2bottom = dpdxyzm1bottom;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Print current status
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream as container
//------------------------------------------------------------------------------
void FlatScalarBoundary::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FlatScalarBoundary > ";
  Boundary::Print(stream);
  stream << std::endl;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
