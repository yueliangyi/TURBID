/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <sys/stat.h>
#include <iomanip>
#include "variable_scalar_flat.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatScalarVariable::FlatScalarVariable(const mpi_t& mpi,
                                         const dcp_t& dcp,
                                         Domain      * domain,
                                         SpaceStepper* sstepper,
                                         TimeStepper * tstepper,
                                         const int index) :
  Variable(mpi,dcp,domain,sstepper,tstepper,index) {
  PS_DEBUG_TRACE_ENTER(name_+"FlatScalarVariable(mpi_t,dcp_t)")
    
    uvwnonlinear = new rccube(mpi_,dcp_,3);
    advphi = new rccube(mpi_,dcp_,1);
    hphi = new rccube(mpi_,dcp_,1,false);
    lapphi = new rccube(mpi_,dcp_,1,false);
    rhsphi = new rccube(mpi_,dcp_,1,false);
    
  PS_DEBUG_TRACE_LEAVE(name_+"FlatScalarVariable(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatScalarVariable::~FlatScalarVariable() {
  
  if (uvwnonlinear != NULL) { delete uvwnonlinear; uvwnonlinear = NULL; }
  if (advphi != NULL) { delete advphi; advphi = NULL; }
  if (hphi != NULL) { delete hphi; hphi = NULL; }
  if (lapphi != NULL) { delete lapphi; lapphi = NULL; }
  if (rhsphi != NULL) { delete rhsphi; rhsphi = NULL; }
  
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarVariable::Save(const string subfolder, const ioctrl_t level) {
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
  
  
  switch (level) {
      
    case IOControlLevelThree:
      filestream_ << "IOControlLevelThree" << std::endl;
      // Save extra/dependent variables.
//      advphi->DftBackward(); advphi->SaveReal(filestream_);
//      lapphi->DftBackward(); lapphi->SaveReal(filestream_);
//      rhsphi->DftBackward(); rhsphi->SaveReal(filestream_);
//      hphi  ->DftBackward(); hphi  ->SaveReal(filestream_);
      
    case IOControlLevelTwo:
      filestream_ << "IOControlLevelTwo" << std::endl;
      // Save partial derivatives.
//      dudxyz[index_]->DftBackward(); dudxyz[index_]->SaveReal(filestream_);
//      dvdxyz[index_]->DftBackward(); dvdxyz[index_]->SaveReal(filestream_);
//      dwdxyz[index_]->DftBackward(); dwdxyz[index_]->SaveReal(filestream_);
//      dpdxyz[index_]->DftBackward(); dpdxyz[index_]->SaveReal(filestream_);
      
    case IOControlLevelOne:
      filestream_ << "IOControlLevelOne" << std::endl;
      // Save basic variables.
      // Save hu/y/z is necessary.
//      uvw [index_]->DftBackward(); uvw [index_]->SaveReal(filestream_);
//      p [index_]->DftBackward(); p [index_]->SaveReal(filestream_);
      
      p[index_]  ->SaveCplx(filestream_);
      
      break;
      
    default:
      name_ += "Save(int,ioctrl_t,string) > Not Supported IO Level Type!";
      throw std::runtime_error(name_);
      break;
      
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
bool FlatScalarVariable::Load(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Load(int,ioctrl_t,string)")
  
  string levelname;
  
  int failurelocal = 0;
  int failureglobal = 0;
  
  if (mpi_.isMain()) {
    filestream_.open(directory_+subfolder+filename_,fstream::in|fstream::binary);
    if (!filestream_.is_open()) { ++failurelocal; }
    std::getline(filestream_,levelname);
  }
  
  // Check whether open the input file successfully.
  MPI_Allreduce(&failurelocal,&failureglobal,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if (failureglobal != 0) {
    name_ += "Load(string,ioctrl_t) > Can not Open the Input File!";
    throw std::runtime_error(name_);
  }
  
  
  
  
  switch (level) {
      
    case IOControlLevelThree:
      // Save extra/dependent variables.
      if (mpi_.isMain() && levelname!="IOControlLevelThree") {
        name_ += "Load(int,ioctrl_t,string) > Wrong Request of IO Level Three!";
        throw std::runtime_error(name_);
      }
//      advphi->LoadReal(filestream_); advphi->DftForward();
//      lapphi->LoadReal(filestream_); lapphi->DftForward();
//      rhsphi->LoadReal(filestream_); rhsphi->DftForward();
//      hphi  ->LoadReal(filestream_); hphi  ->DftForward();
      
    case IOControlLevelTwo:

      // Save partial derivatives.
      if (mpi_.isMain() && levelname!="IOControlLevelTwo") {
        name_ += "Load(int,ioctrl_t,string) > Wrong Request of IO Level Two!";
        throw std::runtime_error(name_);
      }
//      dudxyz[index_]->LoadReal(filestream_); dudxyz[index_]->DftForward();
//      dvdxyz[index_]->LoadReal(filestream_); dvdxyz[index_]->DftForward();
//      dwdxyz[index_]->LoadReal(filestream_); dwdxyz[index_]->DftForward();
      
      
    case IOControlLevelOne:

      // Save basic variables.
      // Save hu/y/z is necessary.
      if (mpi_.isMain() &&  levelname!="IOControlLevelOne") {
        name_ += "Load(int,ioctrl_t,string) > Wrong Request of IO Level One!";
        throw std::runtime_error(name_);
      }
//      uvw[index_]->LoadReal(filestream_); uvw[index_]->DftForward();
//      p[index_]  ->LoadReal(filestream_); p[index_]  ->DftForward();
      
      p[index_]  ->LoadCplx(filestream_);
      
      break;
      
    default:
      name_ += "Load(int,ioctrl_t,string) > Not Supported IO Level Type!";
      throw std::runtime_error(name_);
      break;
      
  }
  
  
  if (mpi_.isMain()) {
    filestream_.close();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  
  PS_DEBUG_TRACE_LEAVE(name_+"Load(int,ioctrl_t,string)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Configuration
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
bool FlatScalarVariable::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  
  Variable::Configure();
  
  
//  const int startstep = tstepper_->GetStartStep();
//  if (startstep >= 0) { Load(std::to_string((long long)startstep)+"/"); }
  int inputindex = 0;
  if (tstepper_->isOutputOnStep()) {
    inputindex = tstepper_->GetStartStep();
  } else {
    inputindex = tstepper_->GetStartTimeIndex();
  }
  if (inputindex >= 0) { Load(std::to_string((long long)inputindex)+"/"); }

  sstepper_->Gradient(*p[index_],*dpdxyz[index_]);
  uvwnonlinear->GetCubeCplxJoin() = uvw[index_]->GetCubeCplxJoin();
  uvwnonlinear->Dealiasing();
  uvwnonlinear->DftBackward();
  
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatScalarVariable::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FlatScalarVariable > ";
  Variable::Print(stream);
  stream << std::endl;
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
