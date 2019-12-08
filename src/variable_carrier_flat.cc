/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <sys/stat.h>
#include <iomanip>
#include "variable_carrier_flat.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatCarrierVariable::FlatCarrierVariable(const mpi_t& mpi,
                                         const dcp_t& dcp,
                                         Domain      * domain,
                                         SpaceStepper* sstepper,
                                         TimeStepper * tstepper,
                                         const int index) :
  Variable(mpi,dcp,domain,sstepper,tstepper,index) {
  PS_DEBUG_TRACE_ENTER(name_+"FlatCarrierVariable(mpi_t,dcp_t)")
    
//  dealiasingu = new rccube(mpi_,dcp_);
//  dealiasingv = new rccube(mpi_,dcp_);
//  dealiasingw = new rccube(mpi_,dcp_);
//  
//  
//  advu = new rccube(mpi_,dcp_);
//  advv = new rccube(mpi_,dcp_);
//  advw = new rccube(mpi_,dcp_);
//  
//  hu   = new rccube(mpi_,dcp_);
//  hv   = new rccube(mpi_,dcp_);
//  hw   = new rccube(mpi_,dcp_);
//  
//  us   = new rccube(mpi_,dcp_);
//  vs   = new rccube(mpi_,dcp_);
//  ws   = new rccube(mpi_,dcp_);
//
//  lapu = new rccube(mpi_,dcp_);
//  lapv = new rccube(mpi_,dcp_);
//  lapw = new rccube(mpi_,dcp_);
//  
//  rhsu = new rccube(mpi_,dcp_);
//  rhsv = new rccube(mpi_,dcp_);
//  rhsw = new rccube(mpi_,dcp_);
//  rhsp = new rccube(mpi_,dcp_);
//
//  divs = new rccube(mpi_,dcp_);
    
    
    
    uvwnonlinear = new rccube(mpi_,dcp_,3);
    advuvw = new rccube(mpi_,dcp_,3);
    huvw = new rccube(mpi_,dcp_,3,false);
    suvw = new rccube(mpi_,dcp_,3,false);
    lapuvw = new rccube(mpi_,dcp_,3,false);
    rhsuvw = new rccube(mpi_,dcp_,3,false);
    rhsp = new rccube(mpi_,dcp_,1,false);
    divs = new rccube(mpi_,dcp_,1,false);
    
    
    tmpu1 = new rccube(mpi_,dcp_,1,false);
    tmpv1 = new rccube(mpi_,dcp_,1,false);
    tmpw1 = new rccube(mpi_,dcp_,1,false);
    tmpp1 = new rccube(mpi_,dcp_,1,false);
    tmpu2 = new rccube(mpi_,dcp_,1,false);
    tmpv2 = new rccube(mpi_,dcp_,1,false);
    tmpw2 = new rccube(mpi_,dcp_,1,false);
    tmpp2 = new rccube(mpi_,dcp_,1,false);
  
  PS_DEBUG_TRACE_LEAVE(name_+"FlatCarrierVariable(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatCarrierVariable::~FlatCarrierVariable() {
  
//  if (dealiasingu != NULL) { delete dealiasingu; dealiasingu = NULL; }
//  if (dealiasingv != NULL) { delete dealiasingv; dealiasingv = NULL; }
//  if (dealiasingw != NULL) { delete dealiasingw; dealiasingw = NULL; }
//  
//  if (advu != NULL) { delete advu; advu = NULL; }
//  if (advv != NULL) { delete advv; advv = NULL; }
//  if (advw != NULL) { delete advw; advw = NULL; }
//  
//  if (hu   != NULL) { delete hu  ; hu   = NULL; }
//  if (hv   != NULL) { delete hv  ; hv   = NULL; }
//  if (hw   != NULL) { delete hw  ; hw   = NULL; }
//  
//  if (us   != NULL) { delete us  ; us   = NULL; }
//  if (vs   != NULL) { delete vs  ; vs   = NULL; }
//  if (ws   != NULL) { delete ws  ; ws   = NULL; }
//  
//  if (lapu != NULL) { delete lapu; lapu = NULL; }
//  if (lapv != NULL) { delete lapv; lapv = NULL; }
//  if (lapw != NULL) { delete lapw; lapw = NULL; }
//  
//  if (rhsu != NULL) { delete rhsu; rhsu = NULL; }
//  if (rhsv != NULL) { delete rhsv; rhsv = NULL; }
//  if (rhsw != NULL) { delete rhsw; rhsw = NULL; }
//  if (rhsp != NULL) { delete rhsp; rhsp = NULL; }
//  
//  if (divs != NULL) { delete divs; divs = NULL; }
  
  
  if (uvwnonlinear != NULL) { delete uvwnonlinear; uvwnonlinear = NULL; }
  if (advuvw != NULL) { delete advuvw; advuvw = NULL; }
  if (huvw != NULL) { delete huvw; huvw = NULL; }
  if (suvw != NULL) { delete suvw; suvw = NULL; }
  if (lapuvw != NULL) { delete lapuvw; lapuvw = NULL; }
  if (rhsuvw != NULL) { delete rhsuvw; rhsuvw = NULL; }
  if (rhsp != NULL) { delete rhsp; rhsp = NULL; }
  if (divs != NULL) { delete divs; divs = NULL; }
  
  if (tmpu1 != NULL) { delete tmpu1; tmpu1 = NULL; }
  if (tmpv1 != NULL) { delete tmpv1; tmpv1 = NULL; }
  if (tmpw1 != NULL) { delete tmpw1; tmpw1 = NULL; }
  if (tmpp1 != NULL) { delete tmpp1; tmpp1 = NULL; }
  
  if (tmpu2 != NULL) { delete tmpu2; tmpu2 = NULL; }
  if (tmpv2 != NULL) { delete tmpv2; tmpv2 = NULL; }
  if (tmpw2 != NULL) { delete tmpw2; tmpw2 = NULL; }
  if (tmpp2 != NULL) { delete tmpp2; tmpp2 = NULL; }
  
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierVariable::Save(const string subfolder, const ioctrl_t level) {
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
//      advu->DftBackward(); advu->SaveReal(filestream_);
//      advv->DftBackward(); advv->SaveReal(filestream_);
//      advw->DftBackward(); advw->SaveReal(filestream_);
//      us  ->DftBackward(); us  ->SaveReal(filestream_);
//      vs  ->DftBackward(); vs  ->SaveReal(filestream_);
//      ws  ->DftBackward(); ws  ->SaveReal(filestream_);
//      lapu->DftBackward(); lapu->SaveReal(filestream_);
//      lapv->DftBackward(); lapv->SaveReal(filestream_);
//      lapw->DftBackward(); lapw->SaveReal(filestream_);
//      rhsu->DftBackward(); rhsu->SaveReal(filestream_);
//      rhsv->DftBackward(); rhsv->SaveReal(filestream_);
//      rhsw->DftBackward(); rhsw->SaveReal(filestream_);
//      rhsp->DftBackward(); rhsp->SaveReal(filestream_);
//      divs->DftBackward(); divs->SaveReal(filestream_);
      
//      advuvw->DftBackward(); advuvw->SaveReal(filestream_);
//      suvw  ->DftBackward(); suvw  ->SaveReal(filestream_);
//      lapuvw->DftBackward(); lapuvw->SaveReal(filestream_);
//      rhsuvw->DftBackward(); rhsuvw->SaveReal(filestream_);
//      rhsp  ->DftBackward(); rhsp  ->SaveReal(filestream_);
//      divs  ->DftBackward(); divs  ->SaveReal(filestream_);
      
    case IOControlLevelTwo:
      filestream_ << "IOControlLevelTwo" << std::endl;
      // Save partial derivatives.
//      dudx[index_]->DftBackward(); dudx[index_]->SaveReal(filestream_);
//      dudy[index_]->DftBackward(); dudy[index_]->SaveReal(filestream_);
//      dudz[index_]->DftBackward(); dudz[index_]->SaveReal(filestream_);
//      dvdx[index_]->DftBackward(); dvdx[index_]->SaveReal(filestream_);
//      dvdy[index_]->DftBackward(); dvdy[index_]->SaveReal(filestream_);
//      dvdz[index_]->DftBackward(); dvdz[index_]->SaveReal(filestream_);
//      dwdx[index_]->DftBackward(); dwdx[index_]->SaveReal(filestream_);
//      dwdy[index_]->DftBackward(); dwdy[index_]->SaveReal(filestream_);
//      dwdz[index_]->DftBackward(); dwdz[index_]->SaveReal(filestream_);
//      dpdx[index_]->DftBackward(); dpdx[index_]->SaveReal(filestream_);
//      dpdy[index_]->DftBackward(); dpdy[index_]->SaveReal(filestream_);
//      dpdz[index_]->DftBackward(); dpdz[index_]->SaveReal(filestream_);
      
//      dudxyz[index_]->DftBackward(); dudxyz[index_]->SaveReal(filestream_);
//      dvdxyz[index_]->DftBackward(); dvdxyz[index_]->SaveReal(filestream_);
//      dwdxyz[index_]->DftBackward(); dwdxyz[index_]->SaveReal(filestream_);
//      dpdxyz[index_]->DftBackward(); dpdxyz[index_]->SaveReal(filestream_);
      
    case IOControlLevelOne:
      filestream_ << "IOControlLevelOne" << std::endl;
      // Save basic variables.
      // Save hu/y/z is necessary.
//      u[index_]->DftBackward(); u[index_]->SaveReal(filestream_);
//      v[index_]->DftBackward(); v[index_]->SaveReal(filestream_);
//      w[index_]->DftBackward(); w[index_]->SaveReal(filestream_);
//      p[index_]->DftBackward(); p[index_]->SaveReal(filestream_);
//      hu       ->DftBackward(); hu->       SaveReal(filestream_);
//      hv       ->DftBackward(); hv->       SaveReal(filestream_);
//      hw       ->DftBackward(); hw->       SaveReal(filestream_);
//      uvw [index_]->DftBackward(); uvw [index_]->SaveReal(filestream_);
//      p   [index_]->DftBackward(); p   [index_]->SaveReal(filestream_);
//      huvw        ->DftBackward(); huvw        ->SaveReal(filestream_);
      
      uvw[index_]->SaveCplx(filestream_);
      p  [index_]->SaveCplx(filestream_);
      
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
bool FlatCarrierVariable::Load(const string subfolder, const ioctrl_t level) {
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
      
//      advu->LoadReal(filestream_); advu->DftForward();
//      advv->LoadReal(filestream_); advv->DftForward();
//      advw->LoadReal(filestream_); advw->DftForward();
//      us  ->LoadReal(filestream_); us  ->DftForward();
//      vs  ->LoadReal(filestream_); vs  ->DftForward();
//      ws  ->LoadReal(filestream_); ws  ->DftForward();
//      lapu->LoadReal(filestream_); lapu->DftForward();
//      lapv->LoadReal(filestream_); lapv->DftForward();
//      lapw->LoadReal(filestream_); lapw->DftForward();
//      rhsu->LoadReal(filestream_); rhsu->DftForward();
//      rhsv->LoadReal(filestream_); rhsv->DftForward();
//      rhsw->LoadReal(filestream_); rhsw->DftForward();
//      rhsp->LoadReal(filestream_); rhsp->DftForward();
//      divs->LoadReal(filestream_); divs->DftForward();
      
//      advuvw->LoadReal(filestream_); advuvw->DftForward();
//      suvw  ->LoadReal(filestream_); suvw  ->DftForward();
//      lapuvw->LoadReal(filestream_); lapuvw->DftForward();
//      rhsuvw->LoadReal(filestream_); rhsuvw->DftForward();
//      rhsp  ->LoadReal(filestream_); rhsp  ->DftForward();
//      divs  ->LoadReal(filestream_); divs  ->DftForward();
      
    case IOControlLevelTwo:

      // Save partial derivatives.
      if (mpi_.isMain() && levelname!="IOControlLevelTwo") {
        name_ += "Load(int,ioctrl_t,string) > Wrong Request of IO Level Two!";
        throw std::runtime_error(name_);
      }
      
//      dudx[index_]->LoadReal(filestream_); dudx[index_]->DftForward();
//      dudy[index_]->LoadReal(filestream_); dudy[index_]->DftForward();
//      dudz[index_]->LoadReal(filestream_); dudz[index_]->DftForward();
//      dvdx[index_]->LoadReal(filestream_); dvdx[index_]->DftForward();
//      dvdy[index_]->LoadReal(filestream_); dvdy[index_]->DftForward();
//      dvdz[index_]->LoadReal(filestream_); dvdz[index_]->DftForward();
//      dwdx[index_]->LoadReal(filestream_); dwdx[index_]->DftForward();
//      dwdy[index_]->LoadReal(filestream_); dwdy[index_]->DftForward();
//      dwdz[index_]->LoadReal(filestream_); dwdz[index_]->DftForward();
//      dpdx[index_]->LoadReal(filestream_); dpdx[index_]->DftForward();
//      dpdy[index_]->LoadReal(filestream_); dpdy[index_]->DftForward();
//      dpdz[index_]->LoadReal(filestream_); dpdz[index_]->DftForward();
      
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
      
//      u[index_]->LoadReal(filestream_); u[index_]->DftForward();
//      v[index_]->LoadReal(filestream_); v[index_]->DftForward();
//      w[index_]->LoadReal(filestream_); w[index_]->DftForward();
//      p[index_]->LoadReal(filestream_); p[index_]->DftForward();
//      hu->       LoadReal(filestream_); hu       ->DftForward();
//      hv->       LoadReal(filestream_); hv       ->DftForward();
//      hw->       LoadReal(filestream_); hw       ->DftForward();
//      uvw[index_]->LoadReal(filestream_); uvw[index_]->DftForward();
//      p  [index_]->LoadReal(filestream_); p  [index_]->DftForward();
//      huvw       ->LoadReal(filestream_); huvw       ->DftForward();
      
      uvw[index_]->LoadCplx(filestream_);
      p  [index_]->LoadCplx(filestream_);
      
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
bool FlatCarrierVariable::Configure(void) {
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
  
  
  
//  // 2.1 Get the gradient of pressure.
//  sstepper_->Gradient(*p   [index_],
//                      *dpdx[index_],
//                      *dpdy[index_],
//                      *dpdz[index_]);
//  
//  
//  
//  dealiasingu->GetCubeCplxJoin() = u[index_]->GetCubeCplxJoin();
//  dealiasingv->GetCubeCplxJoin() = v[index_]->GetCubeCplxJoin();
//  dealiasingw->GetCubeCplxJoin() = w[index_]->GetCubeCplxJoin();
//  
//  
//  dealiasingu->Dealiasing();
//  dealiasingv->Dealiasing();
//  dealiasingw->Dealiasing();
//  
//  dealiasingu->DftBackward();
//  dealiasingv->DftBackward();
//  dealiasingw->DftBackward();
  

  sstepper_->Gradient(*p[index_],*dpdxyz[index_]);
  uvwnonlinear->GetCubeCplxJoin() = uvw[index_]->GetCubeCplxJoin();
  uvwnonlinear->DftBackward();
  
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatCarrierVariable::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FlatCarrierVariable > ";
  Variable::Print(stream);
  stream << std::endl;
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
