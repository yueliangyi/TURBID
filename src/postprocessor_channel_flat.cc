/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <iomanip>
#include "postprocessor_channel_flat.h"
#include "domain_fixed.h"
#include "stepper_space_chebyshev.h"
#include "stepper_time_rk3.h"
#include "phase_carrier_flat.h"
#include "phase_scalar_flat.h"

#include "variable_carrier_flat.h"
#include <sys/stat.h>



#include "p3dfft.h"



namespace yueliangyi {
namespace pseudospectral {
namespace leveltwo {
  
  
  using namespace levelbase;
  using namespace levelone;
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatChannelPostprocessor::FlatChannelPostprocessor(int& argc, char**& argv) : System(argc,argv) {
  PS_DEBUG_TRACE_ENTER("FlatChannelPostprocessor::FlatChannelPostprocessor()")
  
  
  filelog_.SetFileName(IO_FILE_LOG_COMPUTATION);
  filelog_.Configure(printlog_,writelog_);
  
  //
  work_index_range = fileini_.GetVectorInt(SECTION_POST_PROCESSOR,KEY_SPOST_WORK_INDEX_RANGE);
  
  
  //
  typedomain_   = ConvertCase(fileini_.GetValueString(SECTION_DOMAIN,KEY_SD_APPLIED_DOMAIN_TYPE),1);
  typesstepper_ = ConvertCase(fileini_.GetValueString(SECTION_STEPPER_SPACE,KEY_SSS_APPLIED_STEPPER_TYPE),1);
  typetstepper_ = ConvertCase(fileini_.GetValueString(SECTION_STEPPER_TIME,KEY_SST_APPLIED_STEPPER_TYPE),1);
  
  //
  if (typedomain_ == "FIXED") {
    
    
    FixedDomain::cfg_t* cfg = new FixedDomain::cfg_t;
    cfg->Load(fileini_);
    
    domain_.cfg = cfg;
    domain_.obj = new FixedDomain(mpi_,dcp_,domain_.cfg);
    
  }
  
  

  //
  if (typesstepper_ == "CHEBYSHEV") {
    
    ChebyshevSpaceStepper::cfg_t* cfg = new ChebyshevSpaceStepper::cfg_t;
    cfg->Load(fileini_);
    
    sstepper_.cfg = cfg;
    sstepper_.obj = new ChebyshevSpaceStepper(mpi_,dcp_,domain_.obj,sstepper_.cfg);
    
  }
  
  
  //
  if (typetstepper_ == "RUNGEKUTTA3") {
    
    RK3TimeStepper::cfg_t* cfg = new RK3TimeStepper::cfg_t;
    cfg->Load(fileini_);
    
    tstepper_.cfg = cfg;
    tstepper_.obj = new RK3TimeStepper(mpi_,dcp_,tstepper_.cfg);
    
  }
  

  for (int idphase = 0; idphase < phasenum_; ++idphase) {
    
    if (idphase == 0) {
      FlatCarrierPhase::cfg_t* cfg = new FlatCarrierPhase::cfg_t;
      cfg->Load(fileini_);
      phase_t phase;
      phase.cfg = cfg;
      phase.obj = new FlatCarrierPhase(mpi_,
                                       dcp_,
                                       domain_.obj,
                                       sstepper_.obj,
                                       tstepper_.obj,
                                       phase.cfg,
                                       idphase);
      phase_.push_back(phase);
    } else {
      FlatScalarPhase::cfg_t* cfg = new FlatScalarPhase::cfg_t;
      cfg->Load(fileini_);
      phase_t phase;
      phase.cfg = cfg;
      phase.obj = new FlatScalarPhase(mpi_,
                                      dcp_,
                                      domain_.obj,
                                      sstepper_.obj,
                                      tstepper_.obj,
                                      phase.cfg,
                                      idphase);
      phase_.push_back(phase);
    }
  }
  
  

  
  
  
  
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelPostprocessor::FlatChannelPostprocessor()")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatChannelPostprocessor::~FlatChannelPostprocessor(void) {
  
  if (domain_  .obj != NULL) { delete domain_  .obj; domain_  .obj = NULL; }
  if (sstepper_.cfg != NULL) { delete sstepper_.cfg; sstepper_.cfg = NULL; }
  if (tstepper_.cfg != NULL) { delete tstepper_.cfg; tstepper_.cfg = NULL; }
  
  
  for (auto it=phase_.begin(); it!=phase_.end(); ++it) {
    delete it->cfg; it->cfg = NULL;
    delete it->obj; it->obj = NULL;
  }
  
  
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatChannelPostprocessor::Configure(void) {
  PS_DEBUG_TRACE_ENTER("FlatChannelPostprocessor::Configure()")
  
  domain_.obj->Configure();
  
  sstepper_.obj->Configure();
  
  
  sstepper_.obj->Update();
  
  domain_.obj->Update();
  
  tstepper_.obj->Configure();

  phase_[0].obj->Configure();

  if (phasenum_>1) {
    phase_[1].obj->Configure();
  }
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelPostprocessor::Configure()")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatChannelPostprocessor::Execute(void) {
  PS_DEBUG_TRACE_ENTER("FlatChannelPostprocessor::Execute(void)")


  ostringstream& message = filelog_.GetMessage();
  PrintLogHead(message);
  
  levelbase::rccube* votex = new rccube(mpi_,dcp_,1);
  
  
  for (int index = work_index_range(0); index <= work_index_range(1); ++index) {
  
    // load result files
    FlatCarrierVariable* p0var = (FlatCarrierVariable*)phase_[0].obj->GetVariable();
    p0var->Load(std::to_string((long long)index)+"/");

    const int p0index = p0var->GetIndex();
    const int rank = mpi_.GetRank();
    const imat& dimension = dcp_.GetDimension();

    sstepper_.obj->Gradient(*p0var->uvw[p0index],*p0var->dudxyz[p0index],0);
    sstepper_.obj->Gradient(*p0var->uvw[p0index],*p0var->dvdxyz[p0index],1);
    sstepper_.obj->Gradient(*p0var->uvw[p0index],*p0var->dwdxyz[p0index],2);
    p0var->dudxyz[p0index]->DftBackward();
    p0var->dvdxyz[p0index]->DftBackward();
    p0var->dwdxyz[p0index]->DftBackward();
    
    
    for (int indz = 0; indz < dimension(10*2+8,rank); ++indz) {
      for (int indy = 0; indy < dimension(10*2+7,rank); ++indy) {
        for (int indx = 0; indx < dimension(10*2+6,rank); ++indx) {

          dmat vgradtensor(3,3);
          vgradtensor(0,0) = p0var->dudxyz[p0index]->GetCubeReal(0)(indx,indy,indz);
          vgradtensor(0,1) = p0var->dudxyz[p0index]->GetCubeReal(1)(indx,indy,indz);
          vgradtensor(0,2) = p0var->dudxyz[p0index]->GetCubeReal(2)(indx,indy,indz);

          vgradtensor(1,0) = p0var->dvdxyz[p0index]->GetCubeReal(0)(indx,indy,indz);
          vgradtensor(1,1) = p0var->dvdxyz[p0index]->GetCubeReal(1)(indx,indy,indz);
          vgradtensor(1,2) = p0var->dvdxyz[p0index]->GetCubeReal(2)(indx,indy,indz);

          vgradtensor(2,0) = p0var->dwdxyz[p0index]->GetCubeReal(0)(indx,indy,indz);
          vgradtensor(2,1) = p0var->dwdxyz[p0index]->GetCubeReal(1)(indx,indy,indz);
          vgradtensor(2,2) = p0var->dwdxyz[p0index]->GetCubeReal(2)(indx,indy,indz);

          
          // Swirling Strength
          arma::cx_vec eigenvalue;
          arma::cx_mat eigenvector;
          eig_gen(eigenvalue,eigenvector,vgradtensor);
          votex->GetCubeReal(0)(indx,indy,indz) = max(unique(abs(imag(eigenvalue))));
          
          // Q-criterion
//          votex->GetCubeReal(0)(indx,indy,indz) = -arma::accu(vgradtensor%vgradtensor.t());

        }
      }
    }
    
    
    struct stat status = {0};
    const mode_t mode = S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH;
    const string path = "./"+std::to_string((long long)index)+"/";
    if (mpi_.isMain() &&
        stat(path.c_str(),&status)==-1 &&
        mkdir(path.c_str(),mode)==-1) {
      name_ += "Save(int,ioctrl_t,string) > Failed to Create Target Folder!";
      throw std::runtime_error(name_);
    }

    fstream filestream;
    if (mpi_.isMain()) {
      filestream.open(path+"phase_00_vortex.dat",fstream::out|fstream::binary);
    }
    votex->SaveReal(filestream);
    if (mpi_.isMain()) {
      filestream.close();
    }
    
    
//    if (mpi_.isMain()) {
//      filestream.open(path+"phase_00_velocity.dat",fstream::out|fstream::binary);
//    }
//    p0var->uvw[p0index]->DftBackward();
//    p0var->uvw[p0index]->SaveReal(filestream);
//    if (mpi_.isMain()) {
//      filestream.close();
//    }
    
    
    
    
    
    
//    _message_.str("");
//    _message_.clear();
    
    
    PrintLogBody(message);
    filelog_.Save();



    MPI_Barrier(MPI_COMM_WORLD);


  }
  
  delete votex;
  
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelPostprocessor::Execute(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatChannelPostprocessor::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("FlatChannelPostprocessor::Print()")

//  mpi_.Print(stream); stream << std::endl;
//  dtb_.Print(stream);
  
  
  
  stream << std::setw(30) << std::right << "FlatChannelPostprocessor > ";
  System::Print(stream);
  stream << std::endl;
  
  
  
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelPostprocessor::Print()")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatChannelPostprocessor::PrintLogHead(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"PrintLogHead(ostringstream)")
  
  
  if (mpi_.isMain()) {
    stream << IO_FILE_SEPERATOR_STAR << std::endl;
    this->Print(stream);
    mpi_.Print(stream);
    dcp_.Print(stream);
    domain_.obj->Print(stream);
    sstepper_.obj->Print(stream);
    tstepper_.obj->Print(stream);
    phase_[0].obj->Print(stream);
    if (phasenum_>1) {
      phase_[1].obj->Print(stream);
    }
    stream << IO_FILE_SEPERATOR_STAR << std::endl;
  }
  
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"PrintLogHead(ostringstream)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatChannelPostprocessor::PrintLogBody(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"PrintLogBody(ostringstream)")
  
  if (mpi_.isMain()) {
    
    const dvec& feature = phase_[0].obj->GetFeatureGlobal();
    stream
    << std::setw(12) << tstepper_.obj->GetStep()
    << std::setw(16) << std::setprecision(10) << tstepper_.obj->GetTime()
    << std::setw(12) << std::setprecision(6)  << tstepper_.obj->GetCfl()
    << std::setw(12) << feature(0)
    << std::setw(12) << feature(1)
    << std::setw(12) << feature(2)
    << std::setw(12) << feature(3)
    << std::setw(12) << feature(6)
    << std::setw(12) << feature(9)
    << std::setw(12) << feature(10)
    << std::setw(12) << feature(11)
//    << std::setw(12) << feature(12)
//    << std::setw(12) << feature(14)
//    << std::setw(12) << feature(15)
    << std::setw(12) << feature(21)
    << std::setw(12) << feature(22)
    << std::setw(12) << feature(23)
    << std::setw(12) << feature(24)
    << std::setw(12) << feature(25)
    << std::setw(12) << feature(28)
    << std::setw(12) << feature(29)
    << std::endl;
    
  }
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"PrintLogBody(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace leveltwo
} // end namespace pseudospectral
} // end namespace yueliangyi
