/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <iomanip>
#include "solver_channel_flat.h"
#include "domain_fixed.h"
#include "stepper_space_chebyshev.h"
#include "stepper_time_rk3.h"
#include "phase_carrier_flat.h"
#include "phase_scalar_flat.h"



#include "p3dfft.h"



namespace yueliangyi {
namespace pseudospectral {
namespace leveltwo {
  
  
  using namespace levelbase;
  using namespace levelone;
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatChannelSolver::FlatChannelSolver(int& argc, char**& argv) : System(argc,argv) {
  PS_DEBUG_TRACE_ENTER("FlatChannelSolver::FlatChannelSolver()")
  
  
  filelog_.SetFileName(IO_FILE_LOG_COMPUTATION);
  filelog_.Configure(printlog_,writelog_);
  
  
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
  
  

  
  
  
  
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelSolver::FlatChannelSolver()")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatChannelSolver::~FlatChannelSolver(void) {
  
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
bool FlatChannelSolver::Configure(void) {
  PS_DEBUG_TRACE_ENTER("FlatChannelSolver::Configure()")
  
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
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelSolver::Configure()")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatChannelSolver::Execute(void) {
  PS_DEBUG_TRACE_ENTER("FlatChannelSolver::Execute(void)")
  
  
  ostringstream& message = filelog_.GetMessage();
  PrintLogHead(message);
  
  
  const dvec& feature = phase_[0].obj->GetFeatureGlobal();

  
  
  
  if (tstepper_.obj->isOutput()) {
//    phase_[0].obj->Save(step);
  }
  
  
  
  
  while (tstepper_.obj->isForward()) {
    
    
    if (tstepper_.obj->isOutputOnStep()) {
    tstepper_.obj->Update(feature(19));
    } else {
      tstepper_.obj->Update2(feature(19));
      phase_[0].obj->Update();
      if (phasenum_ > 1) {
        phase_[1].obj->Update();
      }
    }
    
    
    
    
    const int step = tstepper_.obj->GetStep();
    
    
    // 2. Forward the whole system by one step.
    //    Since the carrier phase is always the first one being proceeded,
    //    information of carrier is updated earlier than the others.
    //    The domain and space stepper is updated during each stage if necessary.
    
    for (int stage = 0; stage < 3; ++stage) {
      
      
      
      phase_[0].obj->Predict(stage);
      phase_[0].obj->Correct(stage);
      
      
      if (phasenum_ > 1) {
        phase_[1].obj->Predict(stage);
        phase_[1].obj->Correct(stage);
      }

      
      
      
      // If there are some sediment phases.
      switch (phasenum_) {
        case 3:
        case 2:
          break;
        default:
          break;
      }
      
      
      
    } // end of for
    
    
    // update alternation flag
    tstepper_.obj->SetAlternationFlag();
    
    
    if (tstepper_.obj->isOutput()) {
      
      int outputindex = 0;
      if (tstepper_.obj->isOutputOnStep()) {
        outputindex = step;
      } else {
        outputindex = tstepper_.obj->GetOutputTimeIndex();
      }
      

//      if (mpi_.isMain()) {
//        std::cout << outputindex << std::endl;
//      }
      phase_[0].obj->Save(std::to_string((long long)outputindex)+"/");
      if (phasenum_ > 1) {
        phase_[1].obj->Save(std::to_string((long long)outputindex)+"/");
      }
      
      
      _message_.str("");
      _message_.clear();
      
      
    }
    
    
    
    PrintLogBody(message);
    filelog_.Save();
    
    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
  }
  
  
  
  
//  filelog_.Save();
  
  
  
  
  
//  const string outputdir = "./output/";
//  fstream file;
//  file.open (outputdir+"rank"+std::to_string(mpi_.GetRank()),fstream::out);
//  file << _message_.str();
//  file.close();
  
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelSolver::Execute(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatChannelSolver::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("FlatChannelSolver::Print()")

//  mpi_.Print(stream); stream << std::endl;
//  dtb_.Print(stream);
  
  
  
  stream << std::setw(30) << std::right << "FlatChannelSolver > ";
  System::Print(stream);
  stream << std::endl;
  
  
  
  
  PS_DEBUG_TRACE_LEAVE("FlatChannelSolver::Print()")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatChannelSolver::PrintLogHead(ostringstream& stream) const {
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
void FlatChannelSolver::PrintLogBody(ostringstream& stream) const {
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
    
    << std::setw(12) << feature(4)
    << std::setw(12) << feature(5)
    << std::setw(12) << feature(6)
    << std::setw(12) << feature(7)
    << std::setw(12) << feature(8)
    
    << std::setw(12) << feature(10)
    << std::setw(12) << feature(11)
    << std::setw(12) << feature(12)
    << std::setw(12) << feature(13)
    << std::setw(12) << feature(14)
    << std::setw(12) << feature(15)
    
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
