/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "phase_carrier_flat.h"
#include <iomanip>
#include <stdexcept>
#include "variable_carrier_flat.h"
#include "boundary_carrier_flat.h"
#include <sys/stat.h>
#include <cstdlib>
#include <math.h>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Configuration::Load(IniFile& fileini) {
  PS_DEBUG_TRACE_ENTER("FlatCarrierPhase::Configuration::Load(IniFile)")
  
  Phase::Configuration::Load(fileini);
  
  reynolds = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_REYNOLDS_NUMBER);
  froude = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_FROUDE_NUMBER);
  averagedconcentration = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_AVERAGED_CONCENTRATION);
  specific_gravity = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_SPECIFIC_GRAVITY);
//  rossby = fileini.GetVectorDouble(SECTION_PHASE_CARRIER,KEY_SPC_ROSSBY_NUMBER);
  affected_by_scalar_phase = fileini.GetValueBool(SECTION_PHASE_CARRIER,KEY_SPC_AFFECTED_BY_SCALAR_PHASE);

  zboundarytypeu1 = fileini.GetVectorDouble(SECTION_PHASE_CARRIER,KEY_SPC_VERTICAL_TYPE_U1);
  zboundarytypeu2 = fileini.GetVectorDouble(SECTION_PHASE_CARRIER,KEY_SPC_VERTICAL_TYPE_U2);
  zboundarytypeu3 = fileini.GetVectorDouble(SECTION_PHASE_CARRIER,KEY_SPC_VERTICAL_TYPE_U3);
  zboundarytypep  = fileini.GetVectorDouble(SECTION_PHASE_CARRIER,KEY_SPC_VERTICAL_TYPE_P);

  string force = ConvertCase(fileini.GetValueString(SECTION_PHASE_CARRIER,KEY_SPC_APPLIED_FORCE_TYPE),1);
  
  if (force == "NONE") { forcetype = kForceNone; }
  else if (force == "CONSTANT") { forcetype = kForceConstant; }
  else if (force == "OSCILLATORY") { forcetype = kForceOscillatory; }
  else if (force == "PULSATILE") { forcetype = kForcePulsatile; }
  else if (force == "DAMPING") { forcetype = kForceDamping; }
  else {
    throw std::runtime_error(string("FlatCarrierPhase::Configuration::")+
                             "Load(IniFile) > Unsupported Force Type!");
  }
  forcedirection  = fileini.GetVectorInt(SECTION_PHASE_CARRIER,KEY_SPC_FORCE_DIRECTION);
  forceamplitude  = fileini.GetVectorDouble(SECTION_PHASE_CARRIER,KEY_SPC_FORCE_AMPLITUDE);

  
  PS_DEBUG_TRACE_LEAVE("FlatCarrierPhase::Configuration::Load(IniFile)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatCarrierPhase::FlatCarrierPhase(const mpi_t& mpi,
                                   const dcp_t& dcp,
                                   Domain      * domain,
                                   SpaceStepper* sstepper,
                                   TimeStepper * tstepper,
                                   const Base::cfg_t* cfg,
                                   const int index) :
  Phase(mpi,dcp,domain,sstepper,tstepper,cfg,index) {
  PS_DEBUG_TRACE_ENTER(name_+"FlatCarrierPhase(mpi_t,dcp_t,cfg_t)")
    
    
    variable_ = new FlatCarrierVariable(mpi,dcp,domain,sstepper,tstepper,index_);
    boundary_ = new FlatCarrierBoundary(mpi,dcp,domain,sstepper,tstepper,variable_,index_);
    
    
//    velocitylhs_ = new rdcube(mpi_,dcp_,3);
//    velocityrhs_ = new rdcube(mpi_,dcp_,3);
    
    
  PS_DEBUG_TRACE_LEAVE(name_+"FlatCarrierPhase(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatCarrierPhase::~FlatCarrierPhase(void) {
  if (variable_ != NULL) { delete variable_; variable_ = NULL; }
  if (boundary_ != NULL) { delete boundary_; boundary_ = NULL; }
//  if (velocitylhs_ != NULL) { delete velocitylhs_; velocitylhs_ = NULL; }
//  if (velocityrhs_ != NULL) { delete velocityrhs_; velocityrhs_ = NULL; }
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Save(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Save(string,ioctrl_t)")
  
  
  
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
  
  variable_->Save(subfolder);
  boundary_->Save(subfolder);
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Save(string,ioctrl_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Load(const string subfolder, const ioctrl_t level) {
  PS_DEBUG_TRACE_ENTER(name_+"Load(string,ioctrl_t)")
  
  
  if (mpi_.isMain()) {
    filestream_.open(directory_+subfolder+filename_,fstream::in|fstream::binary);
  }
  
  
  if (mpi_.isMain()) {
    filestream_.close();
  }
  
  variable_->Load(subfolder);
  boundary_->Load(subfolder);
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Load(string,ioctrl_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  cfg_t* pcfg = (cfg_t*) cfg_;
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  
  
  // Construct the stage coefficients for Helmholtz solution.
  const double dt  = tstepper_->GetDt();
  dvec stepcoefficient(3);
  stepcoefficient(0) = -pcfg->reynolds/dt/tstepper_->GetCoefficients(2,0);
  stepcoefficient(1) = -pcfg->reynolds/dt/tstepper_->GetCoefficients(2,1);
  stepcoefficient(2) = -pcfg->reynolds/dt/tstepper_->GetCoefficients(2,2);
//  dz2dcp_.Construct(sstepper_,pcfg->zboundarytype,stepcoefficient,3,true);
  dz2dcpu1_.Construct(sstepper_,pcfg->zboundarytypeu1,stepcoefficient,1,true);
  dz2dcpu2_.Construct(sstepper_,pcfg->zboundarytypeu2,stepcoefficient,1,true);
  dz2dcpu3_.Construct(sstepper_,pcfg->zboundarytypeu3,stepcoefficient,1,true);

  
  stepcoefficient.fill(0);
//  dz2dcppressure_.Construct(sstepper_,"0 1 0 1",stepcoefficient,1,true);
  dz2dcppressure_.Construct(sstepper_,pcfg->zboundarytypep,stepcoefficient,1,true);

  
  variable_->Configure();
  boundary_->Configure();
  boundary_->SetFlagOpenTop(pcfg->zboundarytypeu1);
  
  
  
  
  
//  const ivec& size = dcp_.GetSize();
//  const size_t itop = size(2)-1;
//  const size_t ibottom = 0;
  const int index = variable_->GetIndex();
  
  
  
//  if (mpi_.isMain()) {
////    const double scale = 0.16*20;
//    const double scale = 0.025*1;
//    const dvec& gridz = domain_->GetGridZ();
//    for (int id = 0; id < size(2); ++id) {
//      variable_->uvw[index]->GetCubeCplx(0)(id,2,2) = cos((4.5)*arma::datum::pi*gridz(id))*scale;
//      variable_->uvw[index]->GetCubeCplx(0)(id,4,3) = cos((7.5)*arma::datum::pi*gridz(id))*scale;
//      variable_->uvw[index]->GetCubeCplx(1)(id,2,2) = cos((7.5)*arma::datum::pi*gridz(id))*scale;
//      variable_->uvw[index]->GetCubeCplx(1)(id,4,3) = cos((9.5)*arma::datum::pi*gridz(id))*scale;
//      variable_->uvw[index]->GetCubeCplx(2)(id,2,2) = cos((6.5)*arma::datum::pi*gridz(id))*scale;
//      variable_->uvw[index]->GetCubeCplx(2)(id,4,3) = cos((5.5)*arma::datum::pi*gridz(id))*scale;
//    }
//  }
  
  
  
  pvariable->uvwnonlinear->GetCubeCplxJoin() = pvariable->uvw[index]->GetCubeCplxJoin();
  pvariable->uvwnonlinear->DftBackward();
  
//  pvariable->uvw[index]->GetCMatCplxJoin().row(itop) = pvariable->uvw[index]->GetCMatCplxJoin().row(itop-1);
  
//  pvariable->uvw[index]->GetCMatCplxJoin().row(itop).fill(cx_double(0,0));
//  pvariable->uvw[index]->GetCMatCplxJoin().row(ibottom).fill(cx_double(0,0));
  
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
bool FlatCarrierPhase::Update(void) {
  
  cfg_t* pcfg = (cfg_t*) cfg_;
  const double dt  = tstepper_->GetDt();
  dvec stepcoefficient(3);
  stepcoefficient(0) = -pcfg->reynolds/dt/tstepper_->GetCoefficients(2,0);
  stepcoefficient(1) = -pcfg->reynolds/dt/tstepper_->GetCoefficients(2,1);
  stepcoefficient(2) = -pcfg->reynolds/dt/tstepper_->GetCoefficients(2,2);
//  dz2dcp_.Assemble(sstepper_,stepcoefficient);
  dz2dcpu1_.Assemble(sstepper_,stepcoefficient);
  dz2dcpu2_.Assemble(sstepper_,stepcoefficient);
  dz2dcpu3_.Assemble(sstepper_,stepcoefficient);
  // no need for pressure

  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Predict(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Predict(int)")
  
  iscorrectionstep_ = false;
  
  //
  cfg_t* pcfg = (cfg_t*) cfg_;
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  const int index = pvariable->GetIndex();
//  const int step  = tstepper_->GetStep();
  const double dt = tstepper_->GetDt();
  const double c1 = tstepper_->GetCoefficients(0,stage);
  const double c2 = tstepper_->GetCoefficients(1,stage);
  const double c3 = tstepper_->GetCoefficients(2,stage);
  Domain::cfg_t* pcfgdomain = (Domain::cfg_t*) domain_->GetCfg();
  
  
  // advection
  NonlinearVelocity(stage);
  
  
  // update H
  pvariable->huvw->GetCubeCplxJoin() =
    c1*pvariable->huvw->GetCubeCplxJoin()-
    dt*pvariable->advuvw->GetCubeCplxJoin();
  
  
  
  // add Coriolis force
//  pvariable->huvw->GetCubeCplx(0) =
//  pvariable->huvw->GetCubeCplx(0)+
//  dt/pcfg->rossby(0)*pvariable->uvw[index]->GetCubeCplx(1)+
//  dt/pcfg->rossby(1)*pvariable->uvw[index]->GetCubeCplx(2);
//  pvariable->huvw->GetCubeCplx(1) =
//  pvariable->huvw->GetCubeCplx(1)-
//  dt/pcfg->rossby(0)*pvariable->uvw[index]->GetCubeCplx(0);
//  pvariable->huvw->GetCubeCplx(2) =
//  pvariable->huvw->GetCubeCplx(2)+
//  dt/pcfg->rossby(1)*pvariable->uvw[index]->GetCubeCplx(0);

//  pvariable->huvw->GetCubeCplx(0) =
//  pvariable->huvw->GetCubeCplx(0)-
//  dt/pcfg->rossby(0)*pvariable->uvw[index]->GetCubeCplx(1);
//  pvariable->huvw->GetCubeCplx(1) =
//  pvariable->huvw->GetCubeCplx(1)+
//  dt/pcfg->rossby(0)*pvariable->uvw[index]->GetCubeCplx(1)+
//  dt/pcfg->rossby(1)*pvariable->uvw[index]->GetCubeCplx(2);
//  pvariable->huvw->GetCubeCplx(2) =
//  pvariable->huvw->GetCubeCplx(2)-
//  dt/pcfg->rossby(1)*pvariable->uvw[index]->GetCubeCplx(1);
  
  
  
  if (pvariable->uvw.size()>1 && pcfg->affected_by_scalar_phase) {
    
    pvariable->ptmp[index]->GetCubeCplxJoin() = pvariable->p[1]->GetCubeCplxJoin();
//    if (dcp_.isMeanZ()) {
//      pvariable->ptmp[index]->GetCMatCplx(0).col(0).fill(cx_double(0,0));
//    }
    
    
    
    switch (pcfgdomain->slope_direction) {
      case 1: // x1 direction
        pvariable->huvw->GetCubeCplx(0) =
        pvariable->huvw->GetCubeCplx(0)+
        dt*pcfg->averagedconcentration/pow(pcfg->froude,2)*(pcfg->specific_gravity-1)*pvariable->ptmp[index]->GetCubeCplx(0)*sin(pcfgdomain->slope);
        break;
        
      case 2: // x2 direction
        pvariable->huvw->GetCubeCplx(1) =
        pvariable->huvw->GetCubeCplx(1)+
        dt*pcfg->averagedconcentration/pow(pcfg->froude,2)*(pcfg->specific_gravity-1)*pvariable->ptmp[index]->GetCubeCplx(0)*sin(pcfgdomain->slope);
        break;
        
      default:
        break;
    }
    pvariable->huvw->GetCubeCplx(2) =
    pvariable->huvw->GetCubeCplx(2)-
    dt*pcfg->averagedconcentration/pow(pcfg->froude,2)*(pcfg->specific_gravity-1)*pvariable->ptmp[index]->GetCubeCplx(0)*cos(pcfgdomain->slope);
    
    

  }
  
  
  // Add the driven force term to velocity rccube
  if (dcp_.isMeanZ()) {
    
    double force1 = 0;
    double force2 = 0;
    
    switch (pcfg->forcetype) {
        
        // no force
      case kForceNone:
        force1 = 0;
        force2 = 0;
        break;
        
        // constant
      case kForceConstant:
        switch (pcfg->forcedirection(0)) {
          case 1:
            force1 = pcfg->forceamplitude(0);
            break;
          case 2:
            force2 = pcfg->forceamplitude(0);
            break;
          default:
            break;
        }
        break;
        
        // for waves
      case kForceOscillatory:
      {
        double force = pcfg->forceamplitude(1);
        switch (pcfg->forcedirection(1)) {
          case 1:
            if (tstepper_->isOutputOnStep()) {
              force1 = force*cos(force*(tstepper_->GetTime()-dt/2));
            } else {
              force1 = force*cos(force*(tstepper_->GetTime()+dt/2));
            }
            break;
          case 2:
            if (tstepper_->isOutputOnStep()) {
              force2 = force*cos(force*(tstepper_->GetTime()-dt/2));
            } else {
              force2 = force*cos(force*(tstepper_->GetTime()+dt/2));
            }
            break;
          default:
            break;
        }
      }
        break;
        
        // wave and current
      case kForcePulsatile:
      {
        switch (pcfg->forcedirection(0)) {
          case 1:
            force1 = pcfg->forceamplitude(0);
            break;
          case 2:
            force2 = pcfg->forceamplitude(0);
            break;
          default:
            break;
        }
        double force = pcfg->forceamplitude(1);
        switch (pcfg->forcedirection(1)) {
          case 1:
            if (tstepper_->isOutputOnStep()) {
              force1 = force1+force*cos(force*(tstepper_->GetTime()-dt/2));
            } else {
              force1 = force1+force*cos(force*(tstepper_->GetTime()+dt/2));
            }
            break;
          case 2:
            if (tstepper_->isOutputOnStep()) {
              force2 = force2+force*cos(force*(tstepper_->GetTime()-dt/2));
            } else {
              force2 = force2+force*cos(force*(tstepper_->GetTime()+dt/2));
            }
            break;
          default:
            break;
        }
      }
        break;
        
        // wave damping
      case kForceDamping:
      {
        switch (pcfg->forcedirection(0)) {
          case 1:
            force1 = pcfg->forceamplitude(0);
            break;
          case 2:
            force2 = pcfg->forceamplitude(0);
            break;
          default:
            break;
        }
        double force = pcfg->forceamplitude(1);
        const double T = pcfg->forceamplitude(0);
        const double t0 = pcfg->forceamplitude(1);
        const double lastperiod = pcfg->forceamplitude(2);
        const double tarpercent = pcfg->forceamplitude(3);
        const double dampparameter = -log(tarpercent)/(lastperiod+0.25)/T;
        
        switch (pcfg->forcedirection(1)) {
          case 1:
            if (tstepper_->isOutputOnStep()) {
              force1 = force1+(force*cos(force*(tstepper_->GetTime()-dt/2))-
                               dampparameter*sin(force*(tstepper_->GetTime()-dt/2)))*
              exp(-dampparameter*(tstepper_->GetTime()-t0));
            } else {
              force1 = force1+(force*cos(force*(tstepper_->GetTime()+dt/2))-
                               dampparameter*sin(force*(tstepper_->GetTime()+dt/2)))*
              exp(-dampparameter*(tstepper_->GetTime()-t0));
            }
            break;
          case 2:
            if (tstepper_->isOutputOnStep()) {
              force2 = force2+(force*cos(force*(tstepper_->GetTime()-dt/2))-
                               dampparameter*sin(force*(tstepper_->GetTime()-dt/2)))*
              exp(-dampparameter*(tstepper_->GetTime()-t0));
            } else {
              force2 = force2+(force*cos(force*(tstepper_->GetTime()+dt/2))-
                               dampparameter*sin(force*(tstepper_->GetTime()+dt/2)))*
              exp(-dampparameter*(tstepper_->GetTime()-t0));
            }
            break;
            break;
          default:
            break;
        }
      }
        break;
        
      default:
        name_ += "Predict() > Not Supported Force Type!";
        throw std::logic_error(name_);
        break;
        
    }
    pvariable->huvw->GetCMatCplx(0).col(0) += dt*force1;
    pvariable->huvw->GetCMatCplx(1).col(0) += dt*force2;

    
    
    const dvec& gridz = domain_->GetGridZ();
    dmat bulkvelocity = arma::trapz(gridz,real(pvariable->uvw[index]->GetCMatCplx(0).col(0)))/2;
    bulkvelocity_ = bulkvelocity(0);
    
    
    
  }
  
  
  pvariable->lapuvw->GetCMatCplx(0) =
    1*sstepper_->GetCmDz2()              *pvariable->uvw[index]->GetCMatCplx(0)+
    1*sstepper_->GetCmDx2().GetCMatCplx()%pvariable->uvw[index]->GetCMatCplx(0)+
    1*sstepper_->GetCmDy2().GetCMatCplx()%pvariable->uvw[index]->GetCMatCplx(0);
  pvariable->lapuvw->GetCMatCplx(1) =
    1*sstepper_->GetCmDz2()              *pvariable->uvw[index]->GetCMatCplx(1)+
    1*sstepper_->GetCmDx2().GetCMatCplx()%pvariable->uvw[index]->GetCMatCplx(1)+
    1*sstepper_->GetCmDy2().GetCMatCplx()%pvariable->uvw[index]->GetCMatCplx(1);
  pvariable->lapuvw->GetCMatCplx(2) =
    1*sstepper_->GetCmDz2()              *pvariable->uvw[index]->GetCMatCplx(2)+
    1*sstepper_->GetCmDx2().GetCMatCplx()%pvariable->uvw[index]->GetCMatCplx(2)+
    1*sstepper_->GetCmDy2().GetCMatCplx()%pvariable->uvw[index]->GetCMatCplx(2);
  

  // 4.1 Compute the term vec{rhs}^{m0}.
  const double redtc3 = -pcfg->reynolds/dt/c3;
  pvariable->rhsuvw->GetCubeCplxJoin() =
    redtc3*c2*pvariable->huvw->GetCubeCplxJoin()+
    redtc3*pvariable->uvw[index]->GetCubeCplxJoin()-
    pvariable->lapuvw->GetCubeCplxJoin();
  
  
  
  // 5.1 Boundary condition.
  boundary_->Apply(stage);
  if (dcp_.isMeanZ()) {
    const size_t itop = GetDcp().GetSize(2)-1;
    pvariable->rhsuvw->GetCMatCplx(0).at(itop,0) = pcfg->zboundarytypeu1(4);
    pvariable->rhsuvw->GetCMatCplx(1).at(itop,0) = pcfg->zboundarytypeu2(4);
    pvariable->rhsuvw->GetCMatCplx(2).at(itop,0) = pcfg->zboundarytypeu3(4);
  }
  
  
  // ---------------------------------------------------------------------------
  // 6. Solve the velocity Helmholtz equation for the intermediate velocity
  // ---------------------------------------------------------------------------
  
  Helmholtz(*pvariable->rhsuvw,*pvariable->suvw,stage);
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Predict(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Correct(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Correct(int)")
  
  iscorrectionstep_ = true;
  
  
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  const int index = pvariable->GetIndex();
  const double dt = tstepper_->GetDt();
  const double c3 = tstepper_->GetCoefficients(2,stage);
  
  const size_t itop = GetDcp().GetSize(2)-1;
  const size_t ibottom = 0;
  
  
  //
  NonlinearPressure(stage);
  
  
  // 1.1 Compute the divergence of intermediate velocity.
  sstepper_->Divergence(*pvariable->suvw,*pvariable->divs);
  
  // 1.3 Compute the right hand side for pressure.
  pvariable->rhsp->GetCubeCplxJoin() = 0.5/dt/c3*pvariable->divs->GetCubeCplxJoin();

  
  
  // 1.4 Solve the Helmholtz equation of pressure.
  Helmholtz(*pvariable->rhsp,*pvariable->p[index],stage);
  
  // 2.1 Get the gradient of pressure.
  sstepper_->Gradient(*pvariable->p[index],*pvariable->dpdxyz[index],0);
  
  // 2.1 Get the 2nd gradient of pressure.
  // !!!!!!! be careful for what dudxyztmp represents now!
  if (boundary_->IsOpenTop()) {
    pvariable->ptmp[index]->GetCubeCplxJoin() = pvariable->dpdxyz[index]->GetCubeCplx(2);
    sstepper_->Gradient(*pvariable->ptmp[index],*pvariable->dudxyztmp[index],0);
  }
  
  // 2.2 Update boundary information for next time step.
  boundary_->Update(stage);
  
  

  
  // 2.3 Apply the correction of pressure.
  const double coeffix = 2*dt*c3;
  pvariable->uvw[index]->GetCubeCplxJoin() = pvariable->suvw->GetCubeCplxJoin()-coeffix*pvariable->dpdxyz[index]->GetCubeCplxJoin();
  
  
  // modified for test!!!
  if (dcp_.isMeanZ()) {
    pvariable->uvw[index]->GetCMatCplx(2).col(0).fill(cx_double(0,0));
  }

  
  
  
  pvariable->uvwnonlinear->GetCubeCplxJoin() = pvariable->uvw[index]->GetCubeCplxJoin();
  pvariable->uvwnonlinear->DftBackward();
  
  
  
  const int rank = mpi_.GetRank();
  const imat& dimension = dcp_.GetDimension();
  
  if (stage == 2) {
    
//    featurelocal_.zeros();
    featureglobal_.zeros();
    
    // max velocity
    featurelocal_(0) = abs(pvariable->uvwnonlinear->GetCubeReal(0)).max();
    featurelocal_(1) = abs(pvariable->uvwnonlinear->GetCubeReal(1)).max();
    featurelocal_(2) = abs(pvariable->uvwnonlinear->GetCubeReal(2)).max();
    
    
    if (dcp_.isRealBottom()) {
      featurelocal_(3) = abs(pvariable->uvwnonlinear->GetCubeReal(0).slice(0)).max();
      featurelocal_(4) = abs(pvariable->uvwnonlinear->GetCubeReal(1).slice(0)).max();
      featurelocal_(5) = abs(pvariable->uvwnonlinear->GetCubeReal(2).slice(0)).max();
    } else {
      featurelocal_(3) = 0;
      featurelocal_(4) = 0;
      featurelocal_(5) = 0;
    }
    
    if (dcp_.isRealTop()) {
      featurelocal_(6) = abs(pvariable->uvwnonlinear->GetCubeReal(0).slice(dimension(10*2+8,rank)-1)).max();
      featurelocal_(7) = abs(pvariable->uvwnonlinear->GetCubeReal(1).slice(dimension(10*2+8,rank)-1)).max();
      featurelocal_(8) = abs(pvariable->uvwnonlinear->GetCubeReal(2).slice(dimension(10*2+8,rank)-1)).max();
    } else {
      featurelocal_(6) = 0;
      featurelocal_(7) = 0;
      featurelocal_(8) = 0;
    }
    
    
    
    // find local cfl without dt
    domain_->GetMaxULen() =
      arma::abs(pvariable->uvwnonlinear->GetCubeReal(0))/domain_->GetGridLenX()+
      arma::abs(pvariable->uvwnonlinear->GetCubeReal(1))/domain_->GetGridLenY()+
      arma::abs(pvariable->uvwnonlinear->GetCubeReal(2))/domain_->GetGridLenZ();
    featurelocal_(19) = domain_->GetMaxULen().max();
    
    
    
    MPI_Allreduce(featurelocal_.memptr(),
                  featureglobal_.memptr(),
                  50,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD);
   
    //
    // rewrite this part!!!
    if (GetDcp().isMeanZ()) {
      
      featureglobal_(10) = real(pvariable->dudxyz[index]->GetCMatCplx(2)(ibottom,0));
      featureglobal_(11) = real(pvariable->dvdxyz[index]->GetCMatCplx(2)(ibottom,0));
      featureglobal_(12) = real(pvariable->dwdxyz[index]->GetCMatCplx(2)(ibottom,0));
      featureglobal_(13) = real(pvariable->dudxyz[index]->GetCMatCplx(2)(itop,0));
      featureglobal_(14) = real(pvariable->dvdxyz[index]->GetCMatCplx(2)(itop,0));
      featureglobal_(15) = real(pvariable->dwdxyz[index]->GetCMatCplx(2)(itop,0));

      // load scalar phase
      featureglobal_(21) = featurelocal_(21);
      featureglobal_(22) = featurelocal_(22);
      featureglobal_(23) = featurelocal_(23);
      featureglobal_(24) = featurelocal_(24);
      featureglobal_(25) = featurelocal_(25);
      featureglobal_(28) = featurelocal_(28);
      featureglobal_(29) = featurelocal_(29);
    }
 
  }
  
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Correct(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::NonlinearVelocity(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"NonlinearVelocity(int)")
  
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  const int index = pvariable->GetIndex();
  

  if (tstepper_->GetAlternationFlag()) {
  
  pvariable->dudxyztmp[index]->GetCubeReal(0) =
  pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->uvwnonlinear->GetCubeReal(0);
  pvariable->dudxyztmp[index]->GetCubeReal(1) =
  pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->uvwnonlinear->GetCubeReal(1);
  pvariable->dudxyztmp[index]->GetCubeReal(2) =
  pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->uvwnonlinear->GetCubeReal(2);
  
  pvariable->dvdxyztmp[index]->GetCubeReal(0) =
  pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->uvwnonlinear->GetCubeReal(0);
  pvariable->dvdxyztmp[index]->GetCubeReal(1) =
  pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->uvwnonlinear->GetCubeReal(1);
  pvariable->dvdxyztmp[index]->GetCubeReal(2) =
  pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->uvwnonlinear->GetCubeReal(2);
  
  pvariable->dwdxyztmp[index]->GetCubeReal(0) =
  pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->uvwnonlinear->GetCubeReal(0);
  pvariable->dwdxyztmp[index]->GetCubeReal(1) =
  pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->uvwnonlinear->GetCubeReal(1);
  pvariable->dwdxyztmp[index]->GetCubeReal(2) =
  pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->uvwnonlinear->GetCubeReal(2);
  
  
  pvariable->dudxyztmp[index]->DftForward();
  pvariable->dvdxyztmp[index]->DftForward();
  pvariable->dwdxyztmp[index]->DftForward();
  
  
  sstepper_->GradientTmp(*pvariable->dudxyztmp[index],*pvariable->dudxyz[index]);
  sstepper_->GradientTmp(*pvariable->dvdxyztmp[index],*pvariable->dvdxyz[index]);
  sstepper_->GradientTmp(*pvariable->dwdxyztmp[index],*pvariable->dwdxyz[index]);
  
  
  pvariable->advuvw->GetCubeCplx(0) =
  pvariable->dudxyz[index]->GetCubeCplx(0)+
  pvariable->dudxyz[index]->GetCubeCplx(1)+
  pvariable->dudxyz[index]->GetCubeCplx(2);
  
  pvariable->advuvw->GetCubeCplx(1) =
  pvariable->dvdxyz[index]->GetCubeCplx(0)+
  pvariable->dvdxyz[index]->GetCubeCplx(1)+
  pvariable->dvdxyz[index]->GetCubeCplx(2);
  
  pvariable->advuvw->GetCubeCplx(2) =
  pvariable->dwdxyz[index]->GetCubeCplx(0)+
  pvariable->dwdxyz[index]->GetCubeCplx(1)+
  pvariable->dwdxyz[index]->GetCubeCplx(2);
    
  sstepper_->Gradient(*pvariable->uvw[index],*pvariable->dudxyz[index],0);
  sstepper_->Gradient(*pvariable->uvw[index],*pvariable->dvdxyz[index],1);
  sstepper_->Gradient(*pvariable->uvw[index],*pvariable->dwdxyz[index],2);
  if (pvariable->uvw.size()>1) {
  pvariable->dudxyz[index]->DftBackward();
  pvariable->dvdxyz[index]->DftBackward();
  pvariable->dwdxyz[index]->DftBackward();
  }
    
  } else {
  
  
  sstepper_->Gradient(*pvariable->uvw[index],*pvariable->dudxyz[index],0);
  sstepper_->Gradient(*pvariable->uvw[index],*pvariable->dvdxyz[index],1);
  sstepper_->Gradient(*pvariable->uvw[index],*pvariable->dwdxyz[index],2);
  pvariable->dudxyz[index]->DftBackward();
  pvariable->dvdxyz[index]->DftBackward();
  pvariable->dwdxyz[index]->DftBackward();
  
  
  pvariable->advuvw->GetCubeReal(0) =
    pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->dudxyz[index]->GetCubeReal(0)+
    pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->dudxyz[index]->GetCubeReal(1)+
    pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->dudxyz[index]->GetCubeReal(2);
  pvariable->advuvw->GetCubeReal(1) =
    pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->dvdxyz[index]->GetCubeReal(0)+
    pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->dvdxyz[index]->GetCubeReal(1)+
    pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->dvdxyz[index]->GetCubeReal(2);
  pvariable->advuvw->GetCubeReal(2) =
    pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->dwdxyz[index]->GetCubeReal(0)+
    pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->dwdxyz[index]->GetCubeReal(1)+
    pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->dwdxyz[index]->GetCubeReal(2);
  
  
  pvariable->advuvw->DftForward();

  }
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"NonlinearVelocity(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::NonlinearPressure(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"NonlinearPressure(int)")
  
  PS_DEBUG_TRACE_LEAVE(name_+"NonlinearPressure(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::NonlinearCorrection(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"NonlinearCorrection(int)")
  
  PS_DEBUG_TRACE_LEAVE(name_+"NonlinearCorrection(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatCarrierPhase::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FlatCarrierPhase > ";
  Phase::Print(stream);
//  stream << std::endl;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatCarrierPhase::Helmholtz(rccube& target, rccube& result, const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Helmholtz(int,zbdy_t,rccube,rccube)")
  
  const int sizez = dcp_.GetSize(2);
  FlatCarrierVariable* pvariable = (FlatCarrierVariable*)variable_;
  
  if (!iscorrectionstep_) {
    
    
    // compute u
    const cmat& evecnrmu = dz2dcpu1_.GetEigenVectorNorm();
    const cmat& evecinvu = dz2dcpu1_.GetEigenVectorInverse();
    const cmat& fixtfwdu = dz2dcpu1_.GetFixForwardTop();
    const cmat& fixbfwdu = dz2dcpu1_.GetFixForwardBottom();
    const cmat& fixtbwdu = dz2dcpu1_.GetFixBackwardTop();
    const cmat& fixbbwdu = dz2dcpu1_.GetFixBackwardBottom();
    const rccube& evalrcpu = dz2dcpu1_.GetGetEigenVectorReciprocal(stage);
    
    pvariable->tmpu1->GetCubeCplxJoin() = target.GetCubeCplx(0);
    for (int iz = 1; iz < sizez-1; ++iz) {
      pvariable->tmpu1->GetCMatCplxJoin().row(iz) +=
      fixbfwdu(iz)*pvariable->tmpu1->GetCMatCplxJoin().row(0      )+
      fixtfwdu(iz)*pvariable->tmpu1->GetCMatCplxJoin().row(sizez-1);
    }
    pvariable->tmpu1->Multiply(evecnrmu,evalrcpu,evecinvu,*pvariable->tmpu2);
    pvariable->tmpu2->Multiply(fixbbwdu,0      );
    pvariable->tmpu2->Multiply(fixtbwdu,sizez-1);
    result.GetCubeCplx(0) = pvariable->tmpu2->GetCubeCplx(0);
    
    
    
    // compute v
    const cmat& evecnrmv = dz2dcpu2_.GetEigenVectorNorm();
    const cmat& evecinvv = dz2dcpu2_.GetEigenVectorInverse();
    const cmat& fixtfwdv = dz2dcpu2_.GetFixForwardTop();
    const cmat& fixbfwdv = dz2dcpu2_.GetFixForwardBottom();
    const cmat& fixtbwdv = dz2dcpu2_.GetFixBackwardTop();
    const cmat& fixbbwdv = dz2dcpu2_.GetFixBackwardBottom();
    const rccube& evalrcpv = dz2dcpu2_.GetGetEigenVectorReciprocal(stage);
    
    pvariable->tmpv1->GetCubeCplxJoin() = target.GetCubeCplx(1);
    for (int iz = 1; iz < sizez-1; ++iz) {
      pvariable->tmpv1->GetCMatCplxJoin().row(iz) +=
      fixbfwdv(iz)*pvariable->tmpv1->GetCMatCplxJoin().row(0      )+
      fixtfwdv(iz)*pvariable->tmpv1->GetCMatCplxJoin().row(sizez-1);
    }
    pvariable->tmpv1->Multiply(evecnrmv,evalrcpv,evecinvv,*pvariable->tmpv2);
    pvariable->tmpv2->Multiply(fixbbwdv,0      );
    pvariable->tmpv2->Multiply(fixtbwdv,sizez-1);
    result.GetCubeCplx(1) = pvariable->tmpv2->GetCubeCplx(0);
    
    
    // compute w
    const cmat& evecnrmw = dz2dcpu3_.GetEigenVectorNorm();
    const cmat& evecinvw = dz2dcpu3_.GetEigenVectorInverse();
    const cmat& fixtfwdw = dz2dcpu3_.GetFixForwardTop();
    const cmat& fixbfwdw = dz2dcpu3_.GetFixForwardBottom();
    const cmat& fixtbwdw = dz2dcpu3_.GetFixBackwardTop();
    const cmat& fixbbwdw = dz2dcpu3_.GetFixBackwardBottom();
    const rccube& evalrcpw = dz2dcpu3_.GetGetEigenVectorReciprocal(stage);
    
    pvariable->tmpw1->GetCubeCplxJoin() = target.GetCubeCplx(2);
    for (int iz = 1; iz < sizez-1; ++iz) {
      pvariable->tmpw1->GetCMatCplxJoin().row(iz) +=
      fixbfwdw(iz)*pvariable->tmpw1->GetCMatCplxJoin().row(0      )+
      fixtfwdw(iz)*pvariable->tmpw1->GetCMatCplxJoin().row(sizez-1);
    }
    pvariable->tmpw1->Multiply(evecnrmw,evalrcpw,evecinvw,*pvariable->tmpw2);
    pvariable->tmpw2->Multiply(fixbbwdw,0      );
    pvariable->tmpw2->Multiply(fixtbwdw,sizez-1);
    result.GetCubeCplx(2) = pvariable->tmpw2->GetCubeCplx(0);
    
    
    
//    const cmat& evecnrm = dz2dcp_.GetEigenVectorNorm();
//    const cmat& evecinv = dz2dcp_.GetEigenVectorInverse();
//    const cmat& fixtfwd = dz2dcp_.GetFixForwardTop();
//    const cmat& fixbfwd = dz2dcp_.GetFixForwardBottom();
//    const cmat& fixtbwd = dz2dcp_.GetFixBackwardTop();
//    const cmat& fixbbwd = dz2dcp_.GetFixBackwardBottom();
//    const rccube& evalrcp = dz2dcp_.GetGetEigenVectorReciprocal(stage);
//
//    for (int iz = 1; iz < sizez-1; ++iz) {
//      target.GetCMatCplxJoin().row(iz) +=
//      fixbfwd(iz)*target.GetCMatCplxJoin().row(0      )+
//      fixtfwd(iz)*target.GetCMatCplxJoin().row(sizez-1);
//    }
//
//    target.Multiply(evecnrm,evalrcp,evecinv,result);
//
//    result.Multiply(fixbbwd,0      );
//    result.Multiply(fixtbwd,sizez-1);
    
    
  } else {
    
    const cmat& evecnrm = dz2dcppressure_.GetEigenVectorNorm();
    const cmat& evecinv = dz2dcppressure_.GetEigenVectorInverse();
//    const cmat& fixtfwd = dz2dcppressure_.GetFixForwardTop();
//    const cmat& fixbfwd = dz2dcppressure_.GetFixForwardBottom();
    const cmat& fixtbwd = dz2dcppressure_.GetFixBackwardTop();
    const cmat& fixbbwd = dz2dcppressure_.GetFixBackwardBottom();
    const rccube& evalrcp = dz2dcppressure_.GetGetEigenVectorReciprocal(stage);
    
    
    //
    target.GetCMatCplxJoin().row(0      ).fill(cx_double(0,0));
    target.GetCMatCplxJoin().row(sizez-1).fill(cx_double(0,0));
    
//    const double dt = tstepper_->GetDt();
//    const double c3 = tstepper_->GetCoefficients(2,stage);
//    target.GetCMatCplxJoin().row(0      ) = pvariable->suvw->GetCMatCplx(2).row(0      )/2/dt/c3;
//    target.GetCMatCplxJoin().row(sizez-1) = pvariable->suvw->GetCMatCplx(2).row(sizez-1)/2/dt/c3;
//
//    for (int iz = 1; iz < sizez-1; ++iz) {
//      target.GetCMatCplxJoin().row(iz) +=
//      fixbfwd(iz)*target.GetCMatCplxJoin().row(0      )+
//      fixtfwd(iz)*target.GetCMatCplxJoin().row(sizez-1);
//    }
    
    //
    target.Multiply(evecnrm,evalrcp,evecinv,result);
    
    result.Multiply(fixbbwd,0      );
    result.Multiply(fixtbwd,sizez-1);
    
    
  }
  

  
  PS_DEBUG_TRACE_LEAVE(name_+"Helmholtz(int,zbdy_t,rccube,rccube)")
  return true;
}
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
