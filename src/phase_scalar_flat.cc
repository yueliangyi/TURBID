/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "phase_scalar_flat.h"
#include <iomanip>
#include <stdexcept>
#include "variable_scalar_flat.h"
#include "boundary_scalar_flat.h"
#include <sys/stat.h>
#include <cstdlib>


namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::Configuration::Load(IniFile& fileini) {
  PS_DEBUG_TRACE_ENTER("FlatScalarPhase::Configuration::Load(IniFile)")
  
  Phase::Configuration::Load(fileini);

  reynolds = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_REYNOLDS_NUMBER);
  froude = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_FROUDE_NUMBER);
  averagedconcentration = fileini.GetValueDouble(SECTION_PHASE_CARRIER,KEY_SPC_AVERAGED_CONCENTRATION);
  
  schmidt = fileini.GetValueDouble(SECTION_PHASE_SCALAR,KEY_SPS_SCHMIDT_NUMBER);
  settlingvelocity = fileini.GetValueDouble(SECTION_PHASE_SCALAR,KEY_SPS_SETTLING_VELOCITY);
  erosionrate = fileini.GetValueDouble(SECTION_PHASE_SCALAR,KEY_SPS_EROSION_RATE);

  
  
  string bdystr = ConvertCase(fileini.GetValueString(SECTION_PHASE_SCALAR,KEY_SPS_VERTICAL_BOUNDARY_TYPE),1);
  if (bdystr == "ZEROFLUX") { zboundarytype = kBoundaryZeroFlux; }
  else if (bdystr == "EROSIONDEPOSITION") {
    zboundarytype = kBoundaryErosionDeposition;
    criticalshearstress = fileini.GetVectorDouble(SECTION_PHASE_SCALAR,KEY_SPS_CRITICAL_SHEAR_STRESS);
  } else {
    throw std::runtime_error(string("FlatScalarPhase::Configuration::")+
                             "Load(IniFile) > Unsupported Boundary Type!");
  }
  
  PS_DEBUG_TRACE_LEAVE("FlatScalarPhase::Configuration::Load(IniFile)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatScalarPhase::FlatScalarPhase(const mpi_t& mpi,
                                   const dcp_t& dcp,
                                   Domain      * domain,
                                   SpaceStepper* sstepper,
                                   TimeStepper * tstepper,
                                   const Base::cfg_t* cfg,
                                   const int index) :
  Phase(mpi,dcp,domain,sstepper,tstepper,cfg,index) {
  PS_DEBUG_TRACE_ENTER(name_+"FlatScalarPhase(mpi_t,dcp_t,cfg_t)")
    
    
    variable_ = new FlatScalarVariable(mpi,dcp,domain,sstepper,tstepper,index_);
    boundary_ = new FlatScalarBoundary(mpi,dcp,domain,sstepper,tstepper,variable_,index_);
    
    
    // in z-pencile
    const int rank = mpi_.GetRank();
    const imat& dimension = dcp_.GetDimension();
    localmass_ = arma::zeros<dmat>(dimension(10*1+2,rank),6);
    globalmass_ = arma::zeros<dmat>(dimension(10*1+2,rank),6);

    
    
  PS_DEBUG_TRACE_LEAVE(name_+"FlatScalarPhase(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
FlatScalarPhase::~FlatScalarPhase(void) {
  if (variable_ != NULL) { delete variable_; variable_ = NULL; }
  if (boundary_ != NULL) { delete boundary_; boundary_ = NULL; }
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::Save(const string subfolder, const ioctrl_t level) {
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
bool FlatScalarPhase::Load(const string subfolder, const ioctrl_t level) {
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
bool FlatScalarPhase::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  cfg_t* pcfg = (cfg_t*) cfg_;
  FlatScalarVariable* pvariable = (FlatScalarVariable*)variable_;
  
  
  // Construct the stage coefficients for Helmholtz solution.
  const double reynolds = pcfg->reynolds;
  const double schmidt = pcfg->schmidt;
  const double settlingvelocity = pcfg->settlingvelocity;
  const double dt  = tstepper_->GetDt();
  
  dvec stepcoefficient(3);
  stepcoefficient(0) = -reynolds*schmidt/dt/tstepper_->GetCoefficients(2,0);
  stepcoefficient(1) = -reynolds*schmidt/dt/tstepper_->GetCoefficients(2,1);
  stepcoefficient(2) = -reynolds*schmidt/dt/tstepper_->GetCoefficients(2,2);

  dvec bdycoefficient(4);
  bdycoefficient(0) = reynolds*schmidt*settlingvelocity;
  bdycoefficient(1) = 1;
  bdycoefficient(2) = reynolds*schmidt*settlingvelocity;
  bdycoefficient(3) = 1;
  bdycoefficient = bdycoefficient*BOUNDARY_COEF_SCALE;
  dz2dcp_.Construct(sstepper_,bdycoefficient,stepcoefficient,1,true);
  
  
  
  variable_->Configure();
  boundary_->Configure();
  
  
  

  const int index = variable_->GetIndex();
  
  
  pvariable->uvwnonlinear->GetCubeCplxJoin() = pvariable->uvw[index]->GetCubeCplxJoin();
  pvariable->uvwnonlinear->DftBackward();
  
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
bool FlatScalarPhase::Update(void) {
  
  cfg_t* pcfg = (cfg_t*) cfg_;
  const double reynolds = pcfg->reynolds;
  const double schmidt = pcfg->schmidt;
  const double dt  = tstepper_->GetDt();
  dvec stepcoefficient(3);
  stepcoefficient(0) = -reynolds*schmidt/dt/tstepper_->GetCoefficients(2,0);
  stepcoefficient(1) = -reynolds*schmidt/dt/tstepper_->GetCoefficients(2,1);
  stepcoefficient(2) = -reynolds*schmidt/dt/tstepper_->GetCoefficients(2,2);
  dz2dcp_.Assemble(sstepper_,stepcoefficient);
  
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::Predict(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Predict(int)")
  
  iscorrectionstep_ = false;
  
  
  cfg_t* pcfg = (cfg_t*) cfg_;
  FlatScalarVariable* pvariable = (FlatScalarVariable*)variable_;
  const int index = pvariable->GetIndex();
  
  
  pvariable->uvw[index]->GetCubeCplxJoin() = pvariable->uvw[0]->GetCubeCplxJoin();
  
  if (dcp_.isMeanZ()) {
    Domain::cfg_t* pcfgdomain = (Domain::cfg_t*) domain_->GetCfg();
    switch (pcfgdomain->slope_direction) {
      case 1: // x1 direction
        pvariable->uvw[index]->GetCMatCplx(0).col(0) =
        pvariable->uvw[index]->GetCMatCplx(0).col(0)+pcfg->settlingvelocity*sin(pcfgdomain->slope);
        break;
        
      case 2: // x2 direction
        pvariable->uvw[index]->GetCMatCplx(1).col(0) =
        pvariable->uvw[index]->GetCMatCplx(1).col(0)+pcfg->settlingvelocity*sin(pcfgdomain->slope);
        break;
        
      default:
        break;
    }
    pvariable->uvw[index]->GetCMatCplx(2).col(0) =
    pvariable->uvw[index]->GetCMatCplx(2).col(0)-pcfg->settlingvelocity*cos(pcfgdomain->slope);
  }
  
  pvariable->uvwnonlinear->GetCubeCplxJoin() = pvariable->uvw[index]->GetCubeCplxJoin();
  pvariable->uvwnonlinear->DftBackward();
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Predict(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::Correct(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Correct(int)")
  
  iscorrectionstep_ = true;
  
  
  cfg_t* pcfg = (cfg_t*) cfg_;
  FlatScalarVariable* pvariable = (FlatScalarVariable*)variable_;
  const int index = pvariable->GetIndex();
  const double dt = tstepper_->GetDt();
  const double c1 = tstepper_->GetCoefficients(0,stage);
  const double c2 = tstepper_->GetCoefficients(1,stage);
  const double c3 = tstepper_->GetCoefficients(2,stage);
  
  
  // advection
  NonlinearPressure(stage);
  
  
  
  
  // update H
  pvariable->hphi->GetCubeCplxJoin() =
    c1*pvariable->hphi->GetCubeCplxJoin()-
    dt*pvariable->advphi->GetCubeCplxJoin();
  
  
  //
  pvariable->lapphi->GetCMatCplx(0) =
    1*sstepper_->GetCmDz2()              *pvariable->p[index]->GetCMatCplx(0)+
    1*sstepper_->GetCmDx2().GetCMatCplx()%pvariable->p[index]->GetCMatCplx(0)+
    1*sstepper_->GetCmDy2().GetCMatCplx()%pvariable->p[index]->GetCMatCplx(0);
  
  // 4.1 Compute the term vec{rhs}^{m0}.
  const double redtc3 = -pcfg->reynolds*pcfg->schmidt/dt/c3;
  pvariable->rhsphi->GetCubeCplxJoin() =
    redtc3*c2*pvariable->hphi->GetCubeCplxJoin()+
    redtc3*pvariable->p[index]->GetCubeCplxJoin()-
    pvariable->lapphi->GetCubeCplxJoin();
  
  
  // 5.1 Boundary condition.
  boundary_->Apply(stage);
  
  
  switch (pcfg->zboundarytype) {
    case kBoundaryZeroFlux:
      break;
      
    case kBoundaryErosionDeposition:
    {
      double tau_crit = pcfg->criticalshearstress(0);
//      const double susp_coef = 5.4464E-07;
      const double susp_coef = pcfg->erosionrate;
      const size_t ibottom = 0;
      Domain::cfg_t* pcfgdomain = (Domain::cfg_t*) domain_->GetCfg();
      /////////////
      // dimensionless
      pvariable->ptmp[index]->GetCubeReal(0) =
      sqrt(pvariable->dudxyz[0]->GetCubeReal(2)%pvariable->dudxyz[0]->GetCubeReal(2)+
           pvariable->dvdxyz[0]->GetCubeReal(2)%pvariable->dvdxyz[0]->GetCubeReal(2))/pcfg->reynolds;
      
      // NOTE: the '&' is crucial!
      
      pvariable->ptmp[index]->GetCubeReal(0).for_each([tau_crit,susp_coef,pcfg](arma::cube::elem_type& val) {
        if (val < tau_crit) {
          val = 0;
        } else {
          val = susp_coef*(val/tau_crit-1);
        }
      });
      
      pvariable->ptmp[index]->DftForward();
      
      pvariable->rhsphi->GetCMatCplxJoin().row(ibottom) =
      pvariable->ptmp[index]->GetCMatCplxJoin().row(ibottom)-
      pcfg->settlingvelocity*cos(pcfgdomain->slope)*
      pvariable->p[index]->GetCMatCplxJoin().row(ibottom)*
      pcfg->averagedconcentration;
      pvariable->rhsphi->GetCMatCplxJoin().row(ibottom) =
      -1*pvariable->rhsphi->GetCMatCplxJoin().row(ibottom)*
      pcfg->reynolds*pcfg->schmidt/
      pcfg->averagedconcentration*
      BOUNDARY_COEF_SCALE;
      
      // collect erosion and deposition flux
      featurelocal_(28) = real(pvariable->ptmp[index]->GetCMatCplxJoin()(0,0));
      featurelocal_(29) = real(pcfg->settlingvelocity*cos(pcfgdomain->slope)*
                               pvariable->p[index]->GetCMatCplxJoin()(0,0)*
                               pcfg->averagedconcentration);
      
      /////////////
    }
      break;
      
    default:
      break;
  }
  
  
  
  
  
  
  // ---------------------------------------------------------------------------
  // 6. Solve the velocity Helmholtz equation for the intermediate velocity
  // ---------------------------------------------------------------------------
  
  Helmholtz(*pvariable->rhsphi,*pvariable->p[index],stage);
 
  
//  //
//  pvariable->p[index]->DftBackward();
//  pvariable->p[index]->GetCubeReal(0).for_each([](arma::cube::elem_type& val) {
//    if (val < 0) {
//      val = 0;
//    }
//  });
//  pvariable->p[index]->DftForward();
//  //
  
  // remove the negative points
  // reference:
  // J. Bartnicki, A simple filtering procedure for removing negative values
  // from numerical solutions of the advection equation, In Environmental
  // Software, Volume 4, Issue 4, 1989, Pages 187-201, ISSN 0266-9838
  // https://doi.org/10.1016/0266-9838(89)90025-7.
  
  const int rank = mpi_.GetRank();
  const imat& dimension = dcp_.GetDimension();

  // get the the concentration in real space
  pvariable->p[index]->DftBackward();
  
  
  // do iteration
  for (int it = 0; it < 5; ++it) {
  
  localmass_.fill(0);
  globalmass_.fill(0);
  for (int indz = 0; indz < dimension(10*2+8,rank); ++indz) {
    for (int indy = 0; indy < dimension(10*2+7,rank); ++indy) {
      for (int indx = 0; indx < dimension(10*2+6,rank); ++indx) {
        if (pvariable->p[index]->GetCubeReal(0)(indx,indy,indz) > 0) {
          localmass_(dimension(10*2+2,rank)+indz-1,0) += pvariable->p[index]->GetCubeReal(0)(indx,indy,indz);
          localmass_(dimension(10*2+2,rank)+indz-1,1) += 1;
        } else if (pvariable->p[index]->GetCubeReal(0)(indx,indy,indz) < 0) {
          localmass_(dimension(10*2+2,rank)+indz-1,4) += pvariable->p[index]->GetCubeReal(0)(indx,indy,indz);
          localmass_(dimension(10*2+2,rank)+indz-1,5) += 1;
        } else {
          localmass_(dimension(10*2+2,rank)+indz-1,2) += pvariable->p[index]->GetCubeReal(0)(indx,indy,indz);
          localmass_(dimension(10*2+2,rank)+indz-1,3) += 1;
        }
      }
    }
  }

  MPI_Allreduce(localmass_.memptr(),
                globalmass_.memptr(),
                6*dimension(10*1+2,rank),
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);
  
  // apply the filter for removing the negative points
  for (int indz = 0; indz < dimension(10*2+8,rank); ++indz) {
    if (globalmass_(dimension(10*2+2,rank)+indz-1,0)+
        globalmass_(dimension(10*2+2,rank)+indz-1,2)+
        globalmass_(dimension(10*2+2,rank)+indz-1,4) < 0) {
      pvariable->p[index]->GetCubeReal(0).slice(indz).fill(0);
      continue;
    }
    for (int indy = 0; indy < dimension(10*2+7,rank); ++indy) {
      for (int indx = 0; indx < dimension(10*2+6,rank); ++indx) {
        if (pvariable->p[index]->GetCubeReal(0)(indx,indy,indz) > 0) {
          pvariable->p[index]->GetCubeReal(0)(indx,indy,indz) +=
          globalmass_(dimension(10*2+2,rank)+indz-1,4)/globalmass_(dimension(10*2+2,rank)+indz-1,1);
        } else {
          pvariable->p[index]->GetCubeReal(0)(indx,indy,indz) = 0;
        }
      }
    }
  }
    
  } // end of iteration
  
  // transform back to spectral space
  pvariable->p[index]->DftForward();
  
  
  
 
  // Update output info
  const dvec& gridz = domain_->GetGridZ();
  dmat avPhi = arma::trapz(gridz,real(pvariable->p[index]->GetCMatCplx(0).col(0)))/2;
  featurelocal_(21) = avPhi(0);
  Domain::cfg_t* pcfgdomain = (Domain::cfg_t*) domain_->GetCfg();
  const int step = tstepper_->GetStep();
  
  if (stage==2 && pcfgdomain->slope_direction>0) {
    
    // mean velocity in u1 and u2
    if (GetDcp().isMeanZ()) {
      dmat avCurrentSpeed = arma::trapz(gridz,real(pvariable->uvw[index]->GetCMatCplx(0).col(0)))/2;
      featurelocal_(22) = avCurrentSpeed(0);
      avCurrentSpeed = arma::trapz(gridz,real(pvariable->uvw[index]->GetCMatCplx(1).col(0)))/2;
      featurelocal_(24) = avCurrentSpeed(0);
    }
    
    // get flux of phi every 10 steps
    featurelocal_(23) = 0;
    featurelocal_(25) = 0;
    if (step%10 == 0) {
      
      pvariable->ptmp[index]->GetCubeReal(0) =
      pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->p[index]->GetCubeReal(0);
      pvariable->ptmp[index]->DftForward();
      
      if (GetDcp().isMeanZ()) {
        dmat avPhiFlux = arma::trapz(gridz,real(pvariable->ptmp[index]->GetCMatCplx(0).col(0)))/2;
        featurelocal_(23) = avPhiFlux(0);
      }
      
      pvariable->ptmp[index]->GetCubeReal(0) =
      pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->p[index]->GetCubeReal(0);
      pvariable->ptmp[index]->DftForward();
      
      if (GetDcp().isMeanZ()) {
        dmat avPhiFlux = arma::trapz(gridz,real(pvariable->ptmp[index]->GetCMatCplx(0).col(0)))/2;
        featurelocal_(25) = avPhiFlux(0);
      }
      
    }
    
  }
  
  
  
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Correct(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::NonlinearVelocity(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"NonlinearVelocity(int)")
  
  PS_DEBUG_TRACE_LEAVE(name_+"NonlinearVelocity(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::NonlinearPressure(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"NonlinearPressure(int)")
  
  FlatScalarVariable* pvariable = (FlatScalarVariable*)variable_;
  const int index = pvariable->GetIndex();
  
  
  if (tstepper_->GetAlternationFlag()) {
//  const int step = tstepper_->GetStep();
//  if (step%2 != 0) {
  
  pvariable->ptmp[index]->GetCubeCplxJoin() = pvariable->p[index]->GetCubeCplxJoin();
  pvariable->ptmp[index]->DftBackward();
  
  pvariable->dvdxyz[index]->GetCubeReal(0) =
  pvariable->ptmp[index]->GetCubeReal(0)%pvariable->uvwnonlinear->GetCubeReal(0);
  pvariable->dvdxyz[index]->GetCubeReal(1) =
  pvariable->ptmp[index]->GetCubeReal(0)%pvariable->uvwnonlinear->GetCubeReal(1);
  pvariable->dvdxyz[index]->GetCubeReal(2) =
  pvariable->ptmp[index]->GetCubeReal(0)%pvariable->uvwnonlinear->GetCubeReal(2);

  pvariable->dvdxyz[index]->DftForward();
  
  sstepper_->GradientTmp(*pvariable->dvdxyz[index],*pvariable->dudxyz[index]);
  
  pvariable->advphi->GetCMatCplx(0) =
    pvariable->dudxyz[index]->GetCMatCplx(0)+
    pvariable->dudxyz[index]->GetCMatCplx(1)+
    pvariable->dudxyz[index]->GetCMatCplx(2);

  } else {

  sstepper_->Gradient(*pvariable->p[index],*pvariable->dudxyz[index],0);
  pvariable->dudxyz[index]->DftBackward();
  
  pvariable->advphi->GetCubeReal(0) =
    pvariable->uvwnonlinear->GetCubeReal(0)%pvariable->dudxyz[index]->GetCubeReal(0)+
    pvariable->uvwnonlinear->GetCubeReal(1)%pvariable->dudxyz[index]->GetCubeReal(1)+
    pvariable->uvwnonlinear->GetCubeReal(2)%pvariable->dudxyz[index]->GetCubeReal(2);
    
  pvariable->advphi->DftForward();
    
  }
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"NonlinearPressure(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::NonlinearCorrection(const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"NonlinearCorrection(int)")
  
  PS_DEBUG_TRACE_LEAVE(name_+"NonlinearCorrection(int)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void FlatScalarPhase::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "FlatScalarPhase > ";
  Phase::Print(stream);
//  stream << std::endl;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool FlatScalarPhase::Helmholtz(rccube& target, rccube& result, const int stage) {
  PS_DEBUG_TRACE_ENTER(name_+"Helmholtz(int,zbdy_t,rccube,rccube)")
  
  const int sizez = dcp_.GetSize(2);
  
  if (!iscorrectionstep_) {
  } else {
    
    const cmat& evecnrm = dz2dcp_.GetEigenVectorNorm();
    const cmat& evecinv = dz2dcp_.GetEigenVectorInverse();
    const cmat& fixtfwd = dz2dcp_.GetFixForwardTop();
    const cmat& fixbfwd = dz2dcp_.GetFixForwardBottom();
    const cmat& fixtbwd = dz2dcp_.GetFixBackwardTop();
    const cmat& fixbbwd = dz2dcp_.GetFixBackwardBottom();
    const rccube& evalrcp = dz2dcp_.GetGetEigenVectorReciprocal(stage);
    
    for (int iz = 1; iz < sizez-1; ++iz) {
      target.GetCMatCplxJoin().row(iz) +=
      fixbfwd(iz)*target.GetCMatCplxJoin().row(0      )+
      fixtfwd(iz)*target.GetCMatCplxJoin().row(sizez-1);
    }
    
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
