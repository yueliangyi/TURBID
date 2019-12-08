/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_FLAT_CARRIER_PHASE_H_
#define YLY_PS_LVONE_FLAT_CARRIER_PHASE_H_

#include "phase.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FlatCarrierPhase : public Phase {
  
public:
  
  // Type of driven force.
  enum CarrierForceType {
    kForceNone      = 0,
    kForceConstant  = 1,
    kForceOscillatory = 2,
    kForcePulsatile = 3,
    kForceDamping = 4,
  }; typedef CarrierForceType cforce_t;
  
  struct Configuration : public Phase::Configuration {
  public:
    dvec zboundarytypeu1;
    dvec zboundarytypeu2;
    dvec zboundarytypeu3;
    dvec zboundarytypep;
    cforce_t forcetype = kForceConstant;
    ivec forcedirection;
    dvec forceamplitude;
    
    dvec dampingcoefficient;
    
    bool affected_by_scalar_phase = true;
    bool Load(levelbase::IniFile& fileini);
  protected:
  private:
  }; typedef Configuration cfg_t;
  
  
  
  FlatCarrierPhase(const levelbase::mpi_t& mpi,
                   const levelbase::dcp_t& ddbt,
                   Domain      * domain,
                   SpaceStepper* sstepper,
                   TimeStepper * tstepper,
                   const levelbase::Base::cfg_t* cfg,
                   const int index);
  ~FlatCarrierPhase(void);
  
  
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  
  bool Configure(void);
  
  
  bool Predict(const int stage);
  bool Correct(const int stage);

  
  bool Update(void);
  
  
  void Print(ostringstream& stream) const;
  
protected:
  
  
  bool NonlinearVelocity(const int stage);
  bool NonlinearPressure(const int stage);
  bool NonlinearCorrection(const int stage);
  
  bool Helmholtz(levelbase::rccube& target,
                 levelbase::rccube& result,
                 const int stage);
  
  
  dz2dc_t dz2dcpu1_;
  dz2dc_t dz2dcpu2_;
  dz2dc_t dz2dcpu3_;
  dz2dc_t dz2dcppressure_;
  
  
//  levelbase::rdcube* velocitylhs_ = NULL;
//  levelbase::rdcube* velocityrhs_ = NULL;
  
  double bulkvelocity_ = 1;
  double applyforce_ = 0.01;
  
  
  cmat chebbackward_;
  cmat chebforward_;
  cmat chebcoemat_;
  cmat chebsolcoemat_;
  
  
  string name_ = "FlatCarrierPhase::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_FLAT_CARRIER_PHASE_H_
