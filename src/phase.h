/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_PHASE_H_
#define YLY_PS_LVONE_PHASE_H_

#include "base.h"
#include "cube_real_cplx.h"
#include "domain.h"
#include "stepper_time.h"
#include "stepper_space.h"
#include "variable.h"
#include "boundary.h"
#include "file_ini.h"
#include "cmdz2_decomposition.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
  
class Phase : public levelbase::Base {
  
public:
  
  
  struct Configuration : public levelbase::Base::Configuration {
  public:
    double reynolds = 0;
    double froude   = 0;
    double specific_gravity = 2.65;
    double averagedconcentration = 0.001;
    dvec   rossby;
    
    bool Load(levelbase::IniFile& fileini);
    
  protected:
  private:
  }; typedef Configuration cfg_t;
  
  
  
  Phase(const levelbase::mpi_t& mpi,
        const levelbase::dcp_t& ddbt,
        Domain      * domain,
        SpaceStepper* sstepper,
        TimeStepper * tstepper,
        const levelbase::Base::cfg_t* cfg,
        const int index);
  virtual ~Phase(void);
  
  
 
  const dvec& GetFeatureGlobal(void) const { return featureglobal_; }
  
  
  virtual bool Configure(void) = 0;
  
  
  virtual bool Predict(const int stage) = 0;
  virtual bool Correct(const int stage) = 0;
  
  virtual bool Update(void) = 0;
  
  virtual void Print(ostringstream& stream) const;
  
  
  Variable* GetVariable(void) { return variable_; }
  Boundary* GetBoundary(void) { return boundary_; }

  
protected:
  
  
  
  // Function for Helmholtz equation solution.
  virtual bool Helmholtz(levelbase::rccube& target,
                         levelbase::rccube& result,
                         const int stage) = 0;
  
  
  
  Domain      * domain_   = NULL;
  SpaceStepper* sstepper_ = NULL;
  TimeStepper * tstepper_ = NULL;
  Variable    * variable_ = NULL;
  Boundary    * boundary_ = NULL;
  
  const int index_ = 0;
  
  
  // Common decomposition of cmdz2 for direct solution
  // for velocity or concentration
  dz2dc_t dz2dcp_;
//  dz2dr_t dz2dcp_;
  
  bool iscorrectionstep_ = false;
  
  
  // max u,v,w...
  static dvec featurelocal_;
  static dvec featureglobal_;
  
  
  
  
  string name_ = "Phase::";

private:
  
};
  
  
  
struct TypePhase {
  Phase::cfg_t* cfg = NULL;
  Phase       * obj = NULL;
}; typedef TypePhase phase_t;
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_PHASE_H_
