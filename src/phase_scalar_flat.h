/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_FLAT_SCALAR_PHASE_H_
#define YLY_PS_LVONE_FLAT_SCALAR_PHASE_H_

#include "phase.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FlatScalarPhase : public Phase {
  
public:
  
  // Type of driven force.
  enum ScalarBounaryType {
    kBoundaryZeroFlux = 0,
    kBoundaryErosionDeposition   = 1,
  }; typedef ScalarBounaryType sbdy_t;
  
  struct Configuration : public Phase::Configuration {
  public:
    double schmidt = 1;
    double settlingvelocity = 1;
    double erosionrate = 1;
    dvec criticalshearstress;
    sbdy_t zboundarytype = kBoundaryZeroFlux;
    bool Load(levelbase::IniFile& fileini);
  protected:
  private:
  }; typedef Configuration cfg_t;
  
  
  
  FlatScalarPhase(const levelbase::mpi_t& mpi,
                   const levelbase::dcp_t& ddbt,
                   Domain      * domain,
                   SpaceStepper* sstepper,
                   TimeStepper * tstepper,
                   const levelbase::Base::cfg_t* cfg,
                   const int index);
  ~FlatScalarPhase(void);
  
  
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
  
  
  // scalar mass M1 (positive), M2 (0), M3 (negative)
  dmat localmass_;
  dmat globalmass_;

  
  string name_ = "FlatScalarPhase::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_FLAT_SCALAR_PHASE_H_
