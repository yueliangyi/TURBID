/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_TIME_STEPPER_H_
#define YLY_PS_LVONE_TIME_STEPPER_H_

#include "base.h"
#include "cube_real_cplx.h"
#include "file_ini.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class TimeStepper : public levelbase::Base {
  
public:
  
  struct Configuration : public levelbase::Base::Configuration {
  public:
    double maxcfl = 0;
    int    output_based_on_step = 1;
    
    int    maxstep       = 0;
    int    writeinterval = 0;
    double timestep     = 0;
    int    startstep     = 0;
    
    double max_marching_time = 0;
    double output_interval = 0;
    double initial_time_step = 0;
    int    start_time_index     = 0;
    dvec   cfl_adjust_scale;
    
    bool   writeinitial  = false;
    
    bool Load(levelbase::IniFile& fileini);
    
  protected:
  private:
  }; typedef Configuration cfg_t;
  
  
  
  TimeStepper(const levelbase::mpi_t& mpi,
              const levelbase::dcp_t& ddbt,
              const levelbase::Base::cfg_t* cfg );
  virtual ~TimeStepper(void) { }
  
  
  int GetStep(void) const { return step_; }
  double GetDt(void) const { return dt_; }
  double GetTime(void) const { return time_; }
  double GetCfl(void) const { return cfl_; }
  int GetStage(void) const { return stage_; }
  int GetStartStep(void) const { return ((cfg_t*)cfg_)->startstep; }
  int GetStartTimeIndex(void) const { return ((cfg_t*)cfg_)->start_time_index; }
  int GetOutputTimeIndex(void) const { return output_time_index_; }

  void SetDt(const double dt) { dt_ = dt; }
  
  
  const dmat& GetCoefficients(void) const { return coefficients_; }
  double GetCoefficients(const int id, const int stage) const { return coefficients_(id,stage); }
  
  
  bool isForward(void) const { return isforward_; }
  bool isOutput(void) const { return isoutput_; }
  bool isOutputOnStep(void) const { return ((cfg_t*)cfg_)->output_based_on_step; }
  
  
  bool GetAlternationFlag(void) const { return alternationflag_; }
  void SetAlternationFlag(void) {
    if (alternationflag_) {
      alternationflag_ = false;
    } else {
      alternationflag_ = true;
    }
    // better to do broadcast
    mpi_.BcastBool(alternationflag_);
  }
  
  
  bool Update(const double maxulen = 0);
  
  bool Update2(const double maxulen = 0);
  
  
  virtual bool Configure(void);
  virtual void Print(ostringstream& stream) const;
  
protected:
  
  
  
  bool   isforward_ = false;
  bool   isoutput_ = false;
  
  int    step_     = 0;
  double dt_       = 0;
  double time_     = 0;
  double cfl_      = 0;
  int    stage_    = 0;
  
  int output_time_index_ = 0;
  
  dmat coefficients_;
  
  bool   alternationflag_ = false;

  string name_ = "TimeStepper::";

private:
  
};
  
  
  
struct TypeTimeStepper {
  TimeStepper::cfg_t* cfg = NULL;
  TimeStepper       * obj = NULL;
}; typedef TypeTimeStepper tstepper_t;
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_TIME_STEPPER_H_
