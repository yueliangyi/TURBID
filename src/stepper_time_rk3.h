/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_TIME_STEPPER_RK3_H_
#define YLY_PS_LVONE_TIME_STEPPER_RK3_H_

#include "stepper_time.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class RK3TimeStepper : public TimeStepper {
  
public:
  
  
  RK3TimeStepper(const levelbase::mpi_t& mpi,
                 const levelbase::dcp_t& ddbt,
                 const levelbase::Base::cfg_t* cfg);
  
  
  
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  
  bool Configure(void) {
    TimeStepper::Configure();
    return true;
  }
  
  
  void Print(ostringstream& stream) const;
  
protected:
  
  string name_ = "RK3TimeStepper::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_TIME_STEPPER_RK3_H_
