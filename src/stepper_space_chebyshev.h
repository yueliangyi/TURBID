/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_CHEBYSHEV_SPACE_STEPPER_H_
#define YLY_PS_LVONE_CHEBYSHEV_SPACE_STEPPER_H_

#include "stepper_space.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class ChebyshevSpaceStepper : public SpaceStepper {
  
public:
  
  ChebyshevSpaceStepper(const levelbase::mpi_t& mpi,
                        const levelbase::dcp_t& ddbt,
                        Domain* domain,
                        const levelbase::Base::cfg_t* cfg = NULL);
  
  
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  
  bool Configure(void);
  void Print(ostringstream& stream) const;
  
protected:
  
  
  string name_ = "ChebyshevSpaceStepper::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_CHEBYSHEV_SPACE_STEPPER_H_
