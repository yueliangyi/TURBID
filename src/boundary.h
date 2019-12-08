/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_BOUNDARY_H_
#define YLY_PS_LVONE_BOUNDARY_H_

#include "domain.h"
#include "variable.h"
#include "stepper_time.h"
#include "stepper_space.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class Boundary : public levelbase::Base {
  
public:
  
  Boundary(const levelbase::mpi_t& mpi,
           const levelbase::dcp_t& dcp,
           Domain      * domain,
           SpaceStepper* sstepper,
           TimeStepper * tstepper,
           Variable    * variable,
           const int index);
  
  virtual ~Boundary(void) { }
  
  virtual bool Update(const int stage) = 0;
  virtual bool Apply (const int stage) = 0;
  
  virtual bool Configure(void) { return true; }
  virtual void Print(ostringstream& stream) const;
  
  bool SetFlagOpenTop(const dvec& zboundarytype);
  bool IsOpenTop(void) const { return isopentop_; };
  void SetIterationCoefficient(const double itercoef) {
    iterationcoefficient_ = itercoef;
  }
    
protected:
  
  Domain      * domain_   = NULL;
  SpaceStepper* sstepper_ = NULL;
  TimeStepper * tstepper_ = NULL;
  Variable    * variable_ = NULL;
  
  const int index_ = 0;
  bool isopentop_ = false;
  double iterationcoefficient_ = 1;
  
  string name_ = "Boundary::";

private:
  
};
  
  
  
struct TypeBoundary {
  Boundary::cfg_t* cfg = NULL;
  Boundary       * obj = NULL;
}; typedef TypeBoundary boundary_t;
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_BOUNDARY_H_
