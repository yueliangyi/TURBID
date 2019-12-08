/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_VARIABLE_H_
#define YLY_PS_LVONE_VARIABLE_H_

#include "domain.h"
#include "stepper_time.h"
#include "stepper_space.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class Variable : public levelbase::Base {
  
public:
  
  
  Variable(const levelbase::mpi_t& mpi,
           const levelbase::dcp_t& ddbt,
           Domain      * domain,
           SpaceStepper* sstepper,
           TimeStepper * tstepper,
           const int index);
  virtual ~Variable(void);
  
  int GetIndex(void) const { return index_; }
  
  
  
  virtual bool Configure(void) { return true; }
  virtual void Print(ostringstream& stream) const;
  
  
  
  //
  // Variables (part of them) for different phases, including
  // 1. velocity and its derivatives
  // 2. pressure and its derivatives
  // 3. concentration and its derivatives
  // Note that, the variable name "p" is shared. For the carrier phase, it means
  // pressure. For the sediment phase, it means the concentration (the first
  // character of phi). Be careful when using these variables.
  //
//  static std::vector<levelbase::rccube*> u;
//  static std::vector<levelbase::rccube*> dudx;
//  static std::vector<levelbase::rccube*> dudy;
//  static std::vector<levelbase::rccube*> dudz;
//  
//  static std::vector<levelbase::rccube*> v;
//  static std::vector<levelbase::rccube*> dvdx;
//  static std::vector<levelbase::rccube*> dvdy;
//  static std::vector<levelbase::rccube*> dvdz;
//  
//  static std::vector<levelbase::rccube*> w;
//  static std::vector<levelbase::rccube*> dwdx;
//  static std::vector<levelbase::rccube*> dwdy;
//  static std::vector<levelbase::rccube*> dwdz;
//  
//  static std::vector<levelbase::rccube*> p;
//  static std::vector<levelbase::rccube*> dpdx;
//  static std::vector<levelbase::rccube*> dpdy;
//  static std::vector<levelbase::rccube*> dpdz;
  
  
  static std::vector<levelbase::rccube*> uvw;
  static std::vector<levelbase::rccube*> dudxyz;
  static std::vector<levelbase::rccube*> dvdxyz;
  static std::vector<levelbase::rccube*> dwdxyz;
  
  static std::vector<levelbase::rccube*> p;
  static std::vector<levelbase::rccube*> dpdxyz;
  
  
  static std::vector<levelbase::rccube*> dudxyztmp;
  static std::vector<levelbase::rccube*> dvdxyztmp;
  static std::vector<levelbase::rccube*> dwdxyztmp;
  static std::vector<levelbase::rccube*> ptmp;
  
  
  static std::vector<levelbase::rccube*> botsstres;

  
protected:
  
  
  Domain      * domain_   = NULL;
  SpaceStepper* sstepper_ = NULL;
  TimeStepper * tstepper_ = NULL;
  
  
  
  // Index for the position in variable vectors.
  // Be careful that the index is shared within all variables, meaning that the
  // variable should have the same size and sequence.
  int index_ = 0;
  
  
  string name_ = "Variable::";

private:
  
};
  
  
  
struct TypeVariable {
  Variable::cfg_t* cfg = NULL;
  Variable       * obj = NULL;
}; typedef TypeVariable variable_t;
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_VARIABLE_H_
