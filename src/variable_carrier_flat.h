/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_FLAT_CARRIER_VARIABLE_H_
#define YLY_PS_LVONE_FLAT_CARRIER_VARIABLE_H_

#include "variable.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FlatCarrierVariable : public Variable {
  
public:
  
  
  FlatCarrierVariable(const levelbase::mpi_t& mpi,
                      const levelbase::dcp_t& ddbt,
                      Domain      * domain,
                      SpaceStepper* sstepper,
                      TimeStepper * tstepper,
                      const int index);
  ~FlatCarrierVariable(void);
  
  //
  // IO operations for input and output.
  // Though these operations are similar to functions of file class, they have
  // a quite different property thus there is no hierarchy between them.
  //
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  
  bool Configure(void);
  void Print(ostringstream& stream) const;
  
  
  
  
  // Varialbes defined for the carrier phase, including
  // 1. advection A[u]
  // 2. intermediate variable H
  // 3. intermediate velocity u^*
  // 4. laplacian L[u]
  // 5. right-hand side term for u^*
  // 6. divergence of u^*
  // 7. right-hand side term for pressrue p
  
  
//  levelbase::rccube* dealiasingu = NULL;
//  levelbase::rccube* dealiasingv = NULL;
//  levelbase::rccube* dealiasingw = NULL;
//  
//  levelbase::rccube* advu = NULL;
//  levelbase::rccube* advv = NULL;
//  levelbase::rccube* advw = NULL;
//  levelbase::rccube* hu   = NULL;
//  levelbase::rccube* hv   = NULL;
//  levelbase::rccube* hw   = NULL;
//  levelbase::rccube* us   = NULL;
//  levelbase::rccube* vs   = NULL;
//  levelbase::rccube* ws   = NULL;
//  levelbase::rccube* lapu = NULL;
//  levelbase::rccube* lapv = NULL;
//  levelbase::rccube* lapw = NULL;
//  levelbase::rccube* rhsu = NULL;
//  levelbase::rccube* rhsv = NULL;
//  levelbase::rccube* rhsw = NULL;
//  levelbase::rccube* rhsp = NULL;
//  levelbase::rccube* divs = NULL;
  
  
  
  levelbase::rccube* uvwnonlinear = NULL;
  levelbase::rccube* advuvw = NULL;
  levelbase::rccube* huvw = NULL;
  levelbase::rccube* suvw = NULL;
  levelbase::rccube* lapuvw = NULL;
  levelbase::rccube* rhsuvw = NULL;
  levelbase::rccube* rhsp   = NULL;
  levelbase::rccube* divs   = NULL;
  
  
  levelbase::rccube* tmpu1   = NULL;
  levelbase::rccube* tmpv1   = NULL;
  levelbase::rccube* tmpw1   = NULL;
  
  levelbase::rccube* tmpu2   = NULL;
  levelbase::rccube* tmpv2   = NULL;
  levelbase::rccube* tmpw2   = NULL;
  
  levelbase::rccube* tmpp1   = NULL;
  levelbase::rccube* tmpp2   = NULL;
  
  
protected:
  
  string name_ = "FlatCarrierVariable::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_FLAT_CARRIER_VARIABLE_H_
