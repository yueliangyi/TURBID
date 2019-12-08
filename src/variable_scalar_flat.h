/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_FLAT_SCALAR_VARIABLE_H_
#define YLY_PS_LVONE_FLAT_SCALAR_VARIABLE_H_

#include "variable.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FlatScalarVariable : public Variable {
  
public:
  
  
  FlatScalarVariable(const levelbase::mpi_t& mpi,
                      const levelbase::dcp_t& ddbt,
                      Domain      * domain,
                      SpaceStepper* sstepper,
                      TimeStepper * tstepper,
                      const int index);
  ~FlatScalarVariable(void);
  
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
  
  levelbase::rccube* uvwnonlinear = NULL;
  levelbase::rccube* advphi = NULL;
  levelbase::rccube* hphi = NULL;
  levelbase::rccube* lapphi = NULL;
  levelbase::rccube* rhsphi   = NULL;
  
protected:
  
  string name_ = "FlatScalarVariable::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_FLAT_SCALAR_VARIABLE_H_
