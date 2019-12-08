/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_FLAT_SCALAR_BOUNDARY_H_
#define YLY_PS_LVONE_FLAT_SCALAR_BOUNDARY_H_

#include "boundary.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FlatScalarBoundary : public Boundary {
  
public:
  
  FlatScalarBoundary(const levelbase::mpi_t& mpi,
                      const levelbase::dcp_t& dcp,
                      Domain      * domain,
                      SpaceStepper* sstepper,
                      TimeStepper * tstepper,
                      Variable    * variable,
                      const int index);
  
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  bool Update(const int stage);
  bool Apply (const int stage);
  
  bool Configure(void);
  void Print(ostringstream& stream) const;
  
  // Variables for collecting previous information.
//  cmat dpdxyzm1top; cmat dpdxyzm1bottom;
//  cmat dpdxyzm2top; cmat dpdxyzm2bottom;
  
  
protected:
  
  string name_ = "FlatScalarBoundary::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_FLAT_SCALAR_BOUNDARY_H_
