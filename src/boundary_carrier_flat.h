/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_FLAT_CARRIER_BOUNDARY_H_
#define YLY_PS_LVONE_FLAT_CARRIER_BOUNDARY_H_

#include "boundary.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FlatCarrierBoundary : public Boundary {
  
public:
  
  FlatCarrierBoundary(const levelbase::mpi_t& mpi,
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
  arma::Row<cx_double> dp1dxyzm1top; arma::Row<cx_double> dp1dxyzm1bottom;
  arma::Row<cx_double> dp2dxyzm1top; arma::Row<cx_double> dp2dxyzm1bottom;
  arma::Row<cx_double> dp1dxyzm2top; arma::Row<cx_double> dp1dxyzm2bottom;
  arma::Row<cx_double> dp2dxyzm2top; arma::Row<cx_double> dp2dxyzm2bottom;
  
  double olddt_ = 0.0001;

protected:
  
  string name_ = "FlatCarrierBoundary::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_FLAT_CARRIER_BOUNDARY_H_
