/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_DOMAIN_FIXED_H_
#define YLY_PS_LVONE_DOMAIN_FIXED_H_

#include "domain.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class FixedDomain : public Domain {
  
public:
  
  FixedDomain(const levelbase::mpi_t& mpi,
              const levelbase::dcp_t& dcp,
              const levelbase::Base::cfg_t* cfg);
  
  
  bool Save(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  bool Load(const string subfolder = "",
            const ioctrl_t level = IOControlLevelOne);
  
  bool Configure(void) {
    Domain::Configure();
    return true;
  }
  
  bool Update(void) { Domain::Update(); return true; }
  void Print(ostringstream& stream) const;
  
protected:
  
  
  string name_ = "FixedDomain::";

private:
  
};
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_DOMAIN_FIXED_H_
