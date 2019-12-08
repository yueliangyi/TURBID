/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVTWO_SOLVER_CHANNEL_FLAT_H_
#define YLY_PS_LVTWO_SOLVER_CHANNEL_FLAT_H_

#include "system.h"

namespace yueliangyi {
namespace pseudospectral {
namespace leveltwo {
  
  
  
class FlatChannelPostprocessor : public System {
  
public:
  
  FlatChannelPostprocessor(int& argc, char**& argv);
  ~FlatChannelPostprocessor(void);
  
  bool Configure(void);
  
  bool Execute(void);
  
  void Print(ostringstream& stream) const;
  void PrintLogHead(ostringstream& stream) const;
  void PrintLogBody(ostringstream& stream) const;
  
  
protected:
  
  ivec work_index_range;
  
  
  string name_ = "FlatChannelPostprocessor::";

private:
  
};
  
  
  
} // end namespace leveltwo
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVTWO_SOLVER_CHANNEL_FLAT_H_
