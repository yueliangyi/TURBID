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
  
  
  
class FlatChannelSolver : public System {
  
public:
  
  FlatChannelSolver(int& argc, char**& argv);
  ~FlatChannelSolver(void);
  
  bool Configure(void);
  
  bool Execute(void);
  
  void Print(ostringstream& stream) const;
  void PrintLogHead(ostringstream& stream) const;
  void PrintLogBody(ostringstream& stream) const;
  
  
protected:
  
  
  string name_ = "FlatChannelSolver::";

private:
  
};
  
  
  
} // end namespace leveltwo
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVTWO_SOLVER_CHANNEL_FLAT_H_
