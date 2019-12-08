/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVTWO_SYSTEM_H_
#define YLY_PS_LVTWO_SYSTEM_H_

#include "phase.h"
#include "file_log.h"
#include "file_ini.h"

namespace yueliangyi {
namespace pseudospectral {
namespace leveltwo {
  
  
  
  
  // The system is based on MPI (Message Passing Interface) and OpenMP (Open
  // Multi-Processing). However, the feature of OpenMP is inherited from the basic
  // libraries (P3DFFT and Armadillo) directly. Thus, the major task here is to
  // processes all the prerequisites of the system. For different threads/thread
  // group (the group shares the same MPI configuration), specific MPIGroup
  // information should be created. All objects created for one thread (group)
  // will also share the same MPIGroup. Remember to create the MPIGroup before
  // creating other fundamental objects.
  
class System {
  
public:
  
  System(int& argc, char**& argv);
  virtual ~System(void) { }
  
  bool isConfigured(void) { return isconfigured_; }
  
  virtual bool Configure(void) = 0;
  
  virtual bool Execute(void) = 0;
  
  virtual void Print(ostringstream& stream) const;
  virtual void PrintLogHead(ostringstream& stream) const = 0;
  virtual void PrintLogBody(ostringstream& stream) const = 0;
  
protected:
  
  // Disable copy constructor with feature in C++11.
  System(const System& other) = delete;
  
  ivec size_;
  ivec mpigrid_;
  
  levelbase::mpi_t mpi_;
  levelbase::dcp_t dcp_;
  
  
  bool isconfigured_ = false;
  bool printlog_ = true;
  bool writelog_ = true;
  int  phasenum_ = 1;
  
  
  // Description of component type.
  string typedomain_  ;
  string typetstepper_;
  string typesstepper_;
  
  
  levelone::domain_t domain_;
  levelone::sstepper_t sstepper_;
  levelone::tstepper_t tstepper_;
  std::vector<levelone::phase_t> phase_;
  
  
  levelbase::LogFile filelog_;
  levelbase::IniFile fileini_;
  
  
  string name_ = "System::";

private:
  
};
  
  
  
} // end namespace leveltwo
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVTWO_SYSTEM_H_
