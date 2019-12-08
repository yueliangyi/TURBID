/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <iomanip>
#include "phase.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  dvec Phase::featurelocal_ = arma::zeros<dvec>(50);
  dvec Phase::featureglobal_ = arma::zeros<dvec>(50);
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool Phase::Configuration::Load(IniFile& fileini) {
  PS_DEBUG_TRACE_ENTER("Phase::Configuration::Load(IniFile)")
  
  PS_DEBUG_TRACE_LEAVE("Phase::Configuration::Load(IniFile)")
  return true;
}
  

  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Phase::Phase(const mpi_t& mpi,
             const dcp_t& dcp,
             Domain      * domain,
             SpaceStepper* sstepper,
             TimeStepper * tstepper,
             const Base::cfg_t* cfg,
             const int index) :
  Base(mpi,dcp,cfg),
  domain_(domain),sstepper_(sstepper),tstepper_(tstepper),index_(index) {
  PS_DEBUG_TRACE_ENTER("Phase::Phase(mpi_t,dcp_t,cfg_t)")
    
    
//    featurelocal_ = arma::zeros<dvec>(50);
//    featureglobal_ = arma::zeros<dvec>(50);
    
    
    ostringstream indexstring;
    indexstring << std::setw(2) << std::setfill('0') << index_;
    filename_ = string("phase_")+indexstring.str()+"_"+IO_FILE_CLASS_PHASE;
    
    
  PS_DEBUG_TRACE_LEAVE("Phase::Phase(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Phase::~Phase(void) {
  domain_   = NULL;
  tstepper_ = NULL;
  sstepper_ = NULL;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void Phase::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("Phase::Print(ostringstream)")
  
  stream << "Index " << index_;
  stream << std::endl;
  
  variable_->Print(stream);
  boundary_->Print(stream);
  
  
  
  
  PS_DEBUG_TRACE_LEAVE("Phase::Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
