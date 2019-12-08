/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <iomanip>
#include "boundary.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
// Function : Constructor
//------------------------------------------------------------------------------
// Parameter: IN  - mpi       > mpi object
//            IN  - dcp      > domain decomposition
//------------------------------------------------------------------------------
Boundary::Boundary(const mpi_t& mpi,
                   const dcp_t& dcp,
                   Domain      * domain,
                   SpaceStepper* sstepper,
                   TimeStepper * tstepper,
                   Variable    * variable,
                   const int index) :
  levelbase::Base(mpi,dcp,NULL),
  domain_(domain),sstepper_(sstepper),tstepper_(tstepper),variable_(variable),index_(index) {
  PS_DEBUG_TRACE_ENTER(name_+"Boundary(mpi_t,dcp_t)")
    
  // Set the name of output file.
    ostringstream indexstring;
    indexstring << std::setw(2) << std::setfill('0') << index_;
    filename_ = string("phase_")+indexstring.str()+"_"+IO_FILE_CLASS_BOUNDARY;
    
  PS_DEBUG_TRACE_LEAVE(name_+"Boundary(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void Boundary::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << "Index " << index_;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool Boundary::SetFlagOpenTop(const dvec& zboundarytype) {
  PS_DEBUG_TRACE_ENTER(name_+"SetFlagOpenTop(dvec)")
  
  if (std::abs(zboundarytype(2)-0)<=COMPARE_NUMBER_SMALL &&
      std::abs(zboundarytype(3)-1)<=COMPARE_NUMBER_SMALL) {
    isopentop_ = true;
  } else {
    isopentop_ = false;
  }
  
  return true;
  PS_DEBUG_TRACE_LEAVE(name_+"SetFlagOpenTop(dvec)")
}
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
