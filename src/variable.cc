/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "variable.h"
#include <iomanip>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
//std::vector<rccube*> Variable::u;
//std::vector<rccube*> Variable::dudx;
//std::vector<rccube*> Variable::dudy;
//std::vector<rccube*> Variable::dudz;
//
//std::vector<rccube*> Variable::v;
//std::vector<rccube*> Variable::dvdx;
//std::vector<rccube*> Variable::dvdy;
//std::vector<rccube*> Variable::dvdz;
//
//std::vector<rccube*> Variable::w;
//std::vector<rccube*> Variable::dwdx;
//std::vector<rccube*> Variable::dwdy;
//std::vector<rccube*> Variable::dwdz;
//
//std::vector<rccube*> Variable::p;
//std::vector<rccube*> Variable::dpdx;
//std::vector<rccube*> Variable::dpdy;
//std::vector<rccube*> Variable::dpdz;
  
  
  std::vector<levelbase::rccube*> Variable::uvw;
  std::vector<levelbase::rccube*> Variable::dudxyz;
  std::vector<levelbase::rccube*> Variable::dvdxyz;
  std::vector<levelbase::rccube*> Variable::dwdxyz;
  
  std::vector<levelbase::rccube*> Variable::p;
  std::vector<levelbase::rccube*> Variable::dpdxyz;
  
  std::vector<levelbase::rccube*> Variable::dudxyztmp;
  std::vector<levelbase::rccube*> Variable::dvdxyztmp;
  std::vector<levelbase::rccube*> Variable::dwdxyztmp;
  std::vector<levelbase::rccube*> Variable::ptmp;
  
  std::vector<levelbase::rccube*> Variable::botsstres;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Variable::Variable(const mpi_t& mpi,
                   const dcp_t& dcp,
                   Domain      * domain,
                   SpaceStepper* sstepper,
                   TimeStepper * tstepper,
                   const int index) :
  Base(mpi,dcp,NULL),
  domain_(domain),sstepper_(sstepper),tstepper_(tstepper),index_(index) {
  PS_DEBUG_TRACE_ENTER("Variable::Variable(mpi_t,dcp_t)")
    

#ifdef PS_DEBUG_CHECK
    if (index_ != int(uvw.size())) {
      name_ += "Variable(mpi_t,dcp_t) > Wrong Index of Variable!";
      throw std::runtime_error(name_);
    }
#endif
    
    
    uvw   .push_back(new rccube(mpi_,dcp_,3));
    dudxyz.push_back(new rccube(mpi_,dcp_,3));
    dvdxyz.push_back(new rccube(mpi_,dcp_,3));
    dwdxyz.push_back(new rccube(mpi_,dcp_,3));

    p     .push_back(new rccube(mpi_,dcp_,1));
    dpdxyz.push_back(new rccube(mpi_,dcp_,3));
    
    dudxyztmp.push_back(new rccube(mpi_,dcp_,3));
    dvdxyztmp.push_back(new rccube(mpi_,dcp_,3));
    dwdxyztmp.push_back(new rccube(mpi_,dcp_,3));
    ptmp.push_back(new rccube(mpi_,dcp_,1));
    
    botsstres.push_back(new rccube(mpi_,dcp_,1));

    ostringstream indexstring;
    indexstring << std::setw(2) << std::setfill('0') << index_;
    filename_ = string("phase_")+indexstring.str()+"_"+IO_FILE_CLASS_VARIABLE;
    
    
  PS_DEBUG_TRACE_LEAVE("Variable::Variable(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Variable::~Variable(void) {
  
  if (uvw   [index_] != NULL) { delete uvw   [index_]; uvw   [index_] = NULL; }
  if (dudxyz[index_] != NULL) { delete dudxyz[index_]; dudxyz[index_] = NULL; }
  if (dvdxyz[index_] != NULL) { delete dvdxyz[index_]; dvdxyz[index_] = NULL; }
  if (dwdxyz[index_] != NULL) { delete dwdxyz[index_]; dwdxyz[index_] = NULL; }
  
  if (p     [index_] != NULL) { delete p     [index_]; p     [index_] = NULL; }
  if (dpdxyz[index_] != NULL) { delete dpdxyz[index_]; dpdxyz[index_] = NULL; }
  
  if (dudxyztmp[index_] != NULL) { delete dudxyztmp[index_]; dudxyztmp[index_] = NULL; }
  if (dvdxyztmp[index_] != NULL) { delete dvdxyztmp[index_]; dvdxyztmp[index_] = NULL; }
  if (dwdxyztmp[index_] != NULL) { delete dwdxyztmp[index_]; dwdxyztmp[index_] = NULL; }
  if (ptmp[index_] != NULL) { delete ptmp[index_]; ptmp[index_] = NULL; }
  
  if (botsstres[index_] != NULL) { delete botsstres[index_]; botsstres[index_] = NULL; }

}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void Variable::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream << "Index " << index_;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
