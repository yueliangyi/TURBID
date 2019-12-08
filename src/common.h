/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_COMMON_H_
#define YLY_PS_LVBASE_COMMON_H_

#include <string>
#include "configure.h"
#include "armadillo"


namespace yueliangyi {
namespace pseudospectral {


//
typedef std::string string;
typedef std::fstream fstream;
//typedef std::stringstream stringstream;
typedef std::ostringstream ostringstream;
  
typedef arma::cx_double cx_double;
typedef arma::Col   <int            > ivec;
typedef arma::Col   <double         > dvec;
typedef arma::Col   <arma::uword    > uvec;
typedef arma::Mat   <int            > imat;
typedef arma::Mat   <arma::uword    > umat;
typedef arma::Mat   <double         > dmat;
typedef arma::Mat   <arma::cx_double> cmat;
typedef arma::Cube  <double         > dcube;
typedef arma::Cube  <arma::cx_double> ccube;


//
extern ostringstream _message_;


  
//#define PI                          (245850922.0/78256779.0)
#define COMPARE_NUMBER_LARGE        (10E+12)
#define COMPARE_NUMBER_SMALL        (10E-12)


#if defined(PS_DEBUG_TRACE)
#define PS_DEBUG_TRACE_ENTER(msg) { _message_ << "DEBUG:Enter > " << msg << std::endl; }
#define PS_DEBUG_TRACE_LEAVE(msg) { _message_ << "DEBUG:Leave > " << msg << std::endl; }
#else
#define PS_DEBUG_TRACE_ENTER(msg) {}
#define PS_DEBUG_TRACE_LEAVE(msg) {}
#endif
  
  

#define BOUNDARY_COEF_SCALE        (10E+6)
#define BOUNDARY_COEF_SCALE_UNEVEN (1.0)

  
extern string IO_FILE_CONFIGURATION;
extern string IO_FILE_LOG_COMPUTATION;

  
  
  
  

} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_COMMON_H_
