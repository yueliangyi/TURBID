/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_CUBE_REAL_CPLX_H_
#define YLY_PS_LVBASE_CUBE_REAL_CPLX_H_

#include <vector>
#include "base.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
class RealCplxCube {
  
public:
  
  // Overload the file type in Armadillo.
  // Refer to the manual of Armadillo for description.
  enum FileType{
    kRawBinary  = arma::raw_binary ,
    kRawAscii   = arma::raw_ascii  ,
    kArmaBinary = arma::arma_binary,
    kArmaAscii  = arma::arma_ascii
  }; typedef FileType ftype_t;
  
  RealCplxCube(const mpi_t& mpi, const dcp_t& dcp, const int dimension = 1, const bool creal = true);
  
  const mpi_t & GetMpi  (void) const { return mpi_  ; }
  const dcp_t& GetDcp (void) const { return dcp_ ; }
  const dcube & GetCubeRealJoin(void) const { return cuberealjoin_; }
  const ccube & GetCubeCplxJoin(void) const { return cubecplxjoin_; }
  const cmat  & GetCMatCplxJoin(void) const { return cmatcplxjoin_; }
  
  const dcube & GetCubeReal(const int dimension = 0) const { return cubereal_[dimension]; }
  const ccube & GetCubeCplx(const int dimension = 0) const { return cubecplx_[dimension]; }
  const cmat  & GetCMatCplx(const int dimension = 0) const { return cmatcplx_[dimension]; }
  
  dcube& GetCubeRealJoin(void) { return cuberealjoin_; }
  ccube& GetCubeCplxJoin(void) { return cubecplxjoin_; }
  cmat & GetCMatCplxJoin(void) { return cmatcplxjoin_; }
  
  dcube& GetCubeReal(const int dimension = 0) { return cubereal_[dimension]; }
  ccube& GetCubeCplx(const int dimension = 0) { return cubecplx_[dimension]; }
  cmat & GetCMatCplx(const int dimension = 0) { return cmatcplx_[dimension]; }
  
  int GetDimension(void) const { return dimension_; }
  
  // Fill data members with given value.
  // Used for tests usually.
  void Fill(const double value = 0, const int dimension = 0);
  
  RealCplxCube& operator=(dcube& rhs);
  RealCplxCube& operator=(ccube& rhs);
  
  dcube& Export(dcube& target);
  dcube  ExportReal(void);
  ccube& Export(ccube& target);
  ccube  ExportCplx(void);

  // The function Multiply works only for z-pencil domain.
  void Multiply(const cmat& lhs1, RealCplxCube& result);
  void Multiply(const cmat& lhs1, const int position);
  void Multiply(const cmat& lhs1,
                const RealCplxCube& lhs2,
                const cmat& lhs3,
                RealCplxCube& result);
  

  bool DftForward(void);
  bool DftBackward(void);
  bool Dealiasing(void);
  bool Dealiasing2(void);

  bool SaveReal(const string& name, const ftype_t type = kArmaBinary);
//  bool SaveReal(const string& name, const ftype_t type = kArmaBinary, const int dimension = 0);
  bool SaveCplx(const string& name, const ftype_t type = kArmaBinary);
  bool LoadReal(const string& name, const ftype_t type = kArmaBinary);
  bool LoadCplx(const string& name, const ftype_t type = kArmaBinary);
  
  bool SaveReal(fstream& stream, const ftype_t type = kArmaBinary);
  bool SaveCplx(fstream& stream, const ftype_t type = kArmaBinary);
  bool LoadReal(fstream& stream, const ftype_t type = kArmaBinary);
  bool LoadCplx(fstream& stream, const ftype_t type = kArmaBinary);
  
protected:

private:
  
  const mpi_t& mpi_;
  const dcp_t& dcp_;
  
  const int dimension_ = 1;
  
  dcube cuberealjoin_;
  ccube cubecplxjoin_;
  cmat  cmatcplxjoin_ ;
  
  std::vector<dcube> cubereal_;
  std::vector<ccube> cubecplx_;
  std::vector<cmat > cmatcplx_;
  
  
  // Scale for plane DFT.
  double dftscale_ = 0;
  
  string name_ = "RealCplxCube::";
  
}; typedef RealCplxCube rccube;
  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_CUBE_REAL_CPLX_H_
