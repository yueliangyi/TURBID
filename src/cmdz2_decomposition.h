/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_CMDZ2_DECOMPOSITION_H_
#define YLY_PS_LVBASE_CMDZ2_DECOMPOSITION_H_

#include "base.h"
#include "cube_real_cplx.h"
#include "cube_real_decp.h"
#include "stepper_space.h"
#include "domain.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class CmDz2DecompositionCplx {
  
public:
  
  bool Construct(const SpaceStepper* sstepper,
                 const dvec& zboundarytype,
                 const dvec& stepcoefficient,
                 const int dimension = 1,
                 const bool includedxy2 = true);

  ~CmDz2DecompositionCplx(void);
  
  const cmat& GetFixForwardTop(void) const { return fixforwardtop_; }
  const cmat& GetFixForwardBottom(void) const { return fixforwardbottom_; }
  const cmat& GetFixBackwardTop(void) const { return fixbackwardtop_; }
  const cmat& GetFixBackwardBottom(void) const { return fixbackwardbottom_; }
  
  const cmat& GetEigenVectorNorm(void) const { return eigenvectornorm_; }
  const cmat& GetEigenVectorInverse(void) const { return eigenvectorinverse_; }
  const levelbase::rccube& GetGetEigenVectorReciprocal(const int id) const {
    return *eigenvaluereciprocal_[id];
  }
  
//  const dmat& GetCoefficientMatrix(void) const { return coefficientmatrix_; }
  
  bool Assemble(const SpaceStepper* sstepper, const dvec& stepcoefficient);
  
protected:
  
  int dimension_ = 1;
  bool includedxy2_ = true;
  
  dvec zboundarytype_;
  dvec stepcoefficient_;
  dmat coefficientmatrix_;
  cmat fixforwardtop_;
  cmat fixforwardbottom_;
  cmat fixbackwardtop_;
  cmat fixbackwardbottom_;
  cmat eigenvalue_;
  cmat eigenvectornorm_;
  cmat eigenvectorinverse_;
  std::vector<levelbase::rccube*> eigenvaluereciprocal_;
  
  dmat cmdx2_;
  dmat cmdy2_;

private:
  
}; typedef CmDz2DecompositionCplx dz2dc_t;
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_CMDZ2_DECOMPOSITION_H_
