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
  
  
  
  
  
class CmDz2DecompositionReal {
  
public:
  
  bool Construct(SpaceStepper* sstepper,
                 const dvec& zboundarytype,
                 const dvec& stepcoefficient,
                 const int dimension = 1,
                 const bool includedxy2 = true);

  ~CmDz2DecompositionReal(void);
  
  const dmat& GetFixForwardTop(void) const { return fixforwardtop_; }
  const dmat& GetFixForwardBottom(void) const { return fixforwardbottom_; }
  const dmat& GetFixBackwardTop(void) const { return fixbackwardtop_; }
  const dmat& GetFixBackwardBottom(void) const { return fixbackwardbottom_; }
  
  const dmat& GetEigenVectorNorm(void) const { return eigenvectornorm_; }
  const dmat& GetEigenVectorInverse(void) const { return eigenvectorinverse_; }
  const levelbase::rdcube& GetGetEigenVectorReciprocal(const int id) const {
    return *eigenvaluereciprocal_[id];
  }
  
  const dmat& GetCoefficientMatrix(void) const { return coefficientmatrix_; }
  const levelbase::rdcube& GetJm3k6(void) const { return *jm3k6; }
  const levelbase::rdcube& GetJm3k45(void) const { return *jm3k45; }
  const dmat& GetCoefficientMatrixDz1(void) const { return coefficientmatrixdz1_; }
  
protected:
  
  int dimension_ = 1;
  bool includedxy2_ = true;
  
  dvec zboundarytype_;
  dvec stepcoefficient_;
  dmat coefficientmatrix_;
  dmat fixforwardtop_;
  dmat fixforwardbottom_;
  dmat fixbackwardtop_;
  dmat fixbackwardbottom_;
  dmat eigenvalue_;
  dmat eigenvectornorm_;
  dmat eigenvectorinverse_;
  std::vector<levelbase::rdcube*> eigenvaluereciprocal_;
  
  levelbase::rdcube* jm3k6 = NULL;
  levelbase::rdcube* jm3k45 = NULL;
  dmat coefficientmatrixdz1_;

  
  
private:
  
}; typedef CmDz2DecompositionReal dz2dr_t;
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_CMDZ2_DECOMPOSITION_H_
