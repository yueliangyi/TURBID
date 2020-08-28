/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_CUBE_REAL_DECP_H_
#define YLY_PS_LVBASE_CUBE_REAL_DECP_H_

#include "cube_real_cplx.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
class RealDecompositionCube {
  
public:
  
  
  RealDecompositionCube(const mpi_t& mpi, const dcp_t& dcp, const int dimension = 1);
  
  
  const mpi_t& GetMpi(void) const { return mpi_; }
  const dcp_t& GetDcp(void) const { return dcp_; }
  
  const dcube& GetCubeRealZJoin(void) const { return cuberealzjoin_; }
  const dcube& GetCubeRealZ(const int dimension = 0) const { return cuberealz_[dimension]; }
  dcube& GetCubeRealZJoin(void) { return cuberealzjoin_; }
  dcube& GetCubeRealZ(const int dimension = 0) { return cuberealz_[dimension]; }
  
  int GetDimension(void) const { return dimension_; }
  
  dmat& GetDmatReal(const int dimension = 0) { return dmatrealz_[dimension]; }
  
  
  
  // Fill data members with given value.
  // Used for tests usually.
  void Fill(const double value = 0, const int dimension = 0);
  
  
  bool RealPencilX2Z(rccube& rccubeobj);
  bool RealPencilZ2X(rccube& rccubeobj);
  
  
  
  bool InplaceTranspose(void);
  bool MatrixForward(void);
  bool MatrixBackward(void);
  
  
  
  // The function Multiply works only for z-pencil domain.
  void Multiply(const dmat& lhs1, RealDecompositionCube& result);
  void Multiply(const dmat& lhs1, const int position);
  void Multiply(const dmat& lhs1,
                const RealDecompositionCube& lhs2,
                const dmat& lhs3,
                RealDecompositionCube& result);
  
  
protected:
  
private:
  
  const mpi_t& mpi_;
  const dcp_t& dcp_;
  
  std::vector<dcube> cuberealz_;
  std::vector<dmat > dmatrealz_;
  
  dcube cuberealzjoin_;
  dcube cuberealxbuffer_;
  dcube cuberealzbuffer_;
  
  dmat  dmatrealzjoin_ ;
  
  const int dimension_ = 1;
  
  string name_ = "RealDecompositionCube::";
  
  
}; typedef RealDecompositionCube rdcube;
  


} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_CUBE_REAL_DECP_H_
