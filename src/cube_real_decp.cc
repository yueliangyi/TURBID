/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/

#include <stdexcept>
#ifdef USE_MKL
# include <mkl_cblas.h>
#else
# include <cblas.h>
#endif
#include "cube_real_decp.h"
#include "p3dfft.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
RealDecompositionCube::RealDecompositionCube(const mpi_t& mpi,
                                             const dcp_t& dcp,
                                             const int dimension) :
  mpi_(mpi),dcp_(dcp),dimension_(dimension) {
  PS_DEBUG_TRACE_ENTER(name_+"RealDecompositionCube(mpi_t,dcp_t,rccube)")
    
    const int rank = mpi_.GetRank();
    const imat& dcpdimension = dcp_.GetDimension();
    
    cuberealzjoin_ = dcube(dcpdimension(10*4+6,rank),
                           dcpdimension(10*4+7,rank),
                           dcpdimension(10*4+8,rank)*dimension_,
                           arma::fill::zeros);
    
    cuberealxbuffer_ = dcube(dcpdimension(10*4+6,rank),
                             dcpdimension(10*4+7,rank),
                             dcpdimension(10*4+8,rank)*dimension_,
                             arma::fill::zeros);
    
    cuberealzbuffer_ = dcube(dcpdimension(10*4+6,rank),
                             dcpdimension(10*4+7,rank),
                             dcpdimension(10*4+8,rank)*dimension_,
                             arma::fill::zeros);
    
    dmatrealzjoin_ = dmat(&cuberealzjoin_(0,0,0),
                         dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                         dcpdimension(10*4+8,rank)*dimension_,
                         false,
                         false);
    
    
    cuberealz_.resize(dimension_);
    dmatrealz_.resize(dimension_);
    
    
    for (int id = 0; id < dimension_; ++id) {
      cuberealz_[id] = dcube(&cuberealzjoin_(0,0,dcpdimension(10*4+8,rank)*id),
                             dcpdimension(10*4+6,rank),
                             dcpdimension(10*4+7,rank),
                             dcpdimension(10*4+8,rank),
                             false,
                             false);
      dmatrealz_[id] = dmat(&cuberealzjoin_(0,0,dcpdimension(10*4+8,rank)*id),
                            dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                            dcpdimension(10*4+8,rank),
                            false,
                            false);
    }
    
    
    
    
    
    
    
  PS_DEBUG_TRACE_LEAVE(name_+"RealDecompositionCube(mpi_t,dcp_t,rccube)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Fill all data members with given value.
//            This function is created for tests mainly. Try to avoid the usage
//            of it in computation. Note that the complex member will have the
//            same value for the real and imaginary part.
//------------------------------------------------------------------------------
// Parameter: IN  - value     > required value for filling
//------------------------------------------------------------------------------
void RealDecompositionCube::Fill(const double value, const int dimension) {
  PS_DEBUG_TRACE_ENTER(name_+"Fill(double)")
  if (dimension <= 0) {
    GetCubeRealZJoin().fill(value);
  } else {
    GetCubeRealZ(dimension).fill(value);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"Fill(double)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool RealDecompositionCube::RealPencilX2Z(rccube& rccubeobj) {
  PS_DEBUG_TRACE_ENTER(name_+"RealPencilX2Z(rccube)")
  
  
  // check
#ifdef PS_DEBUG_CHECK
  if (dimension_ != rccubeobj.GetDimension()) {
    name_ += "RealPencilX2Z(rccube) > Dimension Mismatch!";
    throw std::runtime_error(name_);
  }
#endif
  
  
  ivec dstart(3);
  ivec dend(3);
  ivec dsize(3);
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    Cp3dfft_rtran_x2z(rccubeobj.GetCubeReal(iddim).memptr(),
                      GetCubeRealZ(iddim).memptr(),
                      cuberealxbuffer_.memptr(),
                      cuberealzbuffer_.memptr(),
                      dstart.memptr(),
                      dend.memptr(),
                      dsize.memptr());
  }
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"RealPencilX2Z(rccube)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool RealDecompositionCube::RealPencilZ2X(rccube& rccubeobj) {
  PS_DEBUG_TRACE_ENTER(name_+"RealPencilZ2X(rccube)")
  
  // check
#ifdef PS_DEBUG_CHECK
  if (dimension_ != rccubeobj.GetDimension()) {
    name_ += "RealPencilZ2X(rccube) > Dimension Mismatch!";
    throw std::runtime_error(name_);
  }
#endif
  
  ivec dstart(3);
  ivec dend(3);
  ivec dsize(3);
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    Cp3dfft_rtran_z2x(GetCubeRealZ(iddim).memptr(),
                      rccubeobj.GetCubeReal(iddim).memptr(),
                      cuberealzbuffer_.memptr(),
                      cuberealxbuffer_.memptr(),
                      dstart.memptr(),
                      dend.memptr(),
                      dsize.memptr());
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"RealPencilZ2X(rccube)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool RealDecompositionCube::InplaceTranspose(void) {
  PS_DEBUG_TRACE_ENTER(name_+"InplaceTranspose(void)")
  
  for (int id = 0; id < dimension_; ++id) {
    arma::inplace_trans(dmatrealz_[id]);
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"InplaceTranspose(void)")
  return true;
}
bool RealDecompositionCube::MatrixForward(void) {
  PS_DEBUG_TRACE_ENTER(name_+"MatrixForward(void)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  for (int id = 0; id < dimension_; ++id) {
    dmatrealz_[id] = dmat(&cuberealzjoin_(0,0,dcpdimension(10*4+8,rank)*id),
                          dcpdimension(10*4+8,rank),
                          dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                          false,
                          false);
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"MatrixForward(void)")
  return true;
}
bool RealDecompositionCube::MatrixBackward(void) {
  PS_DEBUG_TRACE_ENTER(name_+"MatrixBackward(void)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();

  for (int id = 0; id < dimension_; ++id) {
    dmatrealz_[id] = dmat(&cuberealzjoin_(0,0,dcpdimension(10*4+8,rank)*id),
                          dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                          dcpdimension(10*4+8,rank),
                          false,
                          false);
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"MatrixBackward(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void RealDecompositionCube::Multiply(const dmat& lhs1, RealDecompositionCube& result) {
  PS_DEBUG_TRACE_ENTER(name_+"Multiply(dmat,RealDecompositionCube)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  for (int id = 0; id < dimension_; ++id) {
    cblas_dgemm(CblasColMajor,
                CblasNoTrans,
                CblasTrans,
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                dcpdimension(10*4+8,rank),
                dcpdimension(10*4+8,rank),
                1.0,
                this->GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                lhs1.memptr(),
                dcpdimension(10*4+8,rank),
                0.0,
                result.GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank));
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"Multiply(dmat,RealDecompositionCube)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void RealDecompositionCube::Multiply(const dmat& lhs1, const int position) {
  PS_DEBUG_TRACE_ENTER(name_+"Multiply(dmat,int)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  dmat tmp = dmat(dcpdimension(10*4+6,rank),dcpdimension(10*4+7,rank));
  
  for (int id = 0; id < dimension_; ++id) {
    cblas_dgemm(CblasColMajor,
                CblasNoTrans,
                CblasTrans,
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                1,
                dcpdimension(10*4+8,rank),
                1.0,
                this->GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                lhs1.memptr(),
                1,
                0.0,
                tmp.memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank));
    this->GetCubeRealZ(id).slice(position) = tmp;
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"Multiply(dmat,int)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void RealDecompositionCube::Multiply(const dmat& lhs1,
                                     const RealDecompositionCube& lhs2,
                                     const dmat& lhs3,
                                     RealDecompositionCube& result) {
  PS_DEBUG_TRACE_ENTER(name_+"Multiply(dmat,RealDecompositionCube,dmat,RealDecompositionCube)")
  
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  for (int id = 0; id < dimension_; ++id) {
    
    cblas_dgemm(CblasColMajor,
                CblasNoTrans,
                CblasTrans,
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                dcpdimension(10*4+8,rank),
                dcpdimension(10*4+8,rank),
                1.0,
                this->GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                lhs3.memptr(),
                dcpdimension(10*4+8,rank),
                0.0,
                result.GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank));
    
    
    this->GetCubeRealZ(id) = lhs2.GetCubeRealZ(id)%result.GetCubeRealZ(id);
    
    
    cblas_dgemm(CblasColMajor,
                CblasNoTrans,
                CblasTrans,
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                dcpdimension(10*4+8,rank),
                dcpdimension(10*4+8,rank),
                1.0,
                this->GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank),
                lhs1.memptr(),
                dcpdimension(10*4+8,rank),
                0.0,
                result.GetCubeRealZ(id).memptr(),
                dcpdimension(10*4+6,rank)*dcpdimension(10*4+7,rank));
    
  }
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Multiply(dmat,RealDecompositionCube,dmat,RealDecompositionCube)")
}
  
  

  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi
