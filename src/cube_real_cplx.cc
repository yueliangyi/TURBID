/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/

#include <stdexcept>
#include "p3dfft.h"
#include "cube_real_cplx.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
//------------------------------------------------------------------------------
// Function : Constructor
//            Initialize/Allocate the data members. Note that the memories are
//            managed by Armadillo, so there is no need of special treatment.
//            Since the out-place DFT will applied, separate memories are
//            allocated, corresponding to the x and z-pencil.
//------------------------------------------------------------------------------
// Parameter: IN  - mpi       > common mpi object
//            IN  - ddbt      > data distribution object
//------------------------------------------------------------------------------
RealCplxCube::RealCplxCube(const mpi_t& mpi, const dcp_t& dcp, const int dimension, const bool creal) :
  mpi_(mpi),dcp_(dcp),dimension_(dimension) {
  PS_DEBUG_TRACE_ENTER(name_+"RealCplxCube(mpi_t,dcp_t)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  dftscale_ = dcpdimension(10*1+3,rank);
    
  
    if (creal) {
    cuberealjoin_ = dcube(dcpdimension(10*2+6,rank),
                          dcpdimension(10*2+7,rank),
                          dcpdimension(10*2+8,rank)*dimension_,
                          arma::fill::zeros);
    }
    
    cubecplxjoin_ = ccube(dcpdimension(10*3+6,rank),
                          dcpdimension(10*3+7,rank),
                          dcpdimension(10*3+8,rank)*dimension_,
                          arma::fill::zeros);
    
    // Expand the cubez_ into a matrix for multiplication.
    cmatcplxjoin_ = cmat(&cubecplxjoin_(0,0,0),
                         dcpdimension(10*3+6,rank),
                         dcpdimension(10*3+7,rank)*dcpdimension(10*3+8,rank)*dimension_,
                         false,
                         false);
    
    
    

    if (creal) {
    cubereal_.resize(dimension_);
    }
    cubecplx_.resize(dimension_);
    cmatcplx_.resize(dimension_);
    
    
    //
    for (int id = 0; id < dimension_; ++id) {
      if (creal) {
      cubereal_[id] = dcube(&cuberealjoin_(0,0,dcpdimension(10*2+8,rank)*id),
                            dcpdimension(10*2+6,rank),
                            dcpdimension(10*2+7,rank),
                            dcpdimension(10*2+8,rank),
                            false,
                            false);
      }
      cubecplx_[id] = ccube(&cubecplxjoin_(0,0,dcpdimension(10*3+8,rank)*id),
                            dcpdimension(10*3+6,rank),
                            dcpdimension(10*3+7,rank),
                            dcpdimension(10*3+8,rank),
                            false,
                            false);
      cmatcplx_[id] = cmat(&cubecplxjoin_(0,0,dcpdimension(10*3+8,rank)*id),
                           dcpdimension(10*3+6,rank),
                           dcpdimension(10*3+7,rank)*dcpdimension(10*3+8,rank),
                           false,
                           false);
    }
    
    
    
  
  PS_DEBUG_TRACE_LEAVE(name_+"RealCplxCube(mpi_t,dcp_t)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Fill all data members with given value.
//            This function is created for tests mainly. Try to avoid the usage
//            of it in computation. Note that the complex member will have the
//            same value for the real and imaginary part.
//------------------------------------------------------------------------------
// Parameter: IN  - value     > required value for filling
//------------------------------------------------------------------------------
void RealCplxCube::Fill(const double value, const int dimension) {
  PS_DEBUG_TRACE_ENTER(name_+"Fill(double)")
  if (dimension <= 0) {
    GetCubeRealJoin().fill(value);
    GetCubeCplxJoin().fill(cx_double(value,value));
  } else {
    GetCubeReal(dimension).fill(value);
    GetCubeCplx(dimension).fill(cx_double(value,value));
  }
  PS_DEBUG_TRACE_LEAVE(name_+"Fill(double)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Fill the x-pencil data members in all MPI processes with given
//            local cube (double value). Usually, the local cube is stored in
//            main MPI with rank=0.
//------------------------------------------------------------------------------
// Parameter: IN  - rhs       > target RHS
//------------------------------------------------------------------------------
RealCplxCube& RealCplxCube::operator=(dcube& rhs) {
  PS_DEBUG_TRACE_ENTER(name_+"operator=(dcube)")
  
  const int rank = mpi_.GetRank();
  int errflag = MPI_SUCCESS;
  const imat& dcpdimension = dcp_.GetDimension();
  
#ifdef PS_DEBUG_CHECK
  if (mpi_.isMain() &&
      dcpdimension(10*1+0,rank) != (int)rhs.n_rows &&
      dcpdimension(10*1+1,rank) != (int)rhs.n_cols &&
      dcpdimension(10*1+2,rank)*dimension_ != (int)rhs.n_slices) {
    name_ += "operator=(dcube) > Size of RHS is not Suitable!";
    throw std::runtime_error(name_);
  }
#endif
  
  // Reorder the data into splitted boxes in x-pencil.
  dcube rhstemp;
  
  if (mpi_.isMain()) {
    size_t index = 0;
    rhstemp.set_size(size(rhs));
    
    for (int iddim = 0; iddim < dimension_; ++iddim) {
      
      for (int idmpi = 0; idmpi < mpi_.GetSize(); ++idmpi) {
        
//        if (dcpdimension(10*2+9,idmpi) <= 0) {
//          continue;
//        }
//        
//        std::cout
//        << dcpdimension(10*2+0,idmpi)-1 << "-"
//        << dcpdimension(10*2+1,idmpi)-1 << "-"
//        << dcpdimension(10*2+2,idmpi)-1+dcpdimension(10*1+2,idmpi)*iddim << " to "
//        << dcpdimension(10*2+3,idmpi)-1 << "-"
//        << dcpdimension(10*2+4,idmpi)-1 << "-"
//        << dcpdimension(10*2+5,idmpi)-1+dcpdimension(10*1+2,idmpi)*iddim
//        << std::endl;
        
        dcube dcubempi = rhs.subcube(dcpdimension(10*2+0,idmpi)-1,
                                     dcpdimension(10*2+1,idmpi)-1,
                                     dcpdimension(10*2+2,idmpi)-1+dcpdimension(10*1+2,idmpi)*iddim,
                                     dcpdimension(10*2+3,idmpi)-1,
                                     dcpdimension(10*2+4,idmpi)-1,
                                     dcpdimension(10*2+5,idmpi)-1+dcpdimension(10*1+2,idmpi)*iddim);
        for (size_t id = 0; id < dcubempi.n_elem; ++id) {
          rhstemp(index) = dcubempi(id); ++index;
        }
      }  // end of for (idmpi)
      
    }
    
  }
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    errflag += MPI_Scatterv(rhstemp.memptr()+dcpdimension(10*1+4,rank)*iddim,
                            dcp_.GetMpiBufferCountsReal().memptr(),
                            dcp_.GetMpiBufferDisplsReal().memptr(),
                            MPI_DOUBLE,
                            cuberealjoin_.memptr()+dcpdimension(10*2+9,rank)*iddim,
                            dcp_.GetMpiBufferCountsReal()(rank),
                            MPI_DOUBLE,
                            0,
                            MPI_COMM_WORLD);
  }

  
  if (errflag != MPI_SUCCESS) {
    name_ += "operator=(dcube) > Failed to Scatter Data to MPI Members!";
    throw std::runtime_error(name_);
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"operator=(dcube)")
  return *this;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Fill the z-pencil data members in all MPI processes with given
//            local cube (complex value). Usually, the local cube is stored in
//            main MPI with rank=0.
//------------------------------------------------------------------------------
// Parameter: IN  - rhs       > target RHS
//------------------------------------------------------------------------------
RealCplxCube& RealCplxCube::operator=(ccube& rhs) {
  PS_DEBUG_TRACE_ENTER(name_+"operator=(ccube)")
  
  const int rank = mpi_.GetRank();
  int errflag = MPI_SUCCESS;
  const imat& dcpdimension = dcp_.GetDimension();
  
#ifdef PS_DEBUG_CHECK
  if (mpi_.isMain() &&
      dcpdimension(10*1+5,rank) != (int)rhs.n_rows &&
      dcpdimension(10*1+6,rank) != (int)rhs.n_cols &&
      dcpdimension(10*1+7,rank)*dimension_ != (int)rhs.n_slices) {
    name_ += "operator=(ccube) > Size of RHS is not Suitable!";
    throw std::runtime_error(name_);
  }
#endif
  
  // Reorder the data into splitted boxes in z-pencil.
  ccube rhstemp;
  if (mpi_.isMain()) {
    size_t index = 0;
    rhstemp.set_size(size(rhs));
    
    for (int iddim = 0; iddim < dimension_; ++iddim) {
      
      for (int idmpi = 0; idmpi < mpi_.GetSize(); ++idmpi) {
        ccube dcubempi = rhs.subcube(dcpdimension(10*3+0,idmpi)-1,
                                     dcpdimension(10*3+1,idmpi)-1,
                                     dcpdimension(10*3+2,idmpi)-1+dcpdimension(10*1+7,rank)*iddim,
                                     dcpdimension(10*3+3,idmpi)-1,
                                     dcpdimension(10*3+4,idmpi)-1,
                                     dcpdimension(10*3+5,idmpi)-1+dcpdimension(10*1+7,rank)*iddim);
        for (size_t id = 0; id < dcubempi.n_elem; ++id) {
          rhstemp(index) = dcubempi(id); ++index;
        }
      } // end of for (idmpi)
      
    }
    
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    errflag += MPI_Scatterv(rhstemp.memptr()+dcpdimension(10*1+9,rank)*iddim,
                            dcp_.GetMpiBufferCountsCplx().memptr(),
                            dcp_.GetMpiBufferDisplsCplx().memptr(),
                            MPI_C_DOUBLE_COMPLEX,
                            cubecplxjoin_.memptr()+dcpdimension(10*3+9,rank)*iddim,
                            dcp_.GetMpiBufferCountsCplx()(rank),
                            MPI_C_DOUBLE_COMPLEX,
                            0,
                            MPI_COMM_WORLD);
  }
  
  if (errflag != MPI_SUCCESS) {
    name_ += "operator=(ccube) > Failed to Scatter Data to MPI Members!";
    throw std::runtime_error(name_);
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"operator=(ccube)")
  return *this;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Collect the x-pencil data from all MPI processes and export it to
//            the target output. Note that only the main MPI will have the exact
//            output. Others will have an empty one.
//------------------------------------------------------------------------------
// Parameter: IN  - target    > target cube
//------------------------------------------------------------------------------
dcube& RealCplxCube::Export(dcube& target) {
  PS_DEBUG_TRACE_ENTER(name_+"Export(dcube)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  int errflag = MPI_SUCCESS;
  dcube targettemp;
  
  if (mpi_.isMain()) {
    target.resize(dcpdimension(10*1+0,rank),
                  dcpdimension(10*1+1,rank),
                  dcpdimension(10*1+2,rank)*dimension_);
    targettemp.resize(size(target));
  } else {
    target.reset();
  }
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    errflag += MPI_Gatherv(cuberealjoin_.memptr()+dcpdimension(10*2+9,rank)*iddim,
                           dcp_.GetMpiBufferCountsReal()(rank),
                           MPI_DOUBLE,
                           targettemp.memptr()+dcpdimension(10*1+4,rank)*iddim,
                           dcp_.GetMpiBufferCountsReal().memptr(),
                           dcp_.GetMpiBufferDisplsReal().memptr(),
                           MPI_DOUBLE,
                           0,
                           MPI_COMM_WORLD);
  }
  
  if (errflag != MPI_SUCCESS) {
    name_ += "Export(dcube) > Failed to Gather Data from MPI Member!";
    throw std::runtime_error(name_);
  }
  
  
  if (mpi_.isMain()) {
    
    for (int iddim = 0; iddim < dimension_; ++iddim) {
      
      for (int idmpi = 0; idmpi < mpi_.GetSize(); ++idmpi) {
        dcube dcubempi = dcube(targettemp.memptr()+dcp_.GetMpiBufferDisplsReal()(idmpi)+dcpdimension(10*1+4,rank)*iddim,
                               dcpdimension(10*2+6,idmpi),
                               dcpdimension(10*2+7,idmpi),
                               dcpdimension(10*2+8,idmpi),
                               false);
        target.subcube(dcpdimension(10*2+0,idmpi)-1,
                       dcpdimension(10*2+1,idmpi)-1,
                       dcpdimension(10*2+2,idmpi)-1+dcpdimension(10*1+2,idmpi)*iddim,
                       dcpdimension(10*2+3,idmpi)-1,
                       dcpdimension(10*2+4,idmpi)-1,
                       dcpdimension(10*2+5,idmpi)-1+dcpdimension(10*1+2,idmpi)*iddim) = dcubempi;
      }
      
      
    }
    
  } // end of if
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Export(dcube)")
  return target;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Collect the x-pencil data from all MPI processes and return it as
//            the result. Note that only the main MPI will have the exact output.
//            Others will have an empty one.
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
dcube RealCplxCube::ExportReal(void) {
  PS_DEBUG_TRACE_ENTER(name_+"ExportReal(void)")
  dcube cubereal;
  this->Export(cubereal);
  PS_DEBUG_TRACE_LEAVE(name_+"ExportReal(void)")
  return cubereal;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Collect the z-pencil data from all MPI processes and export it to
//            the target output. Note that only the main MPI will have the exact
//            output. Others will have an empty one.
//------------------------------------------------------------------------------
// Parameter: IN  - target    > target cube
//------------------------------------------------------------------------------
ccube& RealCplxCube::Export(ccube& target) {
  PS_DEBUG_TRACE_ENTER(name_+"Export(ccube)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  int errflag = MPI_SUCCESS;
  ccube targettemp;
  
  if (mpi_.isMain()) {
    target.set_size(dcpdimension(10*1+5,rank),
                    dcpdimension(10*1+6,rank),
                    dcpdimension(10*1+7,rank)*dimension_);
    targettemp.resize(size(target));
  } else {
    target.reset();
  }
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    errflag += MPI_Gatherv(cubecplxjoin_.memptr()+dcpdimension(10*3+9,rank)*iddim,
                           dcp_.GetMpiBufferCountsCplx()(rank),
                           MPI_C_DOUBLE_COMPLEX,
                           targettemp.memptr()+dcpdimension(10*1+9,rank)*iddim,
                           dcp_.GetMpiBufferCountsCplx().memptr(),
                           dcp_.GetMpiBufferDisplsCplx().memptr(),
                           MPI_C_DOUBLE_COMPLEX,
                           0,
                           MPI_COMM_WORLD);
  }
  
  if (errflag != MPI_SUCCESS) {
    name_ += "Export(ccube) > Failed to Gather Data from MPI Member!";
    throw std::runtime_error(name_);
  }
  
  if (mpi_.isMain()) {
    
    for (int iddim = 0; iddim < dimension_; ++iddim) {
      for (int idmpi = 0; idmpi < mpi_.GetSize(); ++idmpi) {
        ccube dcubempi = ccube(targettemp.memptr()+dcp_.GetMpiBufferDisplsCplx()(idmpi)+dcpdimension(10*1+9,rank)*iddim,
                               dcpdimension(10*3+6,idmpi),
                               dcpdimension(10*3+7,idmpi),
                               dcpdimension(10*3+8,idmpi),
                               false);
        target.subcube(dcpdimension(10*3+0,idmpi)-1,
                       dcpdimension(10*3+1,idmpi)-1,
                       dcpdimension(10*3+2,idmpi)-1+dcpdimension(10*1+7,idmpi)*iddim,
                       dcpdimension(10*3+3,idmpi)-1,
                       dcpdimension(10*3+4,idmpi)-1,
                       dcpdimension(10*3+5,idmpi)-1+dcpdimension(10*1+7,idmpi)*iddim) = dcubempi;
        
      }
    }
    
  } // end of if
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"Export(ccube)")
  return target;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Collect the z-pencil data from all MPI processes and return it as
//            the result. Note that only the main MPI will have the exact output.
//            Others will have an empty one.
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
ccube RealCplxCube::ExportCplx(void) {
  PS_DEBUG_TRACE_ENTER(name_+"ExportCplx(void)")
  ccube cubecplx;
  this->Export(cubecplx);
  PS_DEBUG_TRACE_LEAVE(name_+"ExportCplx(void)")
  return cubecplx;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Apply matrix multiplication to z-pencil data.
//            Result after multiplication will not be stored back.
//------------------------------------------------------------------------------
// Parameter: IN  - lhs       > coefficient matrix
//            IN  - position  > target position
//------------------------------------------------------------------------------
void RealCplxCube::Multiply(const cmat& lhs, RealCplxCube& result) {
  PS_DEBUG_TRACE_ENTER(name_+"Multiply(cmat,RealCplxCube)")
#ifdef PS_DEBUG_CHECK
  if (this == &result) {
    name_ += "Multiply(cmat,RealCplxCube) > Can not Do In-place Multiplication!";
    throw std::runtime_error(name_);
  }
#endif
  result.GetCMatCplxJoin() = lhs*this->GetCMatCplxJoin();
  PS_DEBUG_TRACE_LEAVE(name_+"Multiply(cmat,RealCplxCube)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Apply matrix multiplication to z-pencil data.
//            With input position, the input coefficient matrix should be
//            row-like and the result will be stored back the specified slice.
//------------------------------------------------------------------------------
// Parameter: IN  - lhs       > coefficient matrix
//            IN  - position  > target position
//------------------------------------------------------------------------------
void RealCplxCube::Multiply(const cmat& lhs, const int position) {
  PS_DEBUG_TRACE_ENTER(name_+"Multiply(cmat,int)")
  cmat tmp = lhs*this->GetCMatCplxJoin();
  this->GetCMatCplxJoin().row(position) = tmp;
  PS_DEBUG_TRACE_LEAVE(name_+"Multiply(cmat,int)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Apply matrix multiplication to z-pencil data.
//            Extend the matrix multiplication in a form of
//                  result = lhs1*(lh2%(lh3*this))
//            Usually, this function is used for the solution of linear
//            equations based on the direct method.
//------------------------------------------------------------------------------
// Parameter: IN  - lhs1      > coefficient matrix
//            IN  - lhs2      > target position
//            IN  - lhs3      > coefficient matrix
//------------------------------------------------------------------------------
void RealCplxCube::Multiply(const cmat& lhs1, const RealCplxCube& lhs2, const cmat& lhs3, RealCplxCube& result) {
  PS_DEBUG_TRACE_ENTER(name_+"Multiply(cmat,RealCplxCube,cmat,RealCplxCube)")
#ifdef PS_DEBUG_CHECK
  if (this == &result) {
    name_ += "Multiply(cmat,RealCplxCube,cmat,RealCplxCube) > Can not Do In-place Multiplication!";
    throw std::runtime_error(name_);
  }
#endif
  result.GetCMatCplxJoin() = lhs1*(lhs2.GetCMatCplxJoin()%(lhs3*this->GetCMatCplxJoin()));
  PS_DEBUG_TRACE_LEAVE(name_+"Multiply(cmat,RealCplxCube,cmat,RealCplxCube)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Apply the discrete fourier transform on plane
//            The out-place strategy is applied, meaning that the original
//            x-pencil data will not be destroyed. Note that no action is
//            applied in z-direction.
//------------------------------------------------------------------------------
// Parameter: IN  - flag      > direction of DFT
//------------------------------------------------------------------------------
bool RealCplxCube::DftForward(void) {
  PS_DEBUG_TRACE_ENTER(name_+"DftForward(void)")
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  unsigned char forward[]="ffn";
  Cp3dfft_ftran_r2c_many(GetCubeRealJoin().memptr(),
                         dcpdimension(10*2+9,rank),
                         (double*)GetCubeCplxJoin().memptr(),
                         dcpdimension(10*3+9,rank),
                         dimension_,
                         forward);
  GetCubeCplxJoin() *= 1/dftscale_;
  PS_DEBUG_TRACE_LEAVE(name_+"DftForward(void)")
  return true;
}
  
  
bool RealCplxCube::DftBackward(void) {
  PS_DEBUG_TRACE_ENTER(name_+"DftBackward(void)")
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  unsigned char backward[]="nff";
  Cp3dfft_btran_c2r_many((double*)GetCubeCplxJoin().memptr(),
                         dcpdimension(10*3+9,rank),
                         GetCubeRealJoin().memptr(),
                         dcpdimension(10*2+9,rank),
                         dimension_,
                         backward);
  PS_DEBUG_TRACE_LEAVE(name_+"DftBackward(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Apply the 2/3 rule dealiasing technique to z-pencil data
//            The 2/3 dealiasing rule can be regarded as a cut-off function,
//            removing the high-order wavenumber components. Unfortunately,
//            compared to the 3/2 rule, this kind of technique may require
//            more grid points on plane.
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
bool RealCplxCube::Dealiasing(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Dealiasing(void)")
  
 /* 
  // note
  // +0 -0 not very good
  // +1 -1 maybe
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
    
    double cutx = 0;
    for (int id = 0; id < dcpdimension(10*3+8,rank); ++id) {
      cutx = double(dcpdimension(10*3+2,rank)+id-1+0)/double(dcpdimension(10*1+7,rank)-1);
//      cutx = double(dcpdimension(10*3+2,rank)+id-1+2)/double(dcpdimension(10*1+7,rank)-1);
      if (cutx>2.0/3.0) {
        GetCubeCplxJoin().slice(id+dcpdimension(10*3+8,rank)*iddim).fill(cx_double(0,0));
      }
    }
    
    double cuty = 0;
    for (int id = 0; id < dcpdimension(10*3+7,rank); ++id) {
      
      
      cuty = abs(double(dcpdimension(10*3+1,rank)+id-1-0-dcpdimension(10*1+6,rank)/2))/double(dcpdimension(10*1+6,rank)/2);
      if (cuty < 1.0/3.0) {
//      cuty = double(dcpdimension(10*3+1,rank)+id-1+1)/double(dcpdimension(10*1+6,rank)-1);
//      if (cuty>1.0/3.0 && cuty<2.0/3.0) {
        GetCubeCplxJoin().subcube(0,
                              id,
                              0,
                              dcpdimension(10*3+6,rank)-1,
                              id,
                              dcpdimension(10*3+8,rank)*dimension_-1).fill(cx_double(0,0));
      }
    }
    
  }
*/  
  PS_DEBUG_TRACE_LEAVE(name_+"Dealiasing(void)")
  return true;
}
//------------------------------------------------------------------------------
// Function : Apply the 2/3 rule dealiasing technique to z-pencil data
//            The 2/3 dealiasing rule can be regarded as a cut-off function,
//            removing the high-order wavenumber components. Unfortunately,
//            compared to the 3/2 rule, this kind of technique may require
//            more grid points on plane.
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
bool RealCplxCube::Dealiasing2(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Dealiasing(void)")
  

  // note
  // +0 -0 not very good
  // +1 -1 maybe
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  for (int iddim = 0; iddim < dimension_; ++iddim) {
  
    double cutx = 0;
    for (int id = 0; id < dcpdimension(10*3+8,rank); ++id) {
      cutx = double(dcpdimension(10*3+2,rank)+id-1+0)/double(dcpdimension(10*1+7,rank)-1);
//      cutx = double(dcpdimension(10*3+2,rank)+id-1+2)/double(dcpdimension(10*1+7,rank)-1);
      if (cutx>2.0/3.0) {
        GetCubeCplxJoin().slice(id+dcpdimension(10*3+8,rank)*iddim).fill(cx_double(0,0));
      }
    }
  
    double cuty = 0;
    for (int id = 0; id < dcpdimension(10*3+7,rank); ++id) {
  
  
      cuty = abs(double(dcpdimension(10*3+1,rank)+id-1-0-dcpdimension(10*1+6,rank)/2))/double(dcpdimension(10*1+6,rank)/2);
      if (cuty < 1.0/3.0) {
//      cuty = double(dcpdimension(10*3+1,rank)+id-1+1)/double(dcpdimension(10*1+6,rank)-1);
//      if (cuty>1.0/3.0 && cuty<2.0/3.0) {
        GetCubeCplxJoin().subcube(0,
                              id,
                              0,
                              dcpdimension(10*3+6,rank)-1,
                              id,
                              dcpdimension(10*3+8,rank)*dimension_-1).fill(cx_double(0,0));
      }
    }
  
  }

  PS_DEBUG_TRACE_LEAVE(name_+"Dealiasing(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Export the x-pencil data first. Store the data locally, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - name      > target file name
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::SaveReal(const string& name, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"SaveReal(string,ftype_t)")
  
  dcube cuberealjoin;
  Export(cuberealjoin);
  if (mpi_.isMain()) {
    cuberealjoin.save(name,arma::file_type(type));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"SaveReal(string,ftype_t)")
  return true;
}
//bool RealCplxCube::SaveReal(const string& name, const ftype_t type, const int dimension) {
//  PS_DEBUG_TRACE_ENTER(name_+"SaveReal(string,ftype_t,int)")
//  
//  dcube cuberealjoin;
//  Export(cuberealjoin);
//  if (mpi_.isMain()) {
//    
//    const ivec& size = dcp_.GetSize();
//    dcube cuberealpart = cuberealjoin.slices(size(2)*dimension,size(2)*(dimension+1)-1);
//    cuberealpart.save(name,arma::file_type(type));
//    
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  
//  PS_DEBUG_TRACE_LEAVE(name_+"SaveReal(string,ftype_t,int)")
//  return true;
//}
  
  
  
//------------------------------------------------------------------------------
// Function : Export the z-pencil data first. Store the data locally, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - name      > target file name
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::SaveCplx(const string& name, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"SaveCplx(string,ftype_t)")
  
  ccube cubecplxjoin;
  Export(cubecplxjoin);
  if (mpi_.isMain()) {
    cubecplxjoin.save(name,arma::file_type(type));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"SaveCplx(string,ftype_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Load the local x-pencil data first. Scatter the data, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - name      > target file name
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::LoadReal(const string& name, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"LoadReal(string,ftype_t)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  dcube cuberealjoin;
  
  if (mpi_.isMain()) {
    cuberealjoin.set_size(dcpdimension(10*1+0,rank),
                          dcpdimension(10*1+1,rank),
                          dcpdimension(10*1+2,rank)*dimension_);
    cuberealjoin.load(name,arma::file_type(type));
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  *this = cuberealjoin;
  
  PS_DEBUG_TRACE_LEAVE(name_+"LoadReal(string,ftype_t)")
  return true;
}
  
  

//------------------------------------------------------------------------------
// Function : Load the local z-pencil data first. Scatter the data, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - name      > target file name
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::LoadCplx(const string& name, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"LoadCplx(string,ftype_t)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  ccube cubecplxjoin;
  
  if (mpi_.isMain()) {
    cubecplxjoin.set_size(dcpdimension(10*1+5,rank),
                          dcpdimension(10*1+6,rank),
                          dcpdimension(10*1+7,rank)*dimension_);
    cubecplxjoin.load(name,arma::file_type(type));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  *this = cubecplxjoin;
  
  PS_DEBUG_TRACE_LEAVE(name_+"LoadCplx(string,ftype_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Export the x-pencil data first. Store the data locally, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::SaveReal(fstream& stream, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"SaveReal(ofstream,ftype_t)")
  
  dcube cuberealjoin;
  Export(cuberealjoin);
  if (mpi_.isMain()) {
    cuberealjoin.save(stream,arma::file_type(type));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"SaveReal(ofstream,ftype_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Export the z-pencil data first. Store the data locally, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::SaveCplx(fstream& stream, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"SaveCplx(ofstream,ftype_t)")
  
  ccube cubecplxjoin;
  Export(cubecplxjoin);
  if (mpi_.isMain()) {
    cubecplxjoin.save(stream,arma::file_type(type));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  PS_DEBUG_TRACE_LEAVE(name_+"SaveCplx(ofstream,ftype_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Load the local x-pencil data first. Scatter the data, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::LoadReal(fstream& stream, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"LoadReal(ifstream,ftype_t)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  dcube cuberealjoin;
  
  if (mpi_.isMain()) {
    cuberealjoin.set_size(dcpdimension(10*1+0,rank),
                          dcpdimension(10*1+1,rank),
                          dcpdimension(10*1+2,rank)*dimension_);
    cuberealjoin.load(stream,arma::file_type(type));
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  *this = cuberealjoin;
  
  PS_DEBUG_TRACE_LEAVE(name_+"LoadReal(ifstream,ftype_t)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Load the local z-pencil data first. Scatter the data, then.
//            Do not call this function frequently.
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream
//            IN  - type      > target file type
//------------------------------------------------------------------------------
bool RealCplxCube::LoadCplx(fstream& stream, const ftype_t type) {
  PS_DEBUG_TRACE_ENTER(name_+"LoadCplx(ifstream,ftype_t)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  ccube cubecplxjoin;
  
  if (mpi_.isMain()) {
    cubecplxjoin.set_size(dcpdimension(10*1+5,rank),
                          dcpdimension(10*1+6,rank),
                          dcpdimension(10*1+7,rank)*dimension_);
    cubecplxjoin.load(stream,arma::file_type(type));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  *this = cubecplxjoin;
  
  PS_DEBUG_TRACE_LEAVE(name_+"LoadCplx(ifstream,ftype_t)")
  return true;
}
  
  
  
} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi
