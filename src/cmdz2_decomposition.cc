/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "cmdz2_decomposition.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool CmDz2DecompositionCplx::Construct(const SpaceStepper* sstepper,
                                       const dvec& zboundarytype,
                                       const dvec& stepcoefficient,
                                       const int dimension,
                                       const bool includedxy2) {
  PS_DEBUG_TRACE_ENTER("CmDz2DecompositionCplx::Construct(SpaceStepper,dmat,dmat)")
  
  
  dimension_ = dimension;
  includedxy2_ = includedxy2;
  
  
  const int sizez = sstepper->GetDcp().GetSize(2);
  dmat cmdz1 = arma::real(sstepper->GetCmDz1());
  dmat cmdz2 = arma::real(sstepper->GetCmDz2());
  
  
  // Create the original coefficient matrix with boundary conditions.
  // For inside region:
  //                      partial2(f)
  //                      ----------- + b*f = c
  //                      partial2(z)
  //
  // For boundaries:
  //                      partial(f) |
  //              a*f + b*---------- |         = c|_{z=+-1}
  //                      partial(z) |_{z=+-1}
  //
  cmdz2.row(0      ) = zboundarytype(1)*cmdz1.row(0      );
  cmdz2.row(sizez-1) = zboundarytype(3)*cmdz1.row(sizez-1);
  cmdz2(0      ,0      ) += zboundarytype(0);
  cmdz2(sizez-1,sizez-1) += zboundarytype(2);
  
  
  
  const double a11 = cmdz2.at(0      ,0      );
  const double ann = cmdz2.at(sizez-1,sizez-1);
  const double a1n = cmdz2.at(0      ,sizez-1);
  const double an1 = cmdz2.at(sizez-1,0      );
  const double den = 1.0/(an1*a1n-a11*ann);
  
  
  
  // Generate fixing coefficients for forward and backward steps.
  dmat fixforwardbottom (1,sizez,arma::fill::zeros);
  dmat fixforwardtop    (1,sizez,arma::fill::zeros);
  dmat fixbackwardbottom(1,sizez,arma::fill::zeros);
  dmat fixbackwardtop   (1,sizez,arma::fill::zeros);
  
  
  for (int ix = 1; ix < sizez-1; ++ix) {
    
    fixbackwardbottom(ix) = (ann*cmdz2.at(0,ix)-a1n*cmdz2.at(sizez-1,ix))*den;
    fixbackwardtop   (ix) = (a11*cmdz2.at(sizez-1,ix)-an1*cmdz2.at(0,ix))*den;
    fixforwardbottom (ix) = (ann*cmdz2.at(ix,0)-an1*cmdz2.at(ix,sizez-1))*den;
    fixforwardtop    (ix) = (a11*cmdz2.at(ix,sizez-1)-a1n*cmdz2.at(ix,0))*den;
    
  }
  fixbackwardbottom(0      ) = -1*ann*den;
  fixbackwardbottom(sizez-1) = +1*a1n*den;
  fixbackwardtop   (0      ) = +1*an1*den;
  fixbackwardtop   (sizez-1) = -1*a11*den;
  
  
  
  // Generate the forward coefficient matrix.
  dmat cmdz2forward(sizez-2,sizez-2,arma::fill::zeros);
  for (int iz = 1; iz < sizez-1; ++iz) {
    for (int ix = 1; ix < sizez-1; ++ix) {
      cmdz2forward(ix-1,iz-1) =
        cmdz2.at(ix,iz) +
        cmdz2.at(ix,0      )*fixbackwardbottom(iz) +
        cmdz2.at(ix,sizez-1)*fixbackwardtop(iz);
    }
  }
  
  
  // Get the eigen-decomposition information.
  arma::cx_vec eigenvalue;
  arma::cx_mat eigenvector;
  eig_gen(eigenvalue,eigenvector,cmdz2forward);
  eigenvalue.set_imag(zeros(size(eigenvalue)));
  
  // Get real eigen-vector matrix.
  dmat eigenvectornorm(sizez,sizez,arma::fill::zeros);
  eigenvectornorm.submat(1,1,sizez-2,sizez-2) = arma::real(eigenvector);
  eigenvectornorm(0      ,0      ) = 1.0;
  eigenvectornorm(sizez-1,sizez-1) = 1.0;
  
  // Get real eigen-vector inverse matrix.
  dmat eigenvectorinverse(sizez,sizez,arma::fill::zeros);
  eigenvectorinverse.submat(1,1,sizez-2,sizez-2) = arma::inv(arma::real(eigenvector));
  eigenvectorinverse(0      ,0      ) = 1.0;
  eigenvectorinverse(sizez-1,sizez-1) = 1.0;
  
  
  
  
  //
  zboundarytype_ = zboundarytype;
  stepcoefficient_ = stepcoefficient;
  coefficientmatrix_ = cmdz2;
  
  coefficientmatrix_.zeros();
  coefficientmatrix_.submat(1,1,sizez-2,sizez-2) = cmdz2forward;
  coefficientmatrix_(0,0) = 1;
  coefficientmatrix_(sizez-1,sizez-1) = 1;
  
  
  
  fixforwardtop_ = cmat(fixforwardtop,zeros(size(fixforwardtop)));
  fixforwardbottom_ = cmat(fixforwardbottom,zeros(size(fixforwardbottom)));
  fixbackwardtop_ = cmat(fixbackwardtop,zeros(size(fixbackwardtop)));
  fixbackwardbottom_ = cmat(fixbackwardbottom,zeros(size(fixbackwardbottom)));
  
  eigenvectornorm_ = cmat(eigenvectornorm,zeros(size(eigenvectornorm)));
  eigenvectorinverse_ = cmat(eigenvectorinverse,zeros(size(eigenvectorinverse)));

  
  eigenvalue_.set_size(sizez,1);
  eigenvalue_.rows(1,sizez-2) = eigenvalue;
  
  
  
  cmdx2_ = arma::real(sstepper->GetCmDx2().GetCMatCplx());
  cmdy2_ = arma::real(sstepper->GetCmDy2().GetCMatCplx());
  
  
  const int stepnum = std::max(1,(int)stepcoefficient.n_elem);
  for (int idstep = 0; idstep < stepnum; ++idstep) {
    
//    double stepvalue = 0;
//    if (!stepcoefficient.is_empty()) {
//      stepvalue = stepcoefficient(idstep);
//    }
    
    
    rccube* reciprocal = new rccube(sstepper->GetMpi(),sstepper->GetDcp(),dimension_);
    eigenvaluereciprocal_.push_back(reciprocal);

    
    
    
    
    
//    for (int id = 1; id < sizez-1; ++id) {
////      reciprocal->GetCMatCplx().row(id).fill(cx_double(1/(eigenvalue_(id).real()+stepvalue),0));
//      for (size_t ixy = 0; ixy < cmdx2_.n_cols; ++ixy) {
//        if (includedxy2_) {
//          reciprocal->GetCMatCplx(0)(id,ixy) =
//          cx_double(1/(cmdx2_(id,ixy)+cmdy2_(id,ixy)+eigenvalue_(id).real()+stepvalue),0);
//        } else {
//          reciprocal->GetCMatCplx(0)(id,ixy) =
//          cx_double(1/(eigenvalue_(id).real()+stepvalue),0);
//        }
//      }
//
//
//    }
//
//
//    for (int iddim = 0; iddim < dimension_-1; ++iddim) {
//      reciprocal->GetCMatCplx(iddim+1) = reciprocal->GetCMatCplx(0);
//    }
//
//
//    reciprocal->GetCMatCplxJoin().row(0      ).fill(cx_double(1,0));
//    reciprocal->GetCMatCplxJoin().row(sizez-1).fill(cx_double(1,0));
//
//
//    //
//    if (fabs(stepvalue) < COMPARE_NUMBER_SMALL &&
//        sstepper->GetDcp().isMeanZ()) {
//
//      for (int iddim = 0; iddim < dimension_; ++iddim) {
//
////        reciprocal->GetCMatCplx(iddim).col(0).fill(cx_double(0,0));
//
//
//
//
//
//        // modified for test!!!
//        for (int idz = 1; idz < sizez-1; ++idz) {
//          if (abs(eigenvalue_(idz).real()) < COMPARE_NUMBER_SMALL) {
//            reciprocal->GetCMatCplx(iddim)(idz,0) = cx_double(0,0);
////            std::cout << idz << std::endl;
//          }
//        }
//
//
//      }
//
//    }
    
    
    
  }
  
  
  Assemble(sstepper,stepcoefficient);

  
  
  
  
  PS_DEBUG_TRACE_LEAVE("CmDz2DecompositionCplx::Construct(SpaceStepper,dmat,dmat)")
  return true;
}
  
  
  
//
//
//
bool CmDz2DecompositionCplx::Assemble(const SpaceStepper* sstepper, const dvec& stepcoefficient) {
  
  const int sizez = sstepper->GetDcp().GetSize(2);

  const int stepnum = std::max(1,(int)stepcoefficient.n_elem);
  for (int idstep = 0; idstep < stepnum; ++idstep) {
    
    double stepvalue = 0;
    if (!stepcoefficient.is_empty()) {
      stepvalue = stepcoefficient(idstep);
    }
    
    rccube* reciprocal = eigenvaluereciprocal_[idstep];
  
    for (size_t ixy = 0; ixy < cmdx2_.n_cols; ++ixy) {
      for (int id = 1; id < sizez-1; ++id) {
        if (includedxy2_) {
          reciprocal->GetCMatCplx(0)(id,ixy) =
          cx_double(1/(cmdx2_(id,ixy)+cmdy2_(id,ixy)+eigenvalue_(id).real()+stepvalue),0);
        } else {
          reciprocal->GetCMatCplx(0)(id,ixy) =
          cx_double(1/(eigenvalue_(id).real()+stepvalue),0);
        }
      }
      
      
    }
  
  
    for (int iddim = 0; iddim < dimension_-1; ++iddim) {
      reciprocal->GetCMatCplx(iddim+1) = reciprocal->GetCMatCplx(0);
    }
  
  
    reciprocal->GetCMatCplxJoin().row(0      ).fill(cx_double(1,0));
    reciprocal->GetCMatCplxJoin().row(sizez-1).fill(cx_double(1,0));
  
  
    //
    if (fabs(stepvalue) < COMPARE_NUMBER_SMALL &&
        sstepper->GetDcp().isMeanZ()) {
      
      for (int iddim = 0; iddim < dimension_; ++iddim) {
        
//        reciprocal->GetCMatCplx(iddim).col(0).fill(cx_double(0,0));
        
        // modified for test!!!
        for (int idz = 1; idz < sizez-1; ++idz) {
          if (abs(eigenvalue_(idz).real()) < COMPARE_NUMBER_SMALL) {
            reciprocal->GetCMatCplx(iddim)(idz,0) = cx_double(0,0);
//            std::cout << idz << std::endl;
          }
        }
        
        
      }
      
    }
    
  }
  
  return true;
    
}

  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
CmDz2DecompositionCplx::~CmDz2DecompositionCplx(void) {
  for (auto it=eigenvaluereciprocal_.begin(); it!=eigenvaluereciprocal_.end(); ++it) {
    delete *it; *it = NULL;
  }
}
  
  
  
  
  

  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
