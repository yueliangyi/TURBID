/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "domain.h"
#include <iomanip>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool Domain::Configuration::Load(IniFile& fileini) {
  PS_DEBUG_TRACE_ENTER("Domain::Configuration::Load(IniFile)")
  
  dvec length = fileini.GetVectorDouble(SECTION_DOMAIN,KEY_SD_DIMENSIONLESS_LENGTH);
  lengthx = length(0);
  lengthy = length(1);
  lengthz = length(2);
  
  slope       = fileini.GetValueDouble(SECTION_DOMAIN,KEY_SD_TOP_BOTTOM_SLOPE);
  slope_direction = fileini.GetValueInt(SECTION_DOMAIN,KEY_SD_TOP_BOTTOM_SLOPE_DIRECTION);
  readprofile = fileini.GetValueBool(SECTION_DOMAIN,KEY_SD_READ_PROFILE_FROM_FILE);
  
//  rippleamplitude = fileini.GetValueDouble(SECTION_DOMAIN,KEY_SD_RIPPLE_BED_AMPLITUDE);
  
  PS_DEBUG_TRACE_LEAVE("Domain::Configuration::Load(IniFile)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Domain::Domain(const mpi_t& mpi, const dcp_t& dcp, const Base::cfg_t* cfg) :
  Base(mpi,dcp,cfg)
//  , physicalgrid_(mpi_,dcp_,3)
  {
  PS_DEBUG_TRACE_ENTER(name_+"Domain(mpi_t,dcp_t,cfg_t)")
    
    cfg_t* pcfg = (cfg_t*)cfg_;
    
#ifdef PS_DEBUG_CHECK
  if (pcfg == NULL) {
    name_ += "Domain(mpi_t,dcp_t,cfg_t) > Can not Initialize Domain Without Configuration!";
    throw std::runtime_error(name_);
  }
  if (abs(pcfg->lengthz-2) > COMPARE_NUMBER_SMALL) {
    name_ += "Domain(mpi_t,dcp_t,cfg_t) > Dimensionless z-length Should be Equal to Two!";
    throw std::runtime_error(name_);
  }
#endif
  
  const ivec& size = dcp_.GetSize();
  gridx_ = arma::zeros<dvec>(size(0));
  gridy_ = arma::zeros<dvec>(size(1));
  gridz_ = arma::zeros<dvec>(size(2));
  gridzcluster_ = arma::zeros<dvec>(size(2));

    
    //
    if (mpi_.isMain()) {
      profiletop_    = arma::zeros(size(0),size(1));
      profilebottom_ = arma::zeros(size(0),size(1));
    }
    
    
    etam_.push_back(new rccube(mpi_,dcp_,1));   // etam
    etam_.push_back(new rccube(mpi_,dcp_,1));   // detamdx1
    etam_.push_back(new rccube(mpi_,dcp_,1));   // detamdy1
    etam_.push_back(new rccube(mpi_,dcp_,1));   // detamdt
    etam_.push_back(new rccube(mpi_,dcp_,1));   // detamdx2
    etam_.push_back(new rccube(mpi_,dcp_,1));   // detamdy2
    
    etar_.push_back(new rccube(mpi_,dcp_,1));   // etar
    etar_.push_back(new rccube(mpi_,dcp_,1));   // detardx1
    etar_.push_back(new rccube(mpi_,dcp_,1));   // detardy1
    etar_.push_back(new rccube(mpi_,dcp_,1));   // detardt
    etar_.push_back(new rccube(mpi_,dcp_,1));   // detardx2
    etar_.push_back(new rccube(mpi_,dcp_,1));   // detardy2
    
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // J31
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // J32
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // J33
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // J34
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // d2xi3dx2
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // d2xi3dy2
    jm3k_.push_back(new rccube(mpi_,dcp_,1));   // J31^2+J32^2+J33^2
    
    
    phy_grid_x1_.push_back(new rccube(mpi_,dcp_,1));
    phy_grid_x2_.push_back(new rccube(mpi_,dcp_,1));
    phy_grid_x3_.push_back(new rccube(mpi_,dcp_,1));   // etab
    phy_grid_x3_.push_back(new rccube(mpi_,dcp_,1));   // detabdx1
    phy_grid_x3_.push_back(new rccube(mpi_,dcp_,1));   // detabdy1
    phy_grid_x3_.push_back(new rccube(mpi_,dcp_,1));   // temp

    dvec_normal_.push_back(new rccube(mpi_,dcp_,1));
    dvec_normal_.push_back(new rccube(mpi_,dcp_,1));
    dvec_normal_.push_back(new rccube(mpi_,dcp_,1));

    dvec_tangent1_.push_back(new rccube(mpi_,dcp_,1));
    dvec_tangent1_.push_back(new rccube(mpi_,dcp_,1));
    dvec_tangent1_.push_back(new rccube(mpi_,dcp_,1));
    
    dvec_tangent2_.push_back(new rccube(mpi_,dcp_,1));
    dvec_tangent2_.push_back(new rccube(mpi_,dcp_,1));
    dvec_tangent2_.push_back(new rccube(mpi_,dcp_,1));
    
    
    // Set the name of output file.
    filename_ = string("phase_00_")+IO_FILE_CLASS_DOMAIN;
    
  PS_DEBUG_TRACE_LEAVE(name_+"Domain(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Domain::~Domain(void) {
  for (auto it=etam_.begin(); it!=etam_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=etar_.begin(); it!=etar_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=jm3k_.begin(); it!=jm3k_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=phy_grid_x1_.begin(); it!=phy_grid_x1_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=phy_grid_x2_.begin(); it!=phy_grid_x2_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=phy_grid_x3_.begin(); it!=phy_grid_x3_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=dvec_normal_.begin(); it!=dvec_normal_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=dvec_tangent1_.begin(); it!=dvec_tangent1_.end(); ++it) {
    delete *it; *it = NULL;
  }
  for (auto it=dvec_tangent2_.begin(); it!=dvec_tangent2_.end(); ++it) {
    delete *it; *it = NULL;
  }
}

  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool Domain::Update(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Update(void)")
  
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  
  // construct the others from top and bottom profiles

  
//  // physicalgrid_
//  for (int id = 0; id < dcpdimension(10*2+6,rank); ++id) {
//    physicalgrid_.GetCubeReal(0).tube(id,
//                                      0,
//                                      id,
//                                      dcpdimension(10*2+7,rank)-1).fill(gridx_(dcpdimension(10*2+0,rank)-1+id));
//  }
//  for (int id = 0; id < dcpdimension(10*2+7,rank); ++id) {
//    physicalgrid_.GetCubeReal(1).tube(0,
//                                      id,
//                                      dcpdimension(10*2+6,rank)-1,
//                                      id).fill(gridy_(dcpdimension(10*2+1,rank)-1+id));
//  }
//  for (int id = 0; id < dcpdimension(10*2+8,rank); ++id) {
//    physicalgrid_.GetCubeReal(2).slice(id) =
//      (1+etar_[0]->GetCubeReal().slice(id))*gridz_(dcpdimension(10*2+2,rank)-1+id)+
//      etam_[0]->GetCubeReal().slice(id);
//  }

  // xi3k_
  for (int id = 0; id < dcpdimension(10*2+8,rank); ++id) {
    
    jm3k_[0]->GetCubeReal().slice(id) =
    -(etam_[1]->GetCubeReal().slice(id)+
      etar_[1]->GetCubeReal().slice(id)*gridz_(dcpdimension(10*2+2,rank)-1+id))/
    (1+etar_[0]->GetCubeReal().slice(id));
    
    jm3k_[1]->GetCubeReal().slice(id) =
    -(etam_[2]->GetCubeReal().slice(id)+
      etar_[2]->GetCubeReal().slice(id)*gridz_(dcpdimension(10*2+2,rank)-1+id))/
    (1+etar_[0]->GetCubeReal().slice(id));
    
    jm3k_[2]->GetCubeReal().slice(id) =
    +1.0/
    (1+etar_[0]->GetCubeReal().slice(id));
    
    jm3k_[3]->GetCubeReal().slice(id) =
    -(etam_[3]->GetCubeReal().slice(id)+
      etar_[3]->GetCubeReal().slice(id)*gridz_(dcpdimension(10*2+2,rank)-1+id))/
    (1+etar_[0]->GetCubeReal().slice(id));
    
    jm3k_[4]->GetCubeReal().slice(id) =
    -(etam_[4]->GetCubeReal().slice(id)+
      etar_[4]->GetCubeReal().slice(id)*gridz_(dcpdimension(10*2+2,rank)-1+id)+
      2*jm3k_[0]->GetCubeReal().slice(id)%etar_[1]->GetCubeReal().slice(id))/
    (1+etar_[0]->GetCubeReal().slice(id));
    
    jm3k_[5]->GetCubeReal().slice(id) =
    -(etam_[5]->GetCubeReal().slice(id)+
      etar_[5]->GetCubeReal().slice(id)*gridz_(dcpdimension(10*2+2,rank)-1+id)+
      2*jm3k_[1]->GetCubeReal().slice(id)%etar_[2]->GetCubeReal().slice(id))/
    (1+etar_[0]->GetCubeReal().slice(id));
    
  }
  
  
  
  jm3k_[6]->GetCubeReal() =
    jm3k_[0]->GetCubeReal()%jm3k_[0]->GetCubeReal()+
    jm3k_[1]->GetCubeReal()%jm3k_[1]->GetCubeReal()+
    jm3k_[2]->GetCubeReal()%jm3k_[2]->GetCubeReal();
  
  
  
  
//  jm3k_[0]->SaveReal("jm3k0",rccube::kArmaAscii);
//  jm3k_[6]->SaveReal("jm3k6",rccube::kArmaAscii);
  
  
//  jm3k_[0]->SaveCplx("jm3k0",rccube::kArmaBinary);
//  jm3k_[1]->SaveCplx("jm3k1",rccube::kArmaBinary);
//  jm3k_[2]->SaveCplx("jm3k2",rccube::kArmaBinary);
  
  
  
//  // remove aliasing error
//  // or use values from spectral domain
//  for (int id = 0; id < 7; ++id) {
//    jm3k_[id]->DftForward();
//    jm3k_[id]->DftBackward();
//  }
  
  
  
  // build normal vectors of bottom surface
  phy_grid_x3_[3]->GetCubeReal() = arma::sqrt(phy_grid_x3_[1]->GetCubeReal()%phy_grid_x3_[1]->GetCubeReal()+
                                              phy_grid_x3_[2]->GetCubeReal()%phy_grid_x3_[2]->GetCubeReal()+
                                              1);
  dvec_normal_[0]->GetCubeReal() = -1*phy_grid_x3_[1]->GetCubeReal()/phy_grid_x3_[3]->GetCubeReal();
  dvec_normal_[1]->GetCubeReal() = -1*phy_grid_x3_[2]->GetCubeReal()/phy_grid_x3_[3]->GetCubeReal();
  dvec_normal_[2]->GetCubeReal() =  1/phy_grid_x3_[3]->GetCubeReal();
  
  
  // build tangent vector in x1 of bottom surface
  phy_grid_x3_[3]->GetCubeReal() = arma::sqrt(phy_grid_x3_[1]->GetCubeReal()%phy_grid_x3_[1]->GetCubeReal()+
                                              1);
  dvec_tangent1_[0]->GetCubeReal() = 1/phy_grid_x3_[3]->GetCubeReal();
  dvec_tangent1_[1]->GetCubeReal() = 0/phy_grid_x3_[3]->GetCubeReal();
  dvec_tangent1_[2]->GetCubeReal() = phy_grid_x3_[1]->GetCubeReal()/phy_grid_x3_[3]->GetCubeReal();
  
  
  // build tangent vector in x1 of bottom surface
  phy_grid_x3_[3]->GetCubeReal() = arma::sqrt(phy_grid_x3_[2]->GetCubeReal()%phy_grid_x3_[2]->GetCubeReal()+
                                              1);
  dvec_tangent2_[0]->GetCubeReal() = 0/phy_grid_x3_[3]->GetCubeReal();
  dvec_tangent2_[1]->GetCubeReal() = 1/phy_grid_x3_[3]->GetCubeReal();
  dvec_tangent2_[2]->GetCubeReal() = phy_grid_x3_[2]->GetCubeReal()/phy_grid_x3_[3]->GetCubeReal();
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Update(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool Domain::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  const ivec& size = dcp_.GetSize();
  
  // 2. Compute the location of grid point.
  //    The grid is uniform in x- and y-direction, while the CGL (Chebyshev
  //    Gauss Lobatto) points are used in z-direction:
  //
  //              z_l = cos(l*pi/N_z) with 0 <= l <= N_z
  //
  for (int id = 0; id < size(0); ++id) { gridx_(id) = pcfg->lengthx/size(0)*id; }
  for (int id = 0; id < size(1); ++id) { gridy_(id) = pcfg->lengthy/size(1)*id; }
  for (int id = 0; id < size(2); ++id) { gridz_(id) = cos((size(2)-id-1)*arma::datum::pi/(size(2)-1)); }
  
  // Get the minimum length of doamin grid.
  lengthminx_ = diff(gridx_).min();
  lengthminy_ = diff(gridy_).min();
  lengthminz_ = diff(gridz_).min();
  
  
  // Build the grid length matrix
  const int rank = mpi_.GetRank();
  const imat& dcpdimension = dcp_.GetDimension();
  maxulen_  = arma::zeros(dcpdimension(10*2+6,rank),dcpdimension(10*2+7,rank),dcpdimension(10*2+8,rank));
  gridlenx_ = arma::ones(dcpdimension(10*2+6,rank),dcpdimension(10*2+7,rank),dcpdimension(10*2+8,rank))*lengthminx_;
  gridleny_ = arma::ones(dcpdimension(10*2+6,rank),dcpdimension(10*2+7,rank),dcpdimension(10*2+8,rank))*lengthminy_;
  gridlenz_ = arma::ones(dcpdimension(10*2+6,rank),dcpdimension(10*2+7,rank),dcpdimension(10*2+8,rank));
  for (int id = 0; id < dcpdimension(10*2+8,rank); ++id) {
    if (dcpdimension(10*2+2,rank) == 1 && id == 0) {
      gridlenz_.slice(id).fill((gridz_(1)-gridz_(0))/2);
      continue;
    }
    if (dcpdimension(10*2+5,rank) == size(2) && id == dcpdimension(10*2+8,rank)-1) {
      gridlenz_.slice(id).fill((gridz_(size(2)-1)-gridz_(size(2)-2))/2);
      continue;
    }
    gridlenz_.slice(id).fill(std::min(gridz_(dcpdimension(10*2+2,rank)-1+id+1)-gridz_(dcpdimension(10*2+2,rank)-1+id),
                                      gridz_(dcpdimension(10*2+2,rank)-1+id)-gridz_(dcpdimension(10*2+2,rank)-1+id-1)));
  }
  
  
  
  
  // load top and bottom profiles
  if (mpi_.isMain() && pcfg->readprofile) {
    profilebottom_.load("profile_bottom.dat",arma::arma_binary);
    profiletop_.load("profile_top.dat",arma::arma_binary);
  }
  
  
  isconfigured_ = true;
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Print current status
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream as container
//------------------------------------------------------------------------------
void Domain::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  
  stream
  << "Length of Domain ["
  << std::setw(6) << pcfg->lengthx << " "
  << std::setw(6) << pcfg->lengthy << " "
  << std::setw(6) << pcfg->lengthz << "]"
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Minimum Length of Grid ["
  << std::setw(6) << lengthminx_ << " "
  << std::setw(6) << lengthminy_ << " "
  << std::setw(6) << lengthminz_ << "]"
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Bottom Slope "
  << std::setw(6) << pcfg->slope;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
