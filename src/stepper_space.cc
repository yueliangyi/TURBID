/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "stepper_space.h"
#include <iomanip>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool SpaceStepper::Configuration::Load(IniFile& fileini) {
  PS_DEBUG_TRACE_ENTER(name_+"Configuration::Load(IniFile)")
  
  string function = ConvertCase(fileini.GetValueString(SECTION_STEPPER_SPACE,KEY_SSS_MAPPING_FUNCTION),1);
  
  if (function == "NONE") { mappingfunction = kNone; }
  else if (function == "TWOENDTANGENT") { mappingfunction = kTwoEndTangent; }
  else if (function == "TWOENDSIN") { mappingfunction = kTwoEndSin; }
  else if (function == "ONEENDTANGENT") { mappingfunction = kOneEndTangent; }
  else {
    throw std::runtime_error(string("SpaceStepper::Configuration::")+
                             "Load(IniFile) > Unsupported Mapping Function!");
  }
  
  mappingcoefficient = fileini.GetVectorDouble(SECTION_STEPPER_SPACE,KEY_SSS_MAPPING_COEFFICIENT);
  PS_DEBUG_TRACE_LEAVE(name_+"Configuration::Load(IniFile)")
  return true;
}
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
SpaceStepper::SpaceStepper(const mpi_t& mpi,
                           const dcp_t& dcp,
                           Domain* domain,
                           const Base::cfg_t* cfg) :
  Base(mpi,dcp,cfg),domain_(domain),
  cmdx1_(mpi,dcp),cmdy1_(mpi,dcp),cmdxy1_(mpi,dcp),
  cmdx2_(mpi,dcp),cmdy2_(mpi,dcp),cmdxy2_(mpi,dcp) {
  PS_DEBUG_TRACE_ENTER(name_+"SpaceStepper(mpi_t,dcp_t,Domain,cfg_t)")
    
    
    filename_ = string("phase_00_")+IO_FILE_CLASS_STEPPER_SPACE;
    
    
  PS_DEBUG_TRACE_LEAVE(name_+"SpaceStepper(mpi_t,dcp_t,Domain,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool SpaceStepper::Update(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Update(void)")
  
  const ivec& size = dcp_.GetSize();
  
  // construct the others from top and bottom profiles
  
  // etam_ and etar_
  dcube etam;
  dcube etar;
  if (mpi_.isMain()) {
    etam = dcube(size(0),size(1),size(2));
    etar = dcube(size(0),size(1),size(2));
    etam.each_slice() = (domain_->GetProfileTop()+domain_->GetProfileBottom())/2;
    etar.each_slice() = (domain_->GetProfileTop()-domain_->GetProfileBottom())/2;
  }
  
  domain_->GetEtaM(0) = etam;
  domain_->GetEtaR(0) = etar;
//  domain_->GetEtaM(0).SaveReal("etamreal",rccube::kArmaAscii);
//  domain_->GetEtaR(0).SaveReal("etarreal",rccube::kArmaAscii);
  
  
  domain_->GetEtaM(0).DftForward();
  domain_->GetEtaR(0).DftForward();
//  domain_->GetEtaM(0).SaveCplx("etamcplx",rccube::kArmaAscii);
//  domain_->GetEtaR(0).SaveCplx("etarcplx",rccube::kArmaAscii);
  
  
  DDx1(domain_->GetEtaM(0),domain_->GetEtaM(1));
  DDx1(domain_->GetEtaR(0),domain_->GetEtaR(1));
  domain_->GetEtaM(1).DftBackward();
  domain_->GetEtaR(1).DftBackward();
  
//  domain_->GetEtaM(1).SaveReal("detamdx1",rccube::kArmaAscii);
//  domain_->GetEtaR(1).SaveReal("detardx1",rccube::kArmaAscii);
  
  
  DDy1(domain_->GetEtaM(0),domain_->GetEtaM(2));
  DDy1(domain_->GetEtaR(0),domain_->GetEtaR(2));
  domain_->GetEtaM(2).DftBackward();
  domain_->GetEtaR(2).DftBackward();
  
  DDx2(domain_->GetEtaM(0),domain_->GetEtaM(4));
  DDx2(domain_->GetEtaR(0),domain_->GetEtaR(4));
  domain_->GetEtaM(4).DftBackward();
  domain_->GetEtaR(4).DftBackward();
//  domain_->GetEtaM(4).SaveReal("detamdx2",rccube::kArmaAscii);
//  domain_->GetEtaR(4).SaveReal("detardx2",rccube::kArmaAscii);
  
  DDy2(domain_->GetEtaM(0),domain_->GetEtaM(5));
  DDy2(domain_->GetEtaR(0),domain_->GetEtaR(5));
  domain_->GetEtaM(5).DftBackward();
  domain_->GetEtaR(5).DftBackward();
  
  
  
  dcube phy_grid_x1;
  dcube phy_grid_x2;
  dcube phy_grid_x3;
  if (mpi_.isMain()) {
    
    phy_grid_x1 = dcube(size(0),size(1),size(2));
    phy_grid_x2 = dcube(size(0),size(1),size(2));
    for (int id3 = 0; id3 < size(2); ++id3) {
      for (int id2 = 0; id2 < size(1); ++id2) {
        for (int id1 = 0; id1 < size(0); ++id1) {
          phy_grid_x1(id1,id2,id3) = domain_->GetGridX()(id1);
          phy_grid_x2(id1,id2,id3) = domain_->GetGridY()(id2);
        }
      }
    }
    
    phy_grid_x3 = dcube(size(0),size(1),size(2));
    phy_grid_x3 = 1+etar;
    const dvec& gridz = domain_->GetGridZ();
    for (int id = 0; id < size(2); ++id) {
      phy_grid_x3.slice(id) = phy_grid_x3.slice(id)*gridz(id);
    }
    phy_grid_x3 = phy_grid_x3+etam;
    
  }
  
  domain_->GetPhyGridX1(0) = phy_grid_x1;
  domain_->GetPhyGridX2(0) = phy_grid_x2;
  domain_->GetPhyGridX1(0).DftForward();
  domain_->GetPhyGridX2(0).DftForward();

  domain_->GetPhyGridX3(0) = phy_grid_x3;
  domain_->GetPhyGridX3(0).DftForward();
  
  DDx1(domain_->GetPhyGridX3(0),domain_->GetPhyGridX3(1));
  DDy1(domain_->GetPhyGridX3(0),domain_->GetPhyGridX3(2));
  
  domain_->GetPhyGridX3(1).DftBackward();
  domain_->GetPhyGridX3(2).DftBackward();
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Update(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool SpaceStepper::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  const int rank = mpi_.GetRank();
  const imat& dimension = dcp_.GetDimension();
  const int sizez = dimension(10*1+5,rank);
  const int sizey = dimension(10*1+6,rank);
  const int sizex = dimension(10*1+7,rank);
  
  
  const double kx1 = 2*arma::datum::pi/domain_->GetLengthX(); double kxn = 0;
  const double ky1 = 2*arma::datum::pi/domain_->GetLengthY(); double kym = 0;
  
  
  ccube cmdx1;
  ccube cmdy1;
  ccube cmdx2;
  ccube cmdy2;
  
  if (mpi_.isMain()) {
  
    cmdx1.set_size(sizez,sizey,sizex);
    cmdy1.set_size(sizez,sizey,sizex);
    cmdx2.set_size(sizez,sizey,sizex);
    cmdy2.set_size(sizez,sizey,sizex);
    
    
    for (int iz = 0; iz < sizez; ++iz) {
      for (int iy = 0; iy < sizey; ++iy) {
        for (int ix = 0; ix < sizex; ++ix) {
          
          kxn = kx1*(ix-0);
          cmdx1(iz,iy,ix) = cx_double(0,+pow(kxn,1));
          cmdx2(iz,iy,ix) = cx_double(-pow(kxn,2),0);
          
          if (iy <= sizey/2) { kym = +ky1*(iy-0); }
          else { kym = -ky1*(sizey-iy); }
          cmdy1(iz,iy,ix) = cx_double(0,+pow(kym,1));
          cmdy2(iz,iy,ix) = cx_double(-pow(kym,2),0);
          
        } // end of ix
      } // end of iy
    } // end of iz
    
  } // end of if
  
  
  ccube cmdxy1 = cmdx1+cmdy1;
  ccube cmdxy2 = cmdx2+cmdy2;
  
  cmdx1_ = cmdx1;
  cmdy1_ = cmdy1;
  cmdxy1_ = cmdxy1;
  cmdx2_ = cmdx2;
  cmdy2_ = cmdy2;
  cmdxy2_ = cmdxy2;
  
  
  
//  if (mpi_.isMain()) {
//    cmdx1.save("cmdx1",arma::arma_ascii);
//    cmdy1.save("cmdy1",arma::arma_ascii);
//    cmdx2.save("cmdx2",arma::arma_ascii);
//    cmdy2.save("cmdy2",arma::arma_ascii);
//  }
  
  
  
  const ccube& tmpx1 = cmdx1_.GetCubeCplxJoin();
  const ccube& tmpy1 = cmdy1_.GetCubeCplxJoin();
  
  cmdxxx1_ = join_slices(join_slices(tmpx1,tmpx1),tmpx1);
  cmdyyy1_ = join_slices(join_slices(tmpy1,tmpy1),tmpy1);
  
  const ccube& tmpx2 = cmdx2_.GetCubeCplxJoin();
  const ccube& tmpy2 = cmdy2_.GetCubeCplxJoin();
  
  cmdxxx2_ = join_slices(join_slices(tmpx2,tmpx2),tmpx2);
  cmdyyy2_ = join_slices(join_slices(tmpy2,tmpy2),tmpy2);
  
  
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
  
// Function : Get the first and second order derivatives
// Parameter: IN  - tar     > target
//            OUT - dd**    > result
// Return   : True if success
bool SpaceStepper::DDx1(rccube& target, rccube& ddx1) const {
  PS_DEBUG_TRACE_ENTER(name_+"DDx1(rccube,rccube)")
  for (int id = 0; id < ddx1.GetDimension(); ++id) {
    ddx1.GetCubeCplx(id) = GetCmDx1().GetCubeCplx()%target.GetCubeCplx(id);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"DDx1(rccube,rccube)")
  return true;
}
bool SpaceStepper::DDx2(rccube& target, rccube& ddx2) const {
  PS_DEBUG_TRACE_ENTER(name_+"DDx2(rccube,rccube)")
  for (int id = 0; id < ddx2.GetDimension(); ++id) {
    ddx2.GetCubeCplx(id) = GetCmDx2().GetCubeCplx()%target.GetCubeCplx(id);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"DDx2(rccube,rccube)")
  return true;
}
bool SpaceStepper::DDy1(rccube& target, rccube& ddy1) const {
  PS_DEBUG_TRACE_ENTER(name_+"DDy1(rccube,rccube)")
  for (int id = 0; id < ddy1.GetDimension(); ++id) {
    ddy1.GetCubeCplx(id) = GetCmDy1().GetCubeCplx()%target.GetCubeCplx(id);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"DDy1(rccube,rccube)")
  return true;
}
bool SpaceStepper::DDy2(rccube& target, rccube& ddy2) const {
  PS_DEBUG_TRACE_ENTER(name_+"DDy2(rccube,rccube)")
  for (int id = 0; id < ddy2.GetDimension(); ++id) {
    ddy2.GetCubeCplx(id) = GetCmDy2().GetCubeCplx()%target.GetCubeCplx(id);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"DDy2(rccube,rccube)")
  return true;
}
bool SpaceStepper::DDz1(rccube& target, rccube& ddz1) const {
  PS_DEBUG_TRACE_ENTER(name_+"DDz1(rccube,rccube)")
  target.Multiply(cmdz1_,ddz1);
  PS_DEBUG_TRACE_LEAVE(name_+"DDz1(rccube,rccube)")
  return true;
}
bool SpaceStepper::DDz2(rccube& target, rccube& ddz2) const {
  PS_DEBUG_TRACE_ENTER(name_+"DDz2(rccube,rccube)")
  target.Multiply(cmdz2_,ddz2);
  PS_DEBUG_TRACE_LEAVE(name_+"DDz2(rccube,rccube)")
  return true;
}
  
  
  
  
  
  
// Function : Functions for some basic math operators, including
//            Gradient, Divergence and Laplacian
// Parameter: IN  - tar     > target
//            OUT - dd**    > result
// Return   : True if success
bool SpaceStepper::Gradient(rccube& target,
                            rccube& ddx1,
                            rccube& ddy1,
                            rccube& ddz1) const {
  PS_DEBUG_TRACE_ENTER(name_+"Gradient(rccube,rccube,rccube,rccube)")
  DDx1(target,ddx1);
  DDy1(target,ddy1);
  DDz1(target,ddz1);
  PS_DEBUG_TRACE_LEAVE(name_+"Gradient(rccube,rccube,rccube,rccube)")
  return true;
}
bool SpaceStepper::Gradient2(rccube& target,
                            rccube& ddx2,
                            rccube& ddy2,
                            rccube& ddz2) const {
  PS_DEBUG_TRACE_ENTER(name_+"Gradient2(rccube,rccube,rccube,rccube)")
  DDx2(target,ddx2);
  DDy2(target,ddy2);
  DDz2(target,ddz2);
  PS_DEBUG_TRACE_LEAVE(name_+"Gradient2(rccube,rccube,rccube,rccube)")
  return true;
}
bool SpaceStepper::Divergence(rccube& targetx,
                              rccube& targety,
                              rccube& targetz,
                              rccube& result) const {
  PS_DEBUG_TRACE_ENTER(name_+"Divergence(rccube,rccube,rccube,rccube)")
  for (int id = 0; id < result.GetDimension(); ++id) {
    result.GetCMatCplx(id) = cmdz1_*targetz.GetCMatCplx(id)+
                             GetCmDx1().GetCMatCplx()%targetx.GetCMatCplx(id)+
                             GetCmDy1().GetCMatCplx()%targety.GetCMatCplx(id);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"Divergence(rccube,rccube,rccube,rccube)")
  return true;
}
bool SpaceStepper::Laplacian(rccube& target, rccube& result) const {
  PS_DEBUG_TRACE_ENTER(name_+"Laplacian(rccube,rccube,rccube,rccube)")
  for (int id = 0; id < result.GetDimension(); ++id) {
    result.GetCMatCplx(id) = cmdz2_*result.GetCMatCplx(id)+
                             GetCmDx2().GetCMatCplx()%result.GetCMatCplx(id)+
                             GetCmDy2().GetCMatCplx()%result.GetCMatCplx(id);
  }
  PS_DEBUG_TRACE_LEAVE(name_+"Laplacian(rccube,rccube,rccube,rccube)")
  return true;
}
  
  
  
  
  
  
  
  
bool SpaceStepper::Gradient(rccube& target, rccube& ddxyz1, const int dimension) const {
  PS_DEBUG_TRACE_ENTER(name_+"Gradient(rccube,rccube,rccube,rccube)")
  ddxyz1.GetCubeCplx(0) = GetCmDx1().GetCubeCplx()%target.GetCubeCplx(dimension);
  ddxyz1.GetCubeCplx(1) = GetCmDy1().GetCubeCplx()%target.GetCubeCplx(dimension);
  ddxyz1.GetCMatCplx(2) = cmdz1_*target.GetCMatCplx(dimension);
  PS_DEBUG_TRACE_LEAVE(name_+"Gradient(rccube,rccube,rccube,rccube)")
  return true;
}
bool SpaceStepper::Gradient2(rccube& target, rccube& ddxyz2, const int dimension) const {
  PS_DEBUG_TRACE_ENTER(name_+"Gradient2(rccube,rccube,rccube,rccube)")
  ddxyz2.GetCubeCplx(0) = GetCmDx2().GetCubeCplx()%target.GetCubeCplx(dimension);
  ddxyz2.GetCubeCplx(1) = GetCmDy2().GetCubeCplx()%target.GetCubeCplx(dimension);
  ddxyz2.GetCMatCplx(2) = cmdz2_*target.GetCMatCplx(dimension);
  PS_DEBUG_TRACE_LEAVE(name_+"Gradient2(rccube,rccube,rccube,rccube)")
  return true;
}
bool SpaceStepper::Divergence(rccube& target, rccube& result) const {
  PS_DEBUG_TRACE_ENTER(name_+"Divergence(rccube,rccube,rccube,rccube)")
  result.GetCMatCplx(0) = cmdz1_*target.GetCMatCplx(2)+
                          GetCmDx1().GetCMatCplx()%target.GetCMatCplx(0)+
                          GetCmDy1().GetCMatCplx()%target.GetCMatCplx(1);
  PS_DEBUG_TRACE_LEAVE(name_+"Divergence(rccube,rccube,rccube,rccube)")
  return true;
}
  
  
  
  
bool SpaceStepper::GradientTmp(rccube& target, rccube& ddxyz1) const {
  PS_DEBUG_TRACE_ENTER(name_+"GradientTmp(rccube,rccube)")
  ddxyz1.GetCubeCplx(0) = GetCmDx1().GetCubeCplx()%target.GetCubeCplx(0);
  ddxyz1.GetCubeCplx(1) = GetCmDy1().GetCubeCplx()%target.GetCubeCplx(1);
  ddxyz1.GetCMatCplx(2) = cmdz1_*target.GetCMatCplx(2);
  PS_DEBUG_TRACE_LEAVE(name_+"GradientTmp(rccube,rccube)")
  return true;
}
  
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void SpaceStepper::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  
  
  stream << "Mapping Function \"";
  switch (pcfg->mappingfunction) {
    case kNone:
      stream << "None";
      break;
      
    case kTwoEndTangent:
      stream << "TwoEndTangent";
      break;
      
    case kTwoEndSin:
      stream << "TwoEndSin";
      break;
      
    case kOneEndTangent:
      stream << "OneEndTangent";
      break;
      
    default:
      stream << "Unknown";
      break;
  }
  
  stream
  << "\""
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Mapping Coefficients " << pcfg->mappingcoefficient;
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
