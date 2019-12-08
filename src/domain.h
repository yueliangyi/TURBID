/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_DOMAIN_H_
#define YLY_PS_LVONE_DOMAIN_H_

#include "base.h"
#include "cube_real_cplx.h"
#include "file_ini.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class Domain : public levelbase::Base {
  
public:
  
  struct Configuration : public levelbase::Base::Configuration {
  public:
    double lengthx = 0;
    double lengthy = 0;
    double lengthz = 0;
    double slope   = 0;
    int    slope_direction = 1; // default x
    bool   readprofile = false;
    double rippleamplitude = 0;
    
    bool Load(levelbase::IniFile& fileini);
    
  protected:
  private:
  }; typedef Configuration cfg_t;
  
  Domain(const levelbase::mpi_t& mpi,
         const levelbase::dcp_t& dcp,
         const levelbase::Base::cfg_t* cfg);
  virtual ~Domain(void);
  
  double GetLengthX(void) const { return ((cfg_t*)cfg_)->lengthx; }
  double GetLengthY(void) const { return ((cfg_t*)cfg_)->lengthy; }
  double GetLengthZ(void) const { return ((cfg_t*)cfg_)->lengthz; }
  
  double GetLengthMinX(void) const { return lengthminx_; }
  double GetLengthMinY(void) const { return lengthminy_; }
  double GetLengthMinZ(void) const { return lengthminz_; }
  
  const dvec& GetGridX(void) const { return gridx_; }
  const dvec& GetGridY(void) const { return gridy_; }
  const dvec& GetGridZ(void) const { return gridz_; }
  dcube& GetMaxULen(void) { return maxulen_; }
  const dcube& GetGridLenX(void) const { return gridlenx_; }
  const dcube& GetGridLenY(void) const { return gridleny_; }
  const dcube& GetGridLenZ(void) const { return gridlenz_; }
  const dmat& GetProfileTop   (void) const { return profiletop_   ; }
  const dmat& GetProfileBottom(void) const { return profilebottom_; }
  
  const levelbase::rccube& GetEtaM(const int dimension = 0) const { return *etam_[dimension]; }
  const levelbase::rccube& GetEtaR(const int dimension = 0) const { return *etar_[dimension]; }
  const levelbase::rccube& GetJm3K(const int dimension = 0) const { return *jm3k_[dimension]; }
  const levelbase::rccube& GetPhyGridX1(const int dimension = 0) const { return *phy_grid_x1_[dimension]; }
  const levelbase::rccube& GetPhyGridX2(const int dimension = 0) const { return *phy_grid_x2_[dimension]; }
  const levelbase::rccube& GetPhyGridX3(const int dimension = 0) const { return *phy_grid_x3_[dimension]; }
  const levelbase::rccube& GetDVecNormal(const int dimension = 0) const { return *dvec_normal_[dimension]; }
  const levelbase::rccube& GetDVecTangent1(const int dimension = 0) const { return *dvec_tangent1_[dimension]; }
  const levelbase::rccube& GetDVecTangent2(const int dimension = 0) const { return *dvec_tangent2_[dimension]; }

  dvec& GetGridZCluster(void) { return gridzcluster_; }
  dcube& GetGridLenZ(void) { return gridlenz_; }
  levelbase::rccube& GetEtaM(const int dimension = 0) { return *etam_[dimension]; }
  levelbase::rccube& GetEtaR(const int dimension = 0) { return *etar_[dimension]; }
  levelbase::rccube& GetJm3K(const int dimension = 0) { return *jm3k_[dimension]; }
  levelbase::rccube& GetPhyGridX1(const int dimension = 0) { return *phy_grid_x1_[dimension]; }
  levelbase::rccube& GetPhyGridX2(const int dimension = 0) { return *phy_grid_x2_[dimension]; }
  levelbase::rccube& GetPhyGridX3(const int dimension = 0) { return *phy_grid_x3_[dimension]; }
  levelbase::rccube& GetDVecNormal(const int dimension = 0) { return *dvec_normal_[dimension]; }
  levelbase::rccube& GetDVecTangent1(const int dimension = 0) { return *dvec_tangent1_[dimension]; }
  levelbase::rccube& GetDVecTangent2(const int dimension = 0) { return *dvec_tangent2_[dimension]; }

  virtual bool Update(void);
  virtual bool Configure(void);
  virtual void Print(ostringstream& stream) const;
  
protected:
  
  // Minimum length of grid.
  double lengthminx_  = 0;
  double lengthminy_  = 0;
  double lengthminz_  = 0;
  
  dvec gridx_;
  dvec gridy_;
  dvec gridz_;
  dvec gridzcluster_;
  
  dcube maxulen_;
  dcube gridlenx_;
  dcube gridleny_;
  dcube gridlenz_;

  dmat profiletop_;
  dmat profilebottom_;
  
//  levelbase::rccube physicalgrid_;
  
  std::vector<levelbase::rccube*> etam_;
  std::vector<levelbase::rccube*> etar_;
  std::vector<levelbase::rccube*> jm3k_;
  
  std::vector<levelbase::rccube*> phy_grid_x1_;
  std::vector<levelbase::rccube*> phy_grid_x2_;
  std::vector<levelbase::rccube*> phy_grid_x3_;
  std::vector<levelbase::rccube*> dvec_normal_;
  std::vector<levelbase::rccube*> dvec_tangent1_;
  std::vector<levelbase::rccube*> dvec_tangent2_;
  
  string name_ = "Domain::";

private:
  
};
  
  
  
struct TypeDomain {
  Domain::cfg_t* cfg = NULL;
  Domain       * obj = NULL;
}; typedef TypeDomain domain_t;
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_DOMAIN_H_
