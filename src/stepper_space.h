/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVONE_SPACE_STEPPER_H_
#define YLY_PS_LVONE_SPACE_STEPPER_H_

#include "base.h"
#include "domain.h"
#include "file_ini.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
  
class SpaceStepper : public levelbase::Base {
  
public:
  
  
  enum MappingFunction {
    kNone          = 0,
    kTwoEndTangent = 1,     // z = 2/pi*arctan[a*tan(pi/2*zeta)]
    kTwoEndSin     = 2,     // 
    kOneEndTangent = 3,     // z = 4/pi*arctan[a*tan(pi/4*(zeta-1))]+1
  }; typedef MappingFunction mfun_t;
  
  
  
  struct Configuration : public levelbase::Base::Configuration {
  public:
    mfun_t mappingfunction = kNone;
    dvec mappingcoefficient;
    bool Load(levelbase::IniFile& fileini);
  protected:
  private:
  }; typedef Configuration cfg_t;
  
  
  
  SpaceStepper(const levelbase::mpi_t& mpi,
               const levelbase::dcp_t& ddbt,
               Domain* domain,
               const levelbase::Base::cfg_t* cfg = NULL);
  virtual ~SpaceStepper(void) { }
  
  
  Domain* GetDomain(void) { return domain_; }
  
  const levelbase::rccube& GetCmDx1 (void) const { return cmdx1_ ; }
  const levelbase::rccube& GetCmDy1 (void) const { return cmdy1_ ; }
  const levelbase::rccube& GetCmDxy1(void) const { return cmdxy1_; }
  
  const levelbase::rccube& GetCmDx2 (void) const { return cmdx2_ ; }
  const levelbase::rccube& GetCmDy2 (void) const { return cmdy2_ ; }
  const levelbase::rccube& GetCmDxy2(void) const { return cmdxy2_; }
  
  const cmat& GetCmDz1(void) const { return cmdz1_ ; }
  const cmat& GetCmDz2(void) const { return cmdz2_ ; }

  
  
  // Discretization of partial derivative calculation.
  bool DDx1(levelbase::rccube& target, levelbase::rccube& ddx1) const;
  bool DDx2(levelbase::rccube& target, levelbase::rccube& ddx2) const;
  bool DDy1(levelbase::rccube& target, levelbase::rccube& ddy1) const;
  bool DDy2(levelbase::rccube& target, levelbase::rccube& ddy2) const;
  bool DDz1(levelbase::rccube& target, levelbase::rccube& ddz1) const;
  bool DDz2(levelbase::rccube& target, levelbase::rccube& ddz2) const;
  
  
  
  // Discretization of general mathematic operators.
  bool Gradient  (levelbase::rccube& target ,
                  levelbase::rccube& ddx1,
                  levelbase::rccube& ddy1,
                  levelbase::rccube& ddz1) const;
  bool Gradient2 (levelbase::rccube& target ,
                  levelbase::rccube& ddx2,
                  levelbase::rccube& ddy2,
                  levelbase::rccube& ddz2) const;
  bool Divergence(levelbase::rccube& targetx,
                  levelbase::rccube& targety,
                  levelbase::rccube& targetz,
                  levelbase::rccube& result ) const;
  bool Laplacian (levelbase::rccube& target ,
                  levelbase::rccube& result ) const;
  
  
  
  bool Gradient  (levelbase::rccube& target,
                  levelbase::rccube& ddxyz1,
                  const int dimension = 0) const;
  bool Gradient2 (levelbase::rccube& target,
                  levelbase::rccube& ddxyz2,
                  const int dimension = 0) const;
  bool Divergence(levelbase::rccube& target,
                  levelbase::rccube& result) const;
  
  
  
  bool GradientTmp(levelbase::rccube& target,
                   levelbase::rccube& ddxyz1) const;
  
  
  
  bool Update(void);
  virtual bool Configure(void);
  virtual void Print(ostringstream& stream) const;
  
protected:
  
  Domain* domain_ = NULL;
  
  
  levelbase::rccube cmdx1_;
  levelbase::rccube cmdy1_;
  levelbase::rccube cmdxy1_;
  
  levelbase::rccube cmdx2_;
  levelbase::rccube cmdy2_;
  levelbase::rccube cmdxy2_;
  
  
  cmat cmdz1_;
  cmat cmdz2_;
  cmat tmz1_;
  cmat tmz2_;
  
  
  ccube cmdxxx1_;
  ccube cmdyyy1_;
  
  ccube cmdxxx2_;
  ccube cmdyyy2_;
  
  
  
  
  string name_ = "SpaceStepper::";

private:
  
};
  
  
  
struct TypeSpaceStepper {
  SpaceStepper::cfg_t* cfg = NULL;
  SpaceStepper       * obj = NULL;
}; typedef TypeSpaceStepper sstepper_t;
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVONE_SPACE_STEPPER_H_
