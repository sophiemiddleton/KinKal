#include "KinKal/KTLine.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  //Note : This class inherits a lot from TLine. So, I think we dont need all that LHelix has here as its defined in TLine
  vector<string> KTLine::paramUnits_ = {
  "mm","radians","mm","","ns"};//d0, phi0,z0, cost
  
  std::vector<std::string> const& KTLine::paramUnits() { return paramUnits_; }
  
  std::string const& KTLine::paramUnit(ParamIndex index) { return paramUnits_[static_cast<size_t>(index)];}
  

  KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) : KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}
  KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, int charge, Vec3 const& bnom, TRange const& range) : TTraj(range), KInter(mom0.M(),charge), bnom_(bnom), needsrot_(false) {
    static double twopi = 2*M_PI; // FIXME
    // Transform into the system where Z is along the Bfield.
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    if(fabs(bnom_.Theta()) >1.0e-6){ //TODO - understand this.
      needsrot_ = true;
      Rotation3D rot(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
      pos = rot(pos);
      mom = rot(mom);
      // create inverse rotation
     brot_ = rot.Inverse();
      // check
      auto test = rot(bnom_);
      if(fabs(test.Theta()) > 1.0e-6)throw std::invalid_argument("BField Error");
    }
    
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm; signed by the charge!!!
    double momToRad = 1000.0/(charge_*bnom_.R()*CLHEP::c_light);
    // reduced mass; note sign convention!
    mbar_ = -mass_*momToRad;
    // transverse radius of the helix
    //TODO - in Line
  }

  KTLine::KTLine PDATA const& pdata, double mass, int charge, double bnom, TRange const& range) : KTLine(pdata,mass,charge,Vec3(0.0,0.0,bnom),range) {}
  KTLine::KTLine( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range) : 
    TTraj(range), KInter(mass,charge), pars_(pdata), bnom_(bnom) {}

  void KTLine::position(Vec4& pos) const {
    Vec3 temp;
    position(pos.T(),temp);
    pos.SetXYZT(temp.X(),temp.Y(),temp.Z(),pos.T());
  }
  
  void KTLine::position(double t, Vec3& pos) const {
   
   // pos.SetX();
   // pos.SetY();
   // pos.SetZ();
    if(needsrot_) pos = brot_(pos);
 } 

  void KTLine::momentum(double tval,Mom4& mom) const{
 
   // mom.SetPx();
   // mom.SetPy();
   // mom.SetPz();
   // mom.SetM(mass_);
   
  }

 void KTLine::velocity(double tval,Vec3& vel) const{
    Mom4 mom;
    momentum(tval,mom);
    vel = mom.Vect()*(CLHEP::c_light*fabs(Q()/ebar()));
    if(needsrot_)vel = brot_(vel);
  }

  void KTLine::direction(double tval,Vec3& dir) const{
    Mom4 mom;
    momentum(tval,mom);
    dir = mom.Vect().Unit();
    if(needsrot_)dir = brot_(dir);
  }

  void KTLine::dirVector(MDir mdir,double tval,Vec3& unit) const {
    //TODO
  }

// derivatives of momentum projected along the given basis WRT the 6 parameters
  void KTLine::momDeriv(MDir mdir, double time, PDER& pder) const {
    // compute some useful quantities
  //TODO
  }

  // derivatives of position.Dot(direction) WRT the 6 parameters
  // these are used to apply the continuity constraint at lossy effects
  void KTLine::posDeriv(double time, PDER& pder) const {
   //TODO

  }

  void KTLine::rangeInTolerance(TRange& brange, BField const& bfield, double tol) const {
    // precompute some factors
    double fact = 0.5*sqrt(rad()*tol*bnom().R())/CLHEP::c_light;
    // Limit to this traj's range
    brange.high() = std::min(brange.high(),range().high());
    // compute the BField difference in the middle of the range
    Vec3 midpos,bvec;
    position(brange.mid(),midpos);
    bfield.fieldVect(bvec,midpos);
    auto db = bvec-bnom();
    double dt = fact/sqrt(db.R());
    // truncate the range if necessary
    if(dt < brange.range())brange.high() = brange.low() + dt;
  }
 

  std::ostream& operator <<(std::ostream& ost, KTLine const& ktline) {
    ost << " KTLine parameters: ";
    for(size_t ipar=0;ipar < KTLine::npars_;ipar++){
      ost << KTLine::paramName(static_cast<KTLine::ParamIndex>(ipar) ) << " " << ktline.param(ipar);
      if(ipar < KTLine::npars_-1) ost << " ";
    }
    ost << ktline.range();
    return ost;
  }

} // KinKal namespace
