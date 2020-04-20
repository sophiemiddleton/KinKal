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

    // Transform into the system where Z is along the Bfield.
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm; signed by the charge!!!
    double momToRad = 1000.0/(charge_*bnom_.R()*CLHEP::c_light);
    // reduced mass; note sign convention!
    mbar_ = -mass_*momToRad;
    // transverse radius of the helix
    //TODO - what do we need in Line class here?
  }

  KTLine::KTLine PDATA const& pdata, double mass, int charge, double bnom, TRange const& range) : KTLine(pdata,mass,charge,Vec3(0.0,0.0,bnom),range) {}
  KTLine::KTLine( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range) : 
    TTraj(range), KInter(mass,charge), pars_(pdata), bnom_(bnom) {}

// KTLine inherits from KInter but we also want it be an instance of KKTrk, we need:
//      void position(Vec4& pos) const; -->TLine
//      void position(float time, Vec3& pos) const; -->TLine
//      void velocity(float time, Vec3& vel) const; -->TLine
//      double speed(float time) const; --> TLine
//      void direction(float time, Vec3& dir) const; --> TLine
//      void print(std::ostream& ost, int detail) const; --> TLine but Override
// These are thing things we need to intiate in KTLine:
//      void momentum(double t,Mom4& mom) const; 
//      void momentum(Vec4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); }
//      double momentum(float time) const; // momentum and energy magnitude in MeV/
//      double momentumVar(float time) const; // variance on momentum value
//      double energy(float time) const; 
//      void rangeInTolerance(TRange& range, BField const& bfield, double tol);
//      PDATA const& params() const;
// Many of these are satisifed already within the TLine class.

  double KTLine::momentumVar(float time) const {
    //TODO
    PDATA::DVEC dMomdP(0.0,  0.0, 0.0 ,0.0 , 0.0);
    dMomdP *= mass()/(pbar()*mbar());
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  void KTLine::momentum(double tval, Mom4& mom) const{
   mom.SetPx(mom()*dir().x());
   mom.SetPy(mom()*dir().y());
   mom.SetPz(mom()*dir().z());
   mom.SetM(mass_);
  }

 void KTLine::velocity(double tval,Vec3& vel) const{//TODO - do we need thins?
    Mom4 mom;
    momentum(tval,mom);
    vel = mom.Vect()*(CLHEP::c_light*fabs(Q()/ebar()));
    if(needsrot_)vel = brot_(vel);
  }

/*  The effects for changes in 2 perpendicular directions (theta1 = theta and
  theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
  are uncorrelated. These axes are track specific. as cosmics are not always coming along the same track direction it is necessary to have difference parameterization than that used for the helixa case. 

*/
  void KTLine::dirVector(MDir mdir,double tval,Vec3& unit) const {
    double phival = phi(time); // azimuth at this point
    double norm = 1.0/copysign(pbar(),mbar_); // sign matters!
    switch ( mdir ) {
      case theta1: // purely plar change theta 1 = theta
	      unit.SetX();
	      unit.SetY();
	      unit.SetZ();
	      unit *= norm;
	    break;
        case theta2: // purely transverse - theta2 = phi()*sin(theta)
	        unit.SetX();
	        unit.SetY();
	        unit.SetZ(0.0);
	    break;
        case momdir: // along momentum: sign matters!
	        direction(time,unit);
	    break;
        default:
	        throw std::invalid_argument("Invalid direction");
    }
    if(needsrot_) unit = brot_(unit);
  }

// derivatives of momentum projected along the given basis WRT the 5 parameters
  void KTLine::momDeriv(MDir mdir, double time, PDER& pder) const {
   //TODO
    // compute some useful quantities
    double dt = time-t0();
    double fltlen = ...*dt;
    // cases
    switch ( mdir ) {
      case theta1:
	      // polar bending: only momentum and position are unchanged
	      pder[cost_] = 1;
	      pder[d0_] = 0;
	      pder[phi0_] = 0;
	      pder[z0_] = -1*fltlen*(1/sintheta));
	      pder[t0_] = dt;
	      
	      break;
      case theta2:
	      // Azimuthal bending: R, Lambda, t0 are unchanged
	      pder[cost_] = 0;
	      pder[d0_] = -1*fltlen;
	      pder[phi0_] = 1/sintheta;
	      pder[z0_] = -d0()*(1/sintheta*tantheta);
	      pder[t0_] = dt;
    
	      break;
      case momdir:
	      // fractional momentum change: position and direction are unchanged
	      pder[cost_] = 0;
	      pder[d0_] = 0;
	      pder[phi0_] = 0;
	      pder[z0_] = 0;
	      pder[t0_] = 0;
	    break;
          default:
	    throw std::invalid_argument("Invalid direction");
        }
  }
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
 
    void KTLine::print(std::ostream& ost, int detail) const {
    ost << " KTLine " <<  range() << " parameters: ";
    for(size_t ipar=0;ipar < KTLine::npars_;ipar++){
      ost << KTLine::paramName(static_cast<KTLine::ParamIndex>(ipar) ) << " : " << param(ipar);
      if(ipar < KTLine::npars_-1) ost << " , ";
    }
    ost << endl;
  }


} // KinKal namespace
