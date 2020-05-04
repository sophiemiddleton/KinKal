/*
KTLine inherits from KInter but we also want it be an instance of KKTrk, for that we need:

      void position(Vec4& pos) const; -->TLine
      void position(float time, Vec3& pos) const; -->TLine
      void velocity(float time, Vec3& vel) const; -->TLine
      double speed(float time) const; --> TLine
      void direction(float time, Vec3& dir) const; --> TLine
      void print(std::ostream& ost, int detail) const; --> TLine but Override

 These are thing things we need to intiate in KTLine:
      void momentum(double t,Mom4& mom) const; -->KTLine
      void momentum(Vec4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); } -->no override
      double momentum(float time) const; // momentum and energy magnitude in MeV/
      double momentumVar(float time) const; // variance on momentum value --> this KTLine class
      double energy(float time) const; -->KTLine.hh
      void rangeInTolerance(TRange& range, BField const& bfield, double tol);-->this KTLine class
      PDATA const& params() const;

    s Middleton 2020

*/

#include "KinKal/KTLine.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  //Note : This class inherits a lot from TLine. So, I think we dont need all that LHelix has here as its defined in TLine
  KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, double mass, int charge, double bnom, TRange const& range)
  : KTLine(pos0, mom0, mass, mass, charge, Vec3(0.0,0.0,bnom), range) {
  }
  KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, double mass, int charge, Vec3 const& bnom, TRange const& range)
  : TTraj(range), KInter(mom0.M(),charge), bnom_(bnom), mass_(mass), needsrot_(false) {
    // Transform into the system where Z is along the Bfield.
    Vec4 pos(pos0);
    Mom4 mom(mom0);
  }
  KTLine::KTLine( PDATA const& pdata, double mass, int charge, double bnom, TRange const& range)
  : KTLine(pdata,mass,charge,Vec3(0.0,0.0,bnom),range) {}


  KTLine::KTLine PDATA const& pdata, double mass, int charge, double bnom, TRange const& range)
  : KTLine(pdata, mass, charge, Vec3( 0.0, 0.0, bnom), range) {}
  KTLine::KTLine( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range)
  :TTraj(range), KInter(mass, charge), pars_(pdata), bnom_(bnom) {}

  void KTLine::momentum(double tval, Mom4& mom) const{
   mom.SetPx(momentum(tval) *dir().x());
   mom.SetPy(momentum(tval)  *dir().y());
   mom.SetPz(momentum(tval)  *dir().z());
   mom.SetM(mass_);
  }

 void KTLine::velocity(double tval,Vec3& vel) const{
    vel = dir()*speed();
  }

/*

The effects for changes in 2 perpendicular directions (theta1 = theta and
theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
are uncorrelated. These axes are track specific. as cosmics are not always
coming along the same track direction it is necessary to have difference
parameterization than that used for the helix case.

*/
  void KTLine::dirVector(MDir mdir,double tval,Vec3& unit) const {
    switch ( mdir ) {
      case theta1: // purely plar change theta 1 = theta
	      unit.SetX(-1*cosTheta()*sinPhi());
	      unit.SetY(-1*sinTheta()*sinPhi());
	      unit.SetZ(sinTheta());
	      unit *= norm;
	    break;
        case theta2: // purely transverse theta2 = -phi()*sin(theta)
	        unit.SetX(-cosPhi());
          unit.SetY(sinPhi());
          unit.SetZ(0.0);
	    break;
        case momdir: // along momentum: check.
	        direction(time,unit);
	    break;
        default:
	        throw std::invalid_argument("Invalid direction");
    }
    if(needsrot_) unit = brot_(unit); //TODO - what is this used for?
  }

// derivatives of momentum projected along the given basis WRT the 5 parameters
  void KTLine::momDeriv(MDir mdir, double time, PDER& pder) const {
    // compute some useful quantities
    double dt = time-t0();
    double l = CLHEP::c_light * beta() * (dt);
    // cases
    switch ( mdir ) {
      case theta1:
	      // polar bending: only momentum and position are unchanged
	      pder[cost_] = 1;
	      pder[d0_] = 0;
	      pder[phi0_] = 0;
	      pder[z0_] = (-1*l/sinTheta());
	      pder[t0_] = dt;

	      break;
      case theta2:
	      // change in phi0*costheta
	      pder[cost_] = 0;
	      pder[d0_] = l/sinTheta();
	      pder[phi0_] = -1/sinTheta();
	      pder[z0_] = d0()*(1/sinTheta()*tanTheta());
	      pder[t0_] = dt;

	      break;
      case momdir:
	      // fractional momentum change: position and direction are unchanged
	      direction(time,unit);
	    break;
          default:
	    throw std::invalid_argument("Invalid direction");
        }
  }
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
