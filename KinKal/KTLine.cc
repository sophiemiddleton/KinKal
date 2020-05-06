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
    /*
    KTLine can take in Momentum externally as a 4-vector or calculate it based. You can initialize the line with an origin (pos0) or the trajectory parameters (pdata)
    To make things work it is essential that:
      * mass - set in KInter and in the KTLine constructors
      * speed - set in TLine and in KTline constructors
      * direction - set in TLine and in KTLine constructors

      are all set.

    */ 

    /* ...In thsis contructor mom4 is passed in, speed is passed to TLine */
    KTLine::KTLine(Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) :     
KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}

  KTLine::KTLine(Vec4 const& pos0, Mom4 const& mom0, int charge, Vec3 const& bnom, TRange const& range)
  : KInter(mom0.M(),charge),  TLine(Vec3(pos0.x(),pos0.y(),pos0.z()), Vec3(mom0.Px()/mass_,mom0.Py()/mass_,mom0.Pz()/mass_), pos0.T(), range),  trange_(range), bnom_(bnom), mom_(mom0) {
   /* Vec3 mom3vec;
    mom3vec.SetXYZ(mom0.Px()/mass_,mom0.Py()/mass_,mom0.Pz()/mass_);
    speed_(mom3vec.mag());*/
  }

   /*...In this instance mass and speed are passed in */
  KTLine::KTLine( PDATA const& pdata, double mass, int charge, double bnom, TRange const& range, double speed)
  : KTLine(pdata,mass,charge,Vec3(0.0,0.0,bnom),range,speed){} //speed passed as 

  KTLine::KTLine( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range, double speed)
  : KInter(mass, charge), TLine(pdata), trange_(range), pars_(pdata), bnom_(bnom), mass_(mass)){} //speed

  void KTLine::momentum(double tval, Mom4& mom) const{
   mom.SetPx(momentumMag(tval)*dir().x());
   mom.SetPy(momentumMag(tval)*dir().y());
   mom.SetPz(momentumMag(tval)*dir().z());
   mom.SetM(mass_);
  }

/*

The effects for changes in 2 perpendicular directions (theta1 = theta and
theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
are uncorrelated. These axes are track specific. as cosmics are not always
coming along the same track direction it is necessary to have difference
parameterization than that used for the helix case.

*/
  Vec3  KTLine::direction(double t, LocalBasis::LocDir mdir) const {
    Vec3 u;
    switch ( mdir ) {
      case LocalBasis::perpdir: // purely polar change theta 1 = theta
	      u.SetX(-1*cosTheta()*sinPhi0());
	      u.SetY(-1*sinTheta()*sinPhi0());
	      u.SetZ(sinTheta());
	      return u;
	    break;
        case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
	        u.SetX(-cosPhi0());
          u.SetY(sinPhi0());
          u.SetZ(0.0);
          return u;
	    break;
        case LocalBasis::momdir: // along momentum: check.
	       u.SetX(mom().Px()/mom().mag());
         u.SetY(mom().Py()/mom().mag());
         u.SetZ(mom().Pz()/mom().mag());
         return u;
	    break;
        default:
	        throw std::invalid_argument("Invalid direction");
    }
    if(needsrot_) {
      u = brot_(u); 
      return u;
    }
  }

// derivatives of momentum projected along the given basis WRT the 5 parameters
  void KTLine::momDeriv(double t, LocalBasis::LocDir mdir, DVEC &pder, Vec3& u) const{
    // compute some useful quantities
    double dt = t-t0();
    double l = CLHEP::c_light * beta() * (dt);
    u.SetX(sinTheta()*sinPhi0());
    u.SetY(sinTheta()*cosPhi0());
    u.SetZ(cosTheta());
    // cases
    switch ( mdir ) {
      case LocalBasis::perpdir:
	      // polar bending: only momentum and position are unchanged
	      pder[cost_] = 1;
	      pder[d0_] = 0;
	      pder[phi0_] = 0;
	      pder[z0_] = (-1*l/sinTheta());
	      pder[t0_] = dt;

	      break;
      case LocalBasis::phidir:
	      // change in phi0*costheta
	      pder[cost_] = 0;
	      pder[d0_] = l/sinTheta();
	      pder[phi0_] = -1/sinTheta();
	      pder[z0_] = d0()*(1/sinTheta()*tanTheta());
	      pder[t0_] = dt;

	      break;
      case LocalBasis::momdir:
	      // fractional momentum change: position and direction are unchanged
	      u = direction(t);
	    break;
          default:
	    throw std::invalid_argument("Invalid direction");
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
