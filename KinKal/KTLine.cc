/*
  KTLine inherits from TLine but we also want it be an instance of KKTrk. We dont need all the parameters and functions to be redefined. KTLine follows from TLine.

  Original Author: S Middleton 2020

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
    */ 

  KTLine::KTLine(Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) :     
  KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {std::cout<<" Constructor 1 "<<speed()<<std::endl;}

  KTLine::KTLine(Vec4 const& pos0, Mom4 const& mom0, int charge, Vec3 const& bnom, TRange const& range)
  : TLine(pos0.Vect(), (mom0.Vect()/mom0.E())*CLHEP::c_light, pos0.T(), range),  trange_(range), bnom_(bnom), pos40_(pos0), mom_(mom0), charge_(charge) {
    mass_ = mom0.M();
    cout<<"Contructor 2 "<<endl;

  }

  KTLine::KTLine( PDATA const& pdata, double mass, int charge, double bnom, TRange const& range)
  : KTLine(pdata,mass,charge,Vec3(0.0,0.0,bnom),range){std::cout<<" Constructor 3 "<<speed()<<std::endl;} 

  KTLine::KTLine( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range)
  : KTLine(pdata.parameters(),pdata.covariance(),mass,charge,bnom,range) {std::cout<<" Constructor 4 "<<speed()<<std::endl;}
  
  KTLine::KTLine(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom, TRange const &trange) :  TLine(pvec, pcov), trange_(trange), bnom_(bnom), mass_(mass), charge_(charge), pars_(pvec, pcov){
  //setspeed(297.2);
  std::cout<<" Constructor 5 "<<speed()<<std::endl;

}

KTLine::KTLine(PDATA const& pdata, KTLine const& ktline) : TLine(ktline.pos40_.Vect(), (ktline.mom_.Vect()/ktline.mom_.E())*CLHEP::c_light, ktline.pos40_.T(), ktline.trange_), KTLine(pdata.parameters(),pdata.covariance(),ktline.mass_,ktline.charge_,ktline.bnom_,ktline.trange_) {std::cout<<" Constructor 6 "<<std::endl;};

  string KTLine::trajName_("KTLine");  
  string const& KTLine::trajName() { return trajName_; }

  void KTLine::momentum(double tval, Mom4& mom) const{
    cout<<" mom mag 1 "<<momentumMag(tval)<<" "<<gamma()<<" "<<mass_<<" "<<beta()<<endl;
    mom.SetPx(momentumMag(tval)*sinTheta() * sinPhi0());
    mom.SetPy(momentumMag(tval)*sinTheta() * cosPhi0());
    mom.SetPz(momentumMag(tval)*cosTheta());
    mom.SetM(mass_);

  }

  Mom4 KTLine::momentum(double tval) const{
    Mom4 mom;
    cout<<" mom mag 2 "<<momentumMag(tval)<<" "<<gamma()<<" "<<mass_<<" "<<beta()<<speed()<<endl;
    mom.SetPx(momentumMag(tval)*sinTheta() * sinPhi0());
    mom.SetPy(momentumMag(tval)*sinTheta() * cosPhi0());
    mom.SetPz(momentumMag(tval)*cosTheta());
    mom.SetM(mass_);
    return mom_;
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
      cout<<"Perpdir "<<endl;
      u.SetX(-1*cosTheta()*sinPhi0());
      u.SetY(-1*cosTheta()*cosPhi0());
      u.SetZ(sinTheta());
      cout<<" Unit in perp "<<u<<endl;
      return u;
    break;
      cout<<"phi dir "<<endl;
      case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      u.SetX(-cosPhi0());
      u.SetY(sinPhi0());
      u.SetZ(0.0);
      cout<<" Unit in phi "<<u<<endl;
      return u;
    break;
      case LocalBasis::momdir: // along momentum: check.
      u.SetX(sinTheta() * sinPhi0());
      u.SetY(sinTheta() * cosPhi0());
      u.SetZ(cosTheta());
      cout<<" Unit in mom "<<u<<endl;
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
	      pder[z0_] = (1*l/sinTheta());
	      pder[t0_] = dt;
        cout<<" deriv perpdir "<<pder <<endl;
	      break;
      case LocalBasis::phidir:

	      // change in phi0*costheta
	      pder[cost_] = 0;
	      pder[d0_] = l/sinTheta();
	      pder[phi0_] = -1/sinTheta();
	      pder[z0_] = d0()*(1/sinTheta()*tanTheta());
	      pder[t0_] = dt;
        cout<<"deriv phi dir "<<pder<<endl;
	      break;
      case LocalBasis::momdir:
        cout<<"mom dir "<<u<<endl;
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
      ost << KTLine::paramName(static_cast<KTLine::ParamIndex>(ipar) ) << " : " << paramVal(ipar);
      if(ipar < KTLine::npars_-1) ost << " , ";
    }
    ost << endl;
  }


} // KinKal namespace
