/*
  KTLine is the Linear Trajectory Specialization of KTRAJ - the kinematic trajectory.
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
  vector<string> KTLine::paramTitles_ = {
      "Transverse DOCA to Z Axis (d_{0})", "Azimuth of POCA (#phi_{0})",
      "Z at POCA (z_{0})", "Tan #lambda", "Time at POCA (t_{0})"};

  vector<string> KTLine::paramNames_ = {"d_{0}", "#phi_{0}", "z_{0}",
                                       "Tan #lambda", "t_{0}"};

  vector<string> KTLine::paramUnits_ = {"mm", "radians", "mm", "", "ns"};

  std::vector<std::string> const &KTLine::paramUnits() { return paramUnits_; }
  std::vector<std::string> const &KTLine::paramNames() { return paramNames_; }
  std::vector<std::string> const &KTLine::paramTitles() { return paramTitles_; }

  std::string const &KTLine::paramName(ParamIndex index) {
    return paramNames_[static_cast<size_t>(index)];
  }
  std::string const &KTLine::paramTitle(ParamIndex index) {
    return paramTitles_[static_cast<size_t>(index)];
  }
  std::string const &KTLine::paramUnit(ParamIndex index) {
    return paramUnits_[static_cast<size_t>(index)];
  }

  string KTLine::trajName_("KTLine");
  string const &KTLine::trajName() { return trajName_; }

  KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) : KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}

  KTLine::KTLine(Vec4 const &pos0, Mom4 const &mom0, int charge, Vec3 const &bnom,
  TRange const &trange) :  bnom_(bnom), mass_(mom0.M()), charge_(charge), trange_(trange)
  {

    pos40_ = pos0;
    mom_ = mom0;
    speed_ = (sqrt(((mom0.Vect() / mom0.E()) * CLHEP::c_light).Mag2()));
    dir_ = ((mom0.Vect() / mom0.E()) * CLHEP::c_light).Unit();

    //static const Vec3 zdir(0.0, 0.0, 1.0);
    //double zddot = zdir.Dot(dir_);
   // double theta = acos(zddot);
    param(t0_) = pos0.T() - (pos0.Z() - param(z0_)) / (cosTheta() * CLHEP::c_light * beta());

 // Transform into the system where Z is along the Bfield.  This is a pure rotation about the origin
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    if(fabs(g2l_(bnom_).Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    pos = g2l_(pos);
    mom = g2l_(mom);
    // create inverse rotation; this moves back into the original coordinate system
    l2g_ = g2l_.Inverse();
    double momToRad = 1.0/(BField::cbar()*charge_*bnom_.R());
    mbar_ = -mass_ * momToRad;

    double pt = sqrt(mom.perp2());
    double radius = fabs(pt*momToRad);

    double lambda = -mom.z()*momToRad;//M_PI_2 - theta;//
    double amsign = copysign(1.0, mbar_);

    Vec3 center = Vec3(pos.x() + mom.y()*momToRad, pos.y() - mom.x()*momToRad, 0.0);
    double rcent = sqrt(center.perp2());
    double fcent = center.phi();
    double centerx = rcent*cos(fcent);
    double centery = rcent*sin(fcent);

    param(tanl_) = amsign*(lambda)/radius;
cout<<amsign<<" * "<<rcent<<" - "<<radius<<" = "<<amsign*(rcent - radius)<<std::endl;
    param(d0_) = amsign*(rcent - radius);
    param(phi0_) = atan2(-amsign * centerx, amsign * centery);

    Vec3 pos3 = Vec3(pos.x(), pos.y(), pos.z());
    double fz0 = (pos3 - center).phi() - pos3.z() / lambda;
    deltaPhi(fz0); //TODO - delta Phi
    double refphi = fz0+amsign*M_PI_2;
    double phival = phi0();
    double dphi = deltaPhi(phival, refphi);
    param(z0_) = dphi * cosTheta(); //TODO  pos.Z();
    param(t0_) = pos.T() - (pos.Z() - param(z0_)) / (sinDip() * CLHEP::c_light * beta());
    cout << "In KTLine. Params set to: " << pars_.parameters() << endl;

    vt_ = CLHEP::c_light * pt / mom.E();
    vz_ = CLHEP::c_light * mom.z() / mom.E();
    /*// test position and momentum function
    Vec4 testpos(pos0);
    // std::cout << "Testpos " << testpos << std::endl;
    position(testpos);
    Mom4 testmom = momentum(testpos.T());
    auto dp = testpos.Vect() - pos0.Vect();
    auto dm = testmom.Vect() - mom0.Vect();
    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5) throw invalid_argument("Rotation Error");*/
  }

  double KTLine::deltaPhi(double &phi, double refphi) const
  {
    double dphi = phi - refphi;
    static const double twopi = 2 * M_PI;
    while (dphi > M_PI)
    {
      dphi -= twopi;
      phi -= twopi;
    }
    while (dphi <= -M_PI)
    {
      dphi += twopi;
      phi += twopi;
    }
    return dphi;
  }

  KTLine::KTLine(PDATA const &pdata, KTLine const &other) : KTLine(other) {
    pars_ = pdata;
  }

  void KTLine::position(Vec4 &pos) const {
    Vec3 pos3 = position(pos.T());
    pos.SetXYZT(pos3.X(), pos3.Y(), pos3.Z(), pos.T());
  }

  /*Vec3 KTLine::position(double time) const {
    if (forcerange_){
      range().forceRange(time);
    }
    return (pos0() + ((time - t0()) * speed()) * dir());
  }*/

 Vec3 KTLine::position(double time) const
  {
    double cDip = cosDip();
    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
    double sphi0 = sin(phi0());
    double cphi0 = cos(phi0());
    return l2g_(Vec3(-d0()*sphi0, d0()*cphi0, z0() + l * tanl()));
  }

  Vec4 KTLine::pos4(double time) const {
    Vec3 temp = position(time);
    return Vec4(temp.X(), temp.Y(), temp.Z(), time);
  }

  void KTLine::momentum(double tval, Mom4 &mom) const {
    Vec3 dir = direction(tval);
    mom.SetPx(momentumMag(tval) * dir.x());
    mom.SetPy(momentumMag(tval) * dir.y());
    mom.SetPz(momentumMag(tval) * dir.z());
    mom.SetM(mass_);
  }

  Mom4 KTLine::momentum(double tval) const {
    Mom4 mom;
    Vec3 dir = direction(tval);
    mom.SetPx(momentumMag(tval) * dir.x());
    mom.SetPy(momentumMag(tval) * dir.y());
    mom.SetPz(momentumMag(tval) * dir.z());
    mom.SetM(mass_);
    return mom_;
  }

  /*
  The effects for changes in 2 perpendicular directions (theta1 = theta and
  theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
  are uncorrelated. These axes are track specific. as cosmics are not always
  coming along the same track direction it is necessary to have difference
  parameterization than that used for the helix case.
  alt dir = a test with the "BTrk parameterization" - just changes signs due to
  swithc in cos<->sin
  */
  Vec3  KTLine::direction(double t, LocalBasis::LocDir mdir) const {

    switch ( mdir ) {
    case LocalBasis::perpdir: // purely polar change theta 1 = theta
      //return l2g_(Vec3(cosTheta()*sinPhi0(),cosTheta()*cosPhi0(),-1*sinTheta()));
      return l2g_(Vec3(-1*sinDip()*cosPhi0(),-1*sinDip()*sinPhi0(),cosTheta()));
    break;
      case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      //return l2g_(Vec3(-cosPhi0(),sinPhi0(),0.0));
      return l2g_(Vec3(-sinPhi0(),cosPhi0(),0.0));
    break;
      case LocalBasis::momdir: // along momentum: check.
      //return l2g_(Vec3(sinPhi0()*sinTheta(),cosPhi0()*sinTheta(), cosTheta())); 
      return l2g_(Vec3(cosPhi0(),sinPhi0(), tanl())); 
    break;
      default:
      throw std::invalid_argument("Invalid direction");
    }
  }


  // derivatives of momentum projected along the given basis WRT the 5 parameters
  KTLine::DVEC KTLine::momDeriv(double time, LocalBasis::LocDir mdir) const {
 // compute some useful quantities
   double tanval = cosTheta()/sinTheta();
    double cosval = sinTheta();
/*    double l = translen(CLHEP::c_light * beta() * (time - t0()));
    double d0val = d0();

  DVEC pder;
    // cases
    switch ( mdir ) {
      case LocalBasis::perpdir:
        // polar bending: only momentum and position are unchanged
        pder[d0_] = d0val;//tanval*(1-cos(omval*l))/omval;
        pder[phi0_] = 0;//-tanval * sin(omval * l) / (1 + omval * d0val);
        pder[z0_] = -l/cosDip();//- l - tanval * tanval * sin(omval * l) / (omval * (1 + omval * d0val));
        pder[tanl_] = 1/(cosDip()*cosDip());
        pder[t0_] = pder[z0_] / vz() + pder[tanl_] * (time - t0()) * cosval * cosval / tanval;
        break;
      case LocalBasis::phidir:
        // Azimuthal bending: R, Lambda, t0 are unchanged
        pder[d0_] = -l*sinTheta();//-sin(omval * l) / (omval * cosval);
        pder[phi0_] = 1/sinTheta();//cos(omval * l) / (cosval * (1 + omval * d0val));
        pder[z0_] = -pder[d0_]/(sinTheta()*tanTheta());//-tanval / (omval * cosval) * (1 - cos(omval * l) / (1 + omval * d0val));
        pder[tanl_] = 0;
        pder[t0_] = pder[z0_] / vz();
        break;
      case LocalBasis::momdir:
        // fractional momentum change: position and direction are unchanged
        pder[d0_] = 0;
        pder[phi0_] = 0;
        pder[z0_] = 0;//-tanval * (l - sin(omval * l) / (omval * (1 + omval * d0val)));
        pder[tanl_] = 0;
        pder[t0_] = pder[z0_] / vz();
        break;
      default:
        throw std::invalid_argument("Invalid direction");
    }
    return pder;    
*/

// compute some useful quantities
    //double vz = CLHEP::c_light * mom().z() / mom().E();
    double l = translen(CLHEP::c_light * beta() * (time - t0()));
    KTLine::DVEC pder;
    //cout << "Mom deriv start params " << pder << endl;
    // cases
    switch (mdir) {
    case LocalBasis::perpdir:
      // polar bending: change in Theta
      pder[tanl_] = 1/(cosDip()*cosDip());
      pder[d0_] = 0;
      pder[phi0_] = 0;
      pder[z0_] = -l * cosTheta(); // alt dir =-l*cosTheta();
      pder[t0_] = pder[z0_] / vz_ + pder[tanl_] * (time - t0()) * cosval * cosval / tanval;//pder[z0_] / vz;
      //cout << "Mom deriv perpdir params " << pder << endl;
      break;
    case LocalBasis::phidir:
      // change in dP/dtheta1 = dP/dphi0*(-1/sintheta)
      pder[tanl_] = 0;//GOOD
      pder[d0_] = -l;                 // alt dir = -l;
      pder[phi0_] = 1 / cosDip(); // alt dir = -1/sinTheta(); GOOD
      pder[z0_] = -d0() / (sinTheta() * tanTheta()); // alt dir = -d0()/(sinTheta()*tanTheta());
      pder[t0_] = pder[z0_] / vz_;
      //cout << "Mom deriv phidir params " << pder << endl;
      break;
    case LocalBasis::momdir:
      pder[tanl_] = 0;
      pder[d0_] = 0;
      pder[phi0_] = 0;
      pder[z0_] = 0;
      pder[t0_] = pder[z0_] / vz_;
      //cout << "Mom deriv momdir params " << pder << endl;
      break;

    default:
      throw std::invalid_argument("Invalid direction");
    }
    return pder;
  }

  void KTLine::print(ostream &ost, int detail) const {
    auto perr = params().diagonal();
    ost << " KTLine " << range() << " parameters: ";
    for (size_t ipar = 0; ipar < KTLine::npars_; ipar++) {
      ost << KTLine::paramName(static_cast<ParamIndex>(ipar)) << " "
          << paramVal(ipar) << " +- " << perr(ipar);
      if (ipar < KTLine::npars_ - 1)
        ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream &operator<<(ostream &ost, KTLine const &lhel) {
    lhel.print(ost, 0);
    return ost;
  }

} // namespace KinKal
