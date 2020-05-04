#ifndef KinKal_KTLine_hh
#define KinKal_KTLine_hh
/*
  Class to join the TLine with KInter for momentum. Create fit using KTrk.
  But before that we need to add in momentum interface here.
  As KTLine inherits from TLine, there is no need to define
  position and direction.

  S Middleton 2020
*/

#include "MatEnv/DetMaterial.hh"
#include "KinKal/KInter.hh"
#include "KinKal/TLine.hh"
#include <vector>
#include <stdexcept>
namespace KinKal {

  class KTLine :  public KInter , public TLine {
    public:

      constexpr static ParamIndex t0Index() { return t0_; }
      typedef ROOT::Math::SVector<double,npars_> PDER; // derivative of parameters type

      // This also requires the nominal BField, which can be a vector (3d) or a scalar (B along z)
      KTLine(Vec4 const& pos, Mom4 const& mom, double mass, int charge, Vec3 const& bnom, TRange const& range=TRange());
      KTLine(Vec4 const& pos, Mom4 const& mom, double mass, int charge, double bnom, TRange const& range=TRange());

      // construct from parameters
      KTLine(PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range=TRange())
      : _pars(pdata), _mass(mass), _charge(charge), _bnom( Vec3(0.0,0.0,bnom)), _trange(range);
      KTLine(PDATA const& pdata, double mass, int charge, double bnom, TRange const& range=TRange())
      : _pars(pdata), _mass(mass), _charge(charge), _bnom( Vec3(0.0,0.0,bnom)), _trange(range);

      //destructor:
      virtual ~KTLine() {}

      // particle momentum as a function of time
      void momentum(double t, Mom4& mom) const;

      // scalar momentum and energy in MeV/c units --> Needed for KKTrk:
      double momentum(double time) const  { return  gamma()*mass*beta(); }//in MeV/c
      double momentumVar(float time) const  { return -1.0; }//FIXME!
      double energy(double time) const  { return  sqrt(mass_*mass_ + gamma()*mass*beta()*gamma()*mass*beta()); }//in MeV E=sqrt(p^2 + m^2)

      // speed in mm/ns
      void print(std::ostream& ost, int detail) const override;
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const {};//infinity for striaght line

      // local momentum direction basis
      virtual void dirVector(MDir dir,double time,Vec3& unit) const;

      // momentum change derivatives; this is required to instantiate a KalTrk using this KInter
      void momDeriv(MDir mdir, double time, PDER& der) const;

      // some useful functions
      double ztime(double zpos) const { return t0() + zpos/(speed()*dir.z()); } //time to travel Z
      int charge() const { return charge_; }
      double beta() const { speed()/CLHEP::c_light;}// relativistic beta
      double gamma() const {1/sqrt(1-(speed()/CLHEP::c_light)*(speed()/CLHEP::c_light)));}// relativistic gamma
      Vec3 const& bnom() const { return bnom_; }//TODO - are these needed?


      //For the momentum, a magnitude and 4 mometnum set functions:
      double momMag() const{ //no time needed here.
        momentum_mag_ = gamma()*mass*beta();
        return mometnum_mag_;
      } //p[MeV/c] = gamma*M[MeV/c^2]*beta[c]
      void set_mom(Vec4 p){ mommentum_ = p; } //set p as a 4 vector
      Vec4 mom() const { return momentum_;}

    private :
      TRange trange_;//Do we need this?
      PDATA pars_; // parameters
      Vec3 bnom_; // nominal BField
      bool needsrot_; // logical flag if Bnom is parallel to global Z or not
      ROOT::Math::Rotation3D brot_; // rotation from the internal coordinate system (along B) to the global
      double momentum_mag_; //magnitude of momentum
      Vec4 momentum_; // 4 momentum vector - px,py,pz,m
      double mass_; //mass in MeV/c2
 };

}
#endif

  };

}
