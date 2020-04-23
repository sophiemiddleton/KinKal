#ifndef KinKal_KTLine_hh
#define KinKal_KTLine_hh
/*
  Class to join the TLine with KInter for momentum. Create fit using KTrk. But before that we need to add in momentum interface here. As KTLine inherits from TLine, there is no need to define position and direction.

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
      KTLine(Vec4 const& pos, Mom4 const& mom, int charge, Vec3 const& bnom, TRange const& range=TRange());
      KTLine(Vec4 const& pos, Mom4 const& mom, int charge, double bnom, TRange const& range=TRange());

      // construct from parameters
      KTLine(PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range=TRange());
      KTLine(PDATA const& pdata, double mass, int charge, double bnom, TRange const& range=TRange());

      //destructor:
      virtual ~KTLine() {}

      // particle momentum as a function of time
      void momentum(double t, Mom4& mom) const;

      // scalar momentum and energy in MeV/c units --> Needed for KKTrk:
      double momentum(double time) const  { return  mass_*pbar()/mbar_; }
      double momentumVar(float time) const  { return -1.0; }//FIXME!
      double energy(double time) const  { return  mass_*ebar()/mbar_; }

      // speed in mm/ns
      void print(std::ostream& ost, int detail) const override;
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const;

      // local momentum direction basis
      virtual void dirVector(MDir dir,double time,Vec3& unit) const;

      // momentum change derivatives; this is required to instantiate a KalTrk using this KInter
      void momDeriv(MDir mdir, double time, PDER& der) const;

      //TODO: Any useful functions?
      double ztime(double zpos) const { return t0() + zpos/(speed()*dir.z()); } //time to travel Z
      int charge() const { return charge_; }

      Vec3 const& bnom() const { return bnom_; }//TODO - are these needed?

      //For the momentum, a magnitude:
      double momMag() const{ return momentum_mag_;} ;
      void set_mom(double p){ momentum_mag_ = p; }
      void set_mom(Vec4 p){ _mommentum = p; }
      Vec4 mom() const { return momentum_;}

    private :
      PDATA pars_; // parameters
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      double pbar_; // here momentum is an input parameter, not calculated
      Vec3 bnom_; // nominal BField
      bool needsrot_; // logical flag if Bnom is parallel to global Z or not
      ROOT::Math::Rotation3D brot_; // rotation from the internal coordinate system (along B) to the global
      double momentum_mag_;
      Vec4 momentum_; // 4 momentum vector
 };

}
#endif

  };

}
