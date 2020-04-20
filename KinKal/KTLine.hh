#ifndef KinKal_KTLine_hh
#define KinKal_KTLine_hh
//
//  Class to join the TLine with KIter for momentum. Create fit using KTrk. This is based on LHelix.
//
#include "MatEnv/DetMaterial.hh"
#include "KinKal/KInter.hh"
#include "KinKal/TLine.hh"
//#include "KinKal/KTrk.hh" //need to initiate this on the KIter subclass
#include <vector>
#include <stdexcept>
namespace KinKal {

  class KTLine :  public KInter , public TLine {
    public:

      // As KTLine inherits from TLine, there is no need to define too much here.
      constexpr static ParamIndex t0Index() { return t0_; }
      typedef ROOT::Math::SVector<double,npars_> PDER; // derivative of parameters type 
      
      // construct from momentum, position, and particle properties.
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
      double energy(double time) const  { return  mass_*ebar()/mbar_; }

      // speed in mm/ns
      void print(std::ostream& ost, int detail) const override;
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const;

      // local momentum direction basis
      virtual void dirVector(MDir dir,double time,Vec3& unit) const;

      // momentum change derivatives; this is required to instantiate a KalTrk using this KInter
      void momDeriv(MDir mdir, double time, PDER& der) const;

      //TODO: check these things are useful
      // simple functions; these can be cached if they cause performance problems
      double pbar2() const { return  pbar_*pbar_; } 
      double pbar() const { return  pbar_; } // momentum in mm
      double ebar2() const { return  pbar_*pbar_ + mbar_*mbar_; }
      double ebar() const { return  sqrt(ebar2()); } // energy in mm
      double mbar() const { return mbar_; } // mass in mm; includes charge information!
      double Q() const { return mbar_/mass_; } // reduced charge
     
      double beta() const { return pbar()/ebar(); } // relativistic beta
      double gamma() const { return fabs(ebar()/mbar_); } // relativistic gamma
      double ztime(double zpos) const { return t0() + zpos/(speed()*dir.z()); } //time to travel Z
      int charge() const { return charge_; }


      Vec3 const& bnom() const { return bnom_; }//Needed?

      //For the momentum, a magnitude:
      double momMag() const{ return _mommentum_mag;} ;
      void set_mom(double p){ _mommentum_mag = p; }
      void set_mom(Vec4 p){ _mommentum = p; }
      Vec4 mom() const { return momentum_;}

    private :
      PDATA pars_; // parameters
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      double pbar_; // here momentum is an input parameter, not calculated
      Vec3 bnom_; // nominal BField
      bool needsrot_; // logical flag if Bnom is parallel to global Z or not
      ROOT::Math::Rotation3D brot_; // rotation from the internal coordinate system (along B) to the global
      double _mommentum_mag;
      Vec4 momentum_; // 4 momentum vector
 };
 
}
#endif

  };

}
