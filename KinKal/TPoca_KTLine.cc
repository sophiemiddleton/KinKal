/*

To instantiate KKTrk on KTLine we need to specialize TPoca on the pair
<KTLine, TLine>. This class takes care of that.

S Middleton 2020

*/

#include "KinKal/TPoca.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/KTLine.hh"
#include "KinKal/PKTraj.hh"
// specializations for TPoca
using namespace std;
namespace KinKal {


//*************** KTLine Stuff ***************** //
//The following code is copied from Helix instance and adapted to the KTLine case.
//1) Specialization for KTLine:
template<> TPoca<KTLine,TLine>::TPoca(KTLine const& ktline, TLine const& tline, double precision) : TPocaBase(ktline,tline,precision)  {
    // reset status
    reset();
     float ktltime,ltime;
    // similar for line; this shouldn't matter, since the solution is linear
    if(hint.particleHint_)
      ktltime = hint.particleToca_;
    else
      ktltime = ktline.ztime(tline.z0());
    if(hint.sensorHint_)
      ltime = hint.sensorToca_;
    else
      ltime= tline.t0();

    // use successive linear approximation until desired precision on DOCA is met.
    float dptoca(std::numeric_limits<float>::max()), dstoca(std::numeric_limits<float>::max());

    // use successive linear approximation until desired precision on DOCA is met.
    double doca(0.0);
    static const unsigned maxiter=100;
    unsigned niter(0);
    // ktline speed doesn't change
    double ktspeed = ktline.speed(ktline.t0());
    Vec3 ktdir;
    while(fabs(dpoca) > precision_ || fabs(dpoca) > precision_  && niter++ < maxiter) {
      // find line's local position and direction
      Vec3 ktpos;
      ktline.position(ktltime, ktpos);
      ktline.direction(ktltime,ktdir);
      auto dpos = tline.pos0()-ktpos;
      // dot products
      double ddot = tline.dir().Dot(ktdir);
      double denom = 1.0 - ddot*ddot;

      // check for parallel
      if(denom<1.0e-5){
	      status_ = TPoca::pocafailed;
	      break;
      }
      double ktdd = dpos.Dot(ktdir);
      double ldd = dpos.Dot(tline.dir());
      // compute length from expansion point to POCA and convert to times
      double ktlen = (ktdd - ldd*ddot)/denom;
      double llen = (ktdd*ddot - ldd)/denom;
      dptoca = (ktlen*ktspeed);
      dstoca = tline.t0() + llen/tline.speed(ltime)) - ltime;
      ktltime += dptoca; // ktline time is iterative
      ltime += dstoca; // line time is always WRT t0, since it uses p0

      // compute DOCA
      ktline.position(kttime,ktpos);
      Vec3 lpos;
      tline.position(ltime,lpos);
      double dd2 = (ktpos-lpos).Mag2();
      if(dd2 < 0.0 ){
	      status_ = TPoca::pocafailed;
	      break;
      }
      doca = sqrt(dd2);
       // update convergence test
      if(isnan(doca)){
	      status_ = pocafailed;
	      break;
      }
    }
    // finalize TPoca
    if(status_ != TPoca::pocafailed){
      if(niter < maxiter)
	      status_ = TPoca::converged;
      else
	      status_ = TPoca::unconverged;
        // set the positions
      partPoca_[0].SetE(kttime);
      ktline.position(partPoca_[0]);
      partPoca_[1].SetE(ltime);
      tline.position(sensPoca_[1]);
      // sign doca by angular momentum projected onto difference vector (same as helix)
      double lsign = tline.dir().Cross(ktdir).Dot(partPoca_[1].Vect()-partPoca_[0].Vect());
      float dsign = copysign(1.0,lsign);
      doca_ = doca*dsign;

      // pre-compute some values needed for the derivative calculations
      Vec3 vdoca, ddir, hdir;
      delta(vdoca);
      ddir = vdoca.Unit();// direction vector along D(POCA) from traj 2 to 1 (line to ktline)
      ktline.direction(particlePoca().T(),hdir);
//TODO - look at the BTrk version (TrkMomCalc)
      // derviatives of TOCA and DOCA WRT particle trajectory parameters
      // no t0 dependence, DOCA is purely geometric
      float ktlphi = ktline.dir().Phi(); // local azimuth of ktline direction
      float lphi = tline.dir().Phi(); // line azimuth
      float d = sqrt((-1*ktline.d0()*-1*ktline.d0()) + (ktline.z0()*ktline.z0()))
      //calculated these using BTrk instances - doc db ref ###
      dDdP_[KTLine::d0_] = 1/(2*d);
      dDdP_[KTLine::cost_] = 0;
      dDdP_[KTLine::phi0_] = 0; //cos^2+sin^2 = 1 so phi0 factors out.
      dDdP_[KTLine::z0_] = 1/(2*d);

      // no spatial dependence, DT is purely temporal
      dTdP_[KTLine::t0_] = 1.0; // time is 100% correlated

      // propagate parameter covariance to variance on doca and toca
      docavar_ = ROOT::Math::Similarity(dDdP(),ktline.params().covariance());
      tocavar_ = ROOT::Math::Similarity(dTdP(),ktline.params().covariance());
      // dot product between directions at POCA
      Vec3 pdir, sdir;
      ktline.direction(particleToca(),pdir);
      tline.direction(sensorToca(),sdir);
      ddot_ = pdir.Dot(sdir);


    }
  }


}
