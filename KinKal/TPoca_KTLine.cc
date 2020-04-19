#include "KinKal/TPoca.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/KTLine.hh"
#include "KinKal/PKTraj.hh"
// specializations for TPoca
using namespace std;
namespace KinKal {
  

//*************** KTLine Stuff ***************** //
// S Middleton
// April 2020

//The following code is copied from above and adapted to the KTLine case.
//1) Specialization for KTLine:
template<> TPoca<KTLine,TLine>::TPoca(KTLine const& ktline, TLine const& tline, double precision) : TPocaBase(ktline,tline,precision)  { 
    // reset status
    reset();
    double ktltime,ltime; 
    // initialize the helix time using the Z position of the line
    // this will fail if the line is long and parallel to the z axis as there can be multiple solutions FIXME!
    htime = ktline.ztime(tline.z0());
    ltime = tline.t0();
    // use successive linear approximation until desired precision on DOCA is met.
    double ddoca(1.0e5);
    double doca(0.0);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter(0);
    // helix speed doesn't change
    double ktspeed = ktline.speed(lhelix.t0());
    Vec3 ktdir;
    while(fabs(ddoca) > precision_ && niter++ < maxiter) {
      // find helix local position and direction
      Vec3 ktpos;
      ktline.position(htime,hpos);
      ktline.direction(htime,hdir);
      auto dpos = tline.pos0()-hpos;
      // dot products
      double ddot = tline.dir().Dot(ktdir);
      double denom = 1.0 - ddot*ddot;
      // check for parallel)
      if(denom<1.0e-5){
	      status_ = TPoca::pocafailed;
	      break;
      }
      double ktdd = dpos.Dot(ktdir);
      double ldd = dpos.Dot(tline.dir());
      // compute length from expansion point to POCA and convert to times
      double ktlen = (ktdd - ldd*ddot)/denom;
      double llen = (ktdd*ddot - ldd)/denom;
      kttime += ktlen/ktspeed; // ktline time is iterative
      ltime = tline.t0() + llen/tline.speed(ltime);  // line time is always WRT t0
      // compute DOCA
      ktline.position(kttime,ktpos);
      Vec3 lpos;
      tline.position(ltime,lpos);
      double dd2 = (ktpos-lpos).Mag2();
      if(dd2 < 0.0 ){
	      status_ = TPoca::pocafailed;
	      break;
      }
      ddoca = doca;
      doca = sqrt(dd2);
      ddoca -= doca;
      if(isnan(ddoca)){
	      status_ = TPoca::pocafailed;
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
      poca_[0].SetE(kttime);
      ktline.position(poca_[0]);
      poca_[1].SetE(ltime);
      tline.position(poca_[1]);
      // sign doca by angular momentum projected onto difference vector
      double lsign = tline.dir().Cross(ktdir).Dot(poca_[1].Vect()-poca_[0].Vect());
      doca_ = copysign(doca,lsign);
    }
  }
//TODO: all these instances needed:
template<> TDPoca<KTLine,TLine>::TDPoca(TPoca<KTLine,TLine> const& tpoca) : TPoca<KTLine,TLine>(tpoca){}

template<> TDPoca<KTLine,TLine>::TDPoca(KTLine const& ktline, TLine const& tline, double precision) : TDPoca<KTLine,TLine>(TPoca<KTLine,TLine>(ktline,tline,precision)){}




}
