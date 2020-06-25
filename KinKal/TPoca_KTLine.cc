/*

To instantiate KKTrk on KTLine we need to specialize TPoca on the pair
<KTLine, TLine>. This class takes care of that.

S Middleton 2020

*/

#include "KinKal/TPoca.hh"
#include "KinKal/KTLine.hh"
#include "KinKal/TLine.hh"
#include "KinKal/KTLine.hh"
#include "KinKal/PKTraj.hh"
#include <limits>

// specializations for TPoca
using namespace std;
namespace KinKal {

//1) Specialization for KTLine:
  template<> TPoca<KTLine,TLine>::TPoca(KTLine const& ktline, TLine const& tline, TPocaHint const& hint, double precision) : TPocaBase(precision),ktraj_(&ktline), straj_(&tline) {
    // reset status
    reset();
     double htoca, stoca;
    // similar for line; this shouldn't matter, since the solution is linear
    if(hint.particleHint_){
      htoca = hint.particleToca_;
    }else
     htoca= ktline.ztime(tline.z0());
    if(hint.sensorHint_)
      stoca = hint.sensorToca_;
    else
      stoca=  tline.t0();

    // use successive linear approximation until desired precision on DOCA is met.
    double dptoca(std::numeric_limits<double>::max()), dstoca(std::numeric_limits<double>::max());

    // use successive linear approximation until desired precision on DOCA is met.
    double doca(0.0);
    static const unsigned maxiter=100;
    unsigned niter(0);
    // ktline speed doesn't change
    double ktspeed = ktline.speed(ktline.t0());
    Vec3 ktdir;
    while((fabs(dptoca) > precision_ || fabs(dptoca) > precision_  )&& niter++ < maxiter) {
      // find line's local position and direction
      Vec3 ktpos = ktline.position(htoca);
      ktdir = ktline.direction(htoca);//,ktdir);
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
      dptoca = (ktdd-ldd*ddot)/(denom*ktspeed);

      dstoca = tline.t0() + (ktdd*ddot-ldd)/(denom*tline.speed(stoca)) - stoca;
      htoca += dptoca; // ktline time is iterative
      stoca += dstoca; // line time is always WRT t0, since it uses p0
      // compute DOCA
      ktpos = ktline.position(htoca);
      Vec3 lpos =  tline.position(stoca);
      double dd2 = (ktpos-lpos).Mag2();;
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
      partPoca_.SetE(htoca);
      ktline.position(partPoca_);
      sensPoca_.SetE(stoca);
      tline.position(sensPoca_);
      // sign doca by angular momentum projected onto difference vector (same as helix)
      double lsign =   tline.dir().Cross(ktline.direction(partPoca_.T())).Dot(sensPoca_.Vect()-partPoca_.Vect());
      double dsign = copysign(1.0,lsign);
      doca_ = doca*dsign;

      // pre-compute some values needed for the derivative calculations
      Vec3 ddir;
      ddir = delta().Vect().Unit();// direction vector along D(POCA) from traj 2 to 1 (line to ktline)
      ktline.direction(particlePoca().T());//hidr

      // derviatives of TOCA and DOCA WRT particle trajectory parameters
      // no t0 dependence, DOCA is purely geometric

      //calculated these using BTrk instances - doc db ref ###
      dDdP_[KTLine::phi0_] = -dsign*(ktline.d0()*sin(ktline.phi0())*ddir.x()+ktline.d0()*cos(ktline.phi0())*ddir.y());
      dDdP_[KTLine::cost_] = 0;
      dDdP_[KTLine::d0_] = -dsign*(-1*cos(ktline.phi0())*ddir.x()+sin(ktline.phi0())*ddir.y());
      dDdP_[KTLine::z0_] = -dsign*ddir.z();
      // no spatial dependence, DT is purely temporal
      dTdP_[KTLine::t0_] = -1.0; // time is 100% correlated

      // propagate parameter covariance to variance on doca and toca
      docavar_ = ROOT::Math::Similarity(dDdP(),ktline.params().covariance());
      tocavar_ = ROOT::Math::Similarity(dTdP(),ktline.params().covariance());
      // dot product between directions at POCA
  
      ddot_ = ktline.direction(particleToca()).Dot(tline.direction(sensorToca()));

    }
  }

   // specialization between a piecewise KTLine and a line
  typedef PKTraj<KTLine> PKTLINE;
  template<> TPoca<PKTLINE,TLine>::TPoca(PKTLINE const& pktline, TLine const& tline, TPocaHint const& hint, double precision) : TPocaBase(precision), ktraj_(&pktline), straj_(&tline)  {
    // iteratively find the nearest piece, and POCA for that piece.  Start at hints if availalble, otherwise the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter=0;
    size_t oldindex= pktline.pieces().size(); //TODO --->do we want piecewise line fit?
    size_t index;
    if(hint.particleHint_){
      index = pktline.nearestIndex(hint.particleToca_);
    }else{
      index = size_t(rint(oldindex/2.0));
    }
    status_ = converged; 
    while(status_ == converged && niter++ < maxiter && index != oldindex){
      // call down to KTLine TPoca
      // prepare for the next iteration
      KTLine const& piece = pktline.pieces()[index];
      TPoca<KTLine,TLine> tpoca(piece,tline,hint,precision);
      status_ = tpoca.status();
      if(tpoca.usable()){
	      // copy over the rest of the state
	      partPoca_ = tpoca.particlePoca();
	      sensPoca_ = tpoca.sensorPoca();
	      doca_ = tpoca.doca();
	      dDdP_ = tpoca.dDdP();
	      dTdP_ = tpoca.dTdP();
	      docavar_ = tpoca.docaVar();
	      tocavar_ = tpoca.tocaVar();
	      ddot_ = tpoca.dirDot();
      }
      oldindex = index;
      index = pktline.nearestIndex(tpoca.particlePoca().T());
    }
    if(status_ == converged && niter >= maxiter) status_ = unconverged;
  }

}
