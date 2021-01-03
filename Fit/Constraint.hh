#ifndef KinKal_Constraint_hh
#define KinKal_Constraint_hh
//
//  class represeting a constraint on the fit parameters due to a hit.  A hit is anything that adds information content to the fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Detector/Hit.hh"
#include <ostream>
#include <memory>

namespace KinKal {
  template <class KTRAJ> class Constraint : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;
      
      Chisq chisq(Parameters const& pdata) const override;
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) override;
      void process(FitState& kkdata,TimeDir tdir) override;
      bool active() const override { return hit_->active(); }
      double time() const override { return hit_->time(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual ~Constraint(){}
      // local functions
      // construct from a hit and reference trajectory
      Constraint(HITPTR const& hit, PKTRAJ const& reftraj,double precision=1e-6);
      // the unbiased parameters are the fit parameters not including the information content of this effect
      Parameters unbiasedParameters() const;
      // Total unbiased chisquared for this effect
      Chisq chisq() const;
      // access the contents
      HITPTR const& hit() const { return hit_; }
      Weights const& weightCache() const { return wcache_; }
      Weights const& hitWeight() const { return hitwt_; }
      double precision() const { return precision_; }
    private:
      HITPTR hit_ ; // hit used for this constraint
      Weights wcache_; // sum of processing weights in opposite directions, excluding this hit's information. used to compute unbiased parameters and chisquared
      Weights hitwt_; // weight representation of the hits constraint
      double vscale_; // variance factor due to annealing 'temperature'
      double precision_; // precision used in TCA calcuation
  };

  template<class KTRAJ> Constraint<KTRAJ>::Constraint(HITPTR const& hit, PKTRAJ const& reftraj,double precision) : hit_(hit), vscale_(1.0), precision_(precision) {
    update(reftraj);
  }
 
  template<class KTRAJ> void Constraint<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    // direction is irrelevant for processing hits 
    if(this->active()){
      // cache the processing weights, adding both processing directions
      wcache_ += kkdata.wData();
      // add this effect's information
      kkdata.append(hitwt_);
    }
    KKEFF::setState(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void Constraint<KTRAJ>::update(PKTRAJ const& pktraj) {
    // reset the processing cache
    wcache_ = Weights();
    // update the hit
    hit_->update(pktraj);
    // get the weight from the hit 
    hitwt_ = hit_->weight();
    // scale weight for the temp
    hitwt_ *= 1.0/vscale_;
    // ready for processing!
    KKEFF::updateState();
  }

  template<class KTRAJ> void Constraint<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // reset the annealing temp and hit precision
    vscale_ = miconfig.varianceScale();
    precision_ = miconfig.tprec_;
    // update the hit's internal state; this depends on the configuration parameters
    if(miconfig.updatehits_)hit_->update(pktraj,miconfig );
    // update the state of this object
    update(pktraj);
  }

  template<class KTRAJ> Chisq Constraint<KTRAJ>::chisq(Parameters const& pdata) const {
    if(this->active()) {
      double chi2(0.0);
      for(size_t idof= 0; idof < hit_->nDOF(); idof++){
	double chi = hit_->chi(idof,pdata);
	chi2 += chi*chi;
      }	
      // correct for current variance scaling
      double chi2  /= vscale_;
      return Chisq(chi2,hit_->nDOF());
    } else
      return Chisq();
  }

  template<class KTRAJ> Chisq Constraint<KTRAJ>::chisq() const {
   if(this->active()) {
     Parameters unbiased = unbiasedParameters();
      return chisq(unbiased);
    } else
      return Chisq();
  } 

  template<class KTRAJ> Parameters Constraint<KTRAJ>::unbiasedParameters() const {
  // this function can't be called on an unprocessed effect
    if( !KKEFF::wasProcessed(TimeDir::forwards) || !KKEFF::wasProcessed(TimeDir::backwards))
      throw  std::invalid_argument("Can't compute unbiased parameters for unprocessed constraint");
    // Invert the cache to get unbiased parameters at this constraint
      return Parameters(wcache_);
  }

  template <class KTRAJ> void Constraint<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << "Constraint " << static_cast<Effect<KTRAJ> const&>(*this) << std::endl;
    if(detail > 0){
      hit_->print(ost,detail);    
      ost << " Constraint Weight " << hitwt_ << std::endl;
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Constraint<KTRAJ> const& kkhit) {
    kkhit.print(ost,0);
    return ost;
  }

}
#endif
