#ifndef KinKal_SimpleWireHit_hh
#define KinKal_SimpleWireHit_hh
//
// Simple implementation of a wire hit based on a Straw detector.
//
#include "KinKal/Detector/WireHit.hh"
namespace KinKal {

// struct for updating simple wire hits; a real wire hit updater needs to know about calibrations, etc
  struct SimpleWireHitUpdater {
    double mindoca_; // minimum DOCA value to set an ambiguity
    double maxdoca_; // maximum DOCA to still use a hit
    bool nulltime_; // constrain time when hit has null ambiguity
    SimpleWireHitUpdater(double mindoca,double maxdoca, bool nulltime) : mindoca_(mindoca), maxdoca_(maxdoca), nulltime_(nulltime) {}
  };

  template <class KTRAJ> class SimpleWireHit : public WireHit<KTRAJ> {
    public:
      using WIREHIT = Hit<KTRAJ>;
      using update = WIREHIT::update;
      SimpleWireHit(BFieldMap const& bfield, Line const& wire, EXINGPTR const& dxing, WireHitState const& whstate,
      double driftspeed, double tvar, double rcell);
// WireHit and Hit interface implementations
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const override;
      // specific to SimpleWireHit: this has a constant drift speed
      virtual ~SimpleWireHit(){}
      double timeVariance() const { return tvar_; }
      double radius() const { return rcell_; }
    private:
      double dvel_; // constant drift speed
      double tvar_; // constant time variance
      double rcell_; // straw radius
  };

  template <class KTRAJ> SimpleWireHit<KTRAJ>::SimpleWireHit(BFieldMap const& bfield, Line const& wire, EXINGPTR const& dxing, WireHitState const& whstate,
      double driftspeed, double tvar, double rcell) :
    WireHit(bfield,wire,dxing,whstate), dvel_(driftspeed), tvar_(tvar), rcell_(rcell) {}


  template <class KTRAJ> void SimpleWireHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // set precision
    precision_ = miconfig.tprec_;
    // update to move to the new trajectory
    update(pktraj);
    // find the wire hit updater in the update params.  There should be 0 or 1
    const SimpleWireHitUpdater* whupdater(0);
    for(auto const& uparams : miconfig.hitupdaters_){
      auto const* whu = std::any_cast<SimpleWireHitUpdater>(&uparams);
      if(whu != 0){
	if(whupdater !=0) throw std::invalid_argument("Multiple SimpleWireHitUpdaters found");
	whupdater = whu;
      }
    }
    // crude updating of ambiguity and activity based on DOCA
    if(whupdater != 0){
      // start with existing state
      WireHitState newstate = WIREHIT::hitstate();
      newstate.nullvar_ = whupdater->mindoca_*whupdater->mindoca_/3.0; // RMS of flat distribution beteween +- mindoca
      double doca = fabs(tresid_.tPoca().doca());
      if(fabs(doca) > whupdater->mindoca_){
	newstate.lrambig_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
	newstate.dimension_ = WireHitState::time;
      } else if( fabs(doca) > whupdater->maxdoca_){
	newstate.dimension_ = WireHitState::none; // disable the hit
      } else {
	newstate.lrambig_ = LRAmbig::null;
	if(whupdater->nulltime_)
	  newstate.dimension_ = WireHitState::both;
	else
	  newstate.dimension_ = WireHitState::distance;
      }
      setHitState(newstate);
      // now update again in case the hit changed
      update(pktraj);
    }
    // OK if no updater is found, hits may be frozen this meta-iteration
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    // simply translate distance to time using the fixed velocity
    dinfo.tdrift_ = drift.R()/dvel_;
    dinfo._dspeed_ = dvel_;
    dinfo.tdriftvar_ = tvar_;
  }

}
#endif
