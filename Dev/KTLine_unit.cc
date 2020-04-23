//
// test basic functions of the KTLine TTraj class
//
#include "KinKal/ktlineix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TAxis3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;

void print_usage() {
  printf("Usage: KTLine  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --tmin f--tmax f --ltime f --By f \n");
}

struct MomVec {
  TPolyLine3D* arrow;
  TPolyMarker3D *start, *end;
  MomVec() : arrow(new TPolyLine3D(2)), start(new TPolyMarker3D(1,21)), end(new TPolyMarker3D(1,22)) {}
};

void drawMom(Vec3 const& start, Vec3 const& momvec,int momcolor,MomVec& mom) {
  mom.arrow->SetPoint(0,start.X(),start.Y(),start.Z());
  auto end = start + momvec;
  mom.arrow->SetPoint(1,end.X(),end.Y(),end.Z());
  mom.arrow->SetLineColor(momcolor);
  mom.arrow->Draw();
  mom.start->SetPoint(1,start.X(),start.Y(),start.Z());
  mom.start->SetMarkerColor(momcolor);
  mom.start->Draw();
  mom.end->SetPoint(1,end.X(),end.Y(),end.Z());
  mom.end->SetMarkerColor(momcolor);
  mom.end->Draw();
}

int main(int argc, char **argv) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass, oz(100.0), ot(0.0);
  double tmin(-5.0), tmax(5.0);
  double ltime(3.0), vprop(0.8), gap(2.0);
  double hlen(500.0); // half-length of the wire
  double By(0.0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 't'  },
    {"tmin",     required_argument, 0, 's'  },
    {"tmax",     required_argument, 0, 'e'  },
    {"ltime",     required_argument, 0, 'l'  },
    {"By",     required_argument, 0, 'y'  },
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'c' : cost = atof(optarg);
		 break;
      case 'a' : phi = atof(optarg);
		 break;
      case 'p' : imass = atoi(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'z' : oz = atof(optarg);
		 break;
      case 't' : ot = atof(optarg);
		 break;
      case 's' : tmin = atof(optarg);
		 break;
      case 'e' : tmax = atof(optarg);
		 break;
      case 'l' : ltime = atof(optarg);
		 break;
      case 'y' : By = atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }

  pmass = masses[imass];

  printf("Testing KTLine with momentum = %f, costheta = %f, phi = %f, mass = %f, charge = %i, z = %f, t = %f \n",mom,cost,phi,pmass,icharge,oz,ot);
// define the BF (tesla)
  Vec3 bnom(0.0,By,1.0);
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(0,1000,0); //TODO - check---->
  //This is how we will send in the momentum (this test is a comsic of energy 1GeV/c coming straight down)
  KTLine ktline(origin,momv,icharge,bnom);
  Mom4 testmom;
  ktline.momentum(ot,testmom);
  cout << "KTLine with momentum " << testmom << " position " << origin << " has parameters: " << ktline << endl;
  Vec3 vel;
  ktline.velocity(ot,vel); //TODO - velocity needs a speed, we need to set speed!
  double dot = vel.Dot(testmom)/CLHEP::c_light;
  cout << "velocity dot mom = " << dot << endl;
  cout << "momentum beta =" << momv.Beta() << " ktline beta = " << ktline.beta() << endl;
  Vec3 mdir;
  ktline.direction(ot,mdir);
  // create the helix at tmin and tmax
  Mom4 tmom;
  Vec4 tpos;
  ktline.momentum(tmax,tmom);
  tpos.SetE(tmax);
  ktline.position(tpos);
  ktlineix ktlinemax(tpos,tmom,icharge,bnom); //TODO
  ktline.momentum(tmin,tmom);
  tpos.SetE(tmin);
  ktline.position(tpos);
  ktlineix ktlinemin(tpos,tmom,icharge,bnom); //TODO

  cout << "KTLine at tmax has parameters : " << ktlinemax << endl;
  cout << "KTLine at tmin has parameters : " << ktlinemin << endl;

// now graph this as a polyline over the specified time range.
  double tstep = 0.1; // nanoseconds
  double trange = tmax-tmin;
  int nsteps = (int)rint(trange/tstep);
// create Canvase
  TCanvas* hcan = new TCanvas("hcan","KTLine",1000,1000);
//TPolyLine to graph the result
  TPolyLine3D* lin = new TPolyLine3D(nsteps+1);
  Vec4 hpos;
  for(int istep=0;istep<nsteps+1;++istep){
  // compute the position from the time
    hpos.SetE(tmin + tstep*istep);
    ktline.position(hpos); //TODO - check this works!!
    // add these positions to the TPolyLine3D
    lin->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  if(icharge > 0)
    lin->SetLineColor(kBlue);
  else
    lin->SetLineColor(kRed);
  lin->Draw();

  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
  rulers->Draw();

// now draw momentum vectors at reference, start and end
  MomVec mstart,mend,mref;
  Vec3 mompos;
  ktline.position(ot,mompos);
  ktline.direction(ot,mdir);
  Vec3 momvec =mom*mdir;
  drawMom(mompos,momvec,kBlack,mref);
  //
  ktline.position(tmin,mompos);
  ktline.direction(tmin,mdir);
  momvec =mom*mdir;
  drawMom(mompos,momvec,kBlue,mstart);
  //
  ktline.position(tmax,mompos);
  ktline.direction(tmax,mdir);
  momvec =mom*mdir;
  drawMom(mompos,momvec,kGreen,mend);
  //
  TLegend* leg = new TLegend(0.8,0.8,1.0,1.0);
  char title[80];
  snprintf(title,80,"KTLine, q=%1i, mom=%3.1g MeV/c",icharge,mom);
  leg->AddEntry(lin,title,"L");
  snprintf(title,80,"Ref. Momentum, t=%4.2g ns",ot);
  leg->AddEntry(mref.arrow,title,"L");
  snprintf(title,80,"Start Momentum, t=%4.2g ns",ot+tmin);
  leg->AddEntry(mstart.arrow,title,"L");
  snprintf(title,80,"End Momentum, t=%4.2g ns",ot+tmax);
  leg->AddEntry(mend.arrow,title,"L");
  leg->Draw();

  // create a TLine near this KTLine, and draw it and the TPoca vector
  Vec3 pos, dir;
  ktline.position(ltime,pos);
  ktline.direction(ltime,dir);
  // rotate the direction
  double ktline_phi = atan2(dir.Y(),dir.X());
  double pphi = ktline_phi + M_PI/2.0;
  Vec3 pdir(cos(pphi),sin(pphi),0.0);
  double pspeed = CLHEP::c_light*vprop; // vprop is relative to c
  Vec3 pvel = pdir*pspeed;
  // shift the position
  Vec3 perpdir(-sin(phi),cos(phi),0.0);
  Vec3 ppos = pos + gap*perpdir;
// time range;
  TRange prange(ltime-hlen/pspeed, ltime+hlen/pspeed);
  TLine tline(ppos, pvel,ltime,prange);
// find TPoca
  TPoca<KTLine,TLine> tp(ktline,tline); //TODO - check POCA can do this
  cout << "TPoca status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
  if(tp.status() == TPocaBase::converged) {
    // draw the line and TPoca
    TPolyLine3D* line = new TPolyLine3D(2);
    Vec3 plow, phigh;
    tline.position(tline.range().low(),plow);
    tline.position(tline.range().high(),phigh);
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->SetLineColor(kOrange);
    line->Draw();
    TPolyLine3D* poca = new TPolyLine3D(2);
    poca->SetPoint(0,tp.poca0().X() ,tp.poca0().Y() ,tp.poca0().Z());
    poca->SetPoint(1,tp.poca1().X() ,tp.poca1().Y() ,tp.poca1().Z());
    poca->SetLineColor(kBlack);
    poca->Draw();
  }

  snprintf(title,80,"KTLine_m%3.1f_p%3.2f_q%i.root",pmass,mom,icharge);
  cout << "Saving canvas to " << title << endl;
  hcan->SaveAs(title);

  exit(EXIT_SUCCESS);
}
