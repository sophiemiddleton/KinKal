/*
 Test derivatives of the KTLine TTraj class.

 S Middleton 2020
*/
#include "KinKal/KTLine.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TArrow.h"
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

void print_usage() {
  printf("Usage: KTLineDerivs  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --dmin f --dmax f --ttest f --By f \n");
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  // save canvases
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass, oz(100.0), ot(0.0), ttest(5.0);
  double dmin(-5e-2), dmax(5e-2);
  double By(0.0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 'o'  },
    {"dmin",     required_argument, 0, 's'  },
    {"dmax",     required_argument, 0, 'e'  },
    {"ttest",     required_argument, 0, 't'  },
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
      case 'o' : ot = atof(optarg);
		 break;
      case 's' : dmin = atof(optarg);
		 break;
      case 'e' : dmax = atof(optarg);
		 break;
      case 't' : ttest = atof(optarg);
		 break;
      case 'y' : By = atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
  Vec3 bnom(0.0,By,1.0);
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);

  pmass = masses[imass];
  Mom4 momv(mom*sint*sin(phi),mom*sint*cos(phi),mom*cost,pmass); //Derived from direction*mom
  KTLine ref_ktline(origin,momv,pmass,icharge,bnom);
  cout << "Reference " << ref_ktline << endl;
  Vec4 refpos4;
  refpos4.SetE(ttest);
  ref_ktline.position(refpos4);
  cout << "origin position " << origin << " test position " << refpos4 << endl;
  Mom4 refmom;
  ref_ktline.momentum(ttest,refmom);
  int ndel(50);
  // graphs to compare parameter change
  TGraph* d0graph[3];
  TGraph* costgraph[3];
  TGraph* t0graph[3];
  TGraph* z0graph[3];
  TGraph* phi0graph[3];

  // graphs to compare momentum change
  TGraph* mom0graph[3];
  TGraph* mom1graph[3];
  TGraph* mom2graph[3];
  // position change
  TGraph* posgraph[6];
  // gaps
  TGraph* gapgraph[3];
  // canvases
  TCanvas* dhcan[3];
  TCanvas* dmomcan[3];
  TFile ktline_deriv("KTLineDerivs.root","RECREATE");
  // loop over derivative directions
  double del = (dmax-dmin)/(ndel-1);
  for(int idir=0;idir<3;++idir){
    KInter::MDir tdir =static_cast<KInter::MDir>(idir);
    Vec3 dmomdir;
    ref_ktline.dirVector(tdir,ttest,dmomdir);

    // parameter change
    d0graph[idir] = new TGraph(ndel);
    d0graph[idir]->SetTitle("d_{0};exact;1st derivative");
    costgraph[idir] = new TGraph(ndel);
    costgraph[idir]->SetTitle("cost;exact;1st derivative");
    t0graph[idir] = new TGraph(ndel);
    t0graph[idir]->SetTitle("t_{0};exact;1st derivative");
    phi0graph[idir] = new TGraph(ndel);
    phi0graph[idir]->SetTitle("phi_{0};exact;1st derivative");
    z0graph[idir] = new TGraph(ndel);
    z0graph[idir]->SetTitle("z_{0};exact;1st derivative");

    mom0graph[idir] = new TGraph(ndel);
    mom0graph[idir]->SetTitle("Momentum Direction;exact;1st derivative");
    mom1graph[idir] = new TGraph(ndel);
    mom1graph[idir]->SetTitle("Theta Direction;exact;1st derivative");
    mom2graph[idir] = new TGraph(ndel);
    mom2graph[idir]->SetTitle("Phi Direction;exact;1st derivative");
    gapgraph[idir] = new TGraph(ndel);
    gapgraph[idir]->SetTitle("Gap;change;gap value (mm)");

    // scan range of change
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id;
      //  compute exact altered params
      Vec3 newmom = refmom.Vect() + delta*dmomdir*mom;
      Mom4 momv(newmom.X(),newmom.Y(),newmom.Z(),pmass);
      KTLine xktline(refpos4,momv,pmass,icharge,bnom);
      // now, compute 1st order change in parameters
      KTLine::PDER pder;
      ref_ktline.momDeriv(tdir,ttest,pder);
      cout << "derivative vector" << pder << endl;
      auto dvec = ref_ktline.params().parameters() + delta*pder;
      KTLine::PDATA pdata(dvec,ref_ktline.params().covariance());
      KTLine dktline(pdata,ref_ktline.mass(),ref_ktline.charge(),bnom);

      Vec4 xpos, dpos;
      xpos.SetE(ttest);
      dpos.SetE(ttest);
      xktline.position(xpos);
      dktline.position(dpos);
      Mom4 dmom;
      dktline.momentum(ttest,dmom);
      Vec4 gap = dpos - refpos4;
      gapgraph[idir]->SetPoint(id,delta,sqrt(gap.Vect().Mag2()));
      // parameter diff
      d0graph[idir]->SetPoint(id,xktline.d0()-ref_ktline.d0(),dktline.d0()-ref_ktline.d0());
      costgraph[idir]->SetPoint(id,xktline.cost()-ref_ktline.cost(),dktline.cost()-ref_ktline.cost());
      t0graph[idir]->SetPoint(id,xktline.t0()-ref_ktline.t0(),dktline.t0()-ref_ktline.t0());
      phi0graph[idir]->SetPoint(id,xktline.phi0()-ref_ktline.phi0(),dktline.phi0()-ref_ktline.phi0());
      z0graph[idir]->SetPoint(id,xktline.z0()-ref_ktline.z0(),dktline.z0()-ref_ktline.z0());
      // compare momenta after change
      Vec3 dxmom = momv.Vect() - refmom.Vect();
      Vec3 ddmom = dmom.Vect() - refmom.Vect();
      Vec3 changedir;
      ref_ktline.dirVector(KInter::momdir,ttest,changedir);
      mom0graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
      ref_ktline.dirVector(KInter::theta1,ttest,changedir);
      mom1graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
      ref_ktline.dirVector(KInter::theta2,ttest,changedir);
      mom2graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
    }
    char title[80];
    char name[80];
    snprintf(name,80,"dh%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"KTLine Change %s",KInter::directionName(tdir).c_str());
    dhcan[idir] = new TCanvas(name,title,1200,800);
    dhcan[idir]->Divide(3,2);
    dhcan[idir]->cd(1);
    d0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(2);
    costgraph[idir]->Draw("AC*");
    dhcan[idir]->cd(3);
    t0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(4);
    phi0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(5);
    z0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(6);

    dhcan[idir]->Write();

    snprintf(name,80,"dm%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"Mom Change %s",KInter::directionName(tdir).c_str());
    dmomcan[idir] = new TCanvas(name,title,800,800);
    dmomcan[idir]->Divide(2,2);
    dmomcan[idir]->cd(1);
    mom0graph[idir]->Draw("AC*");
    dmomcan[idir]->cd(2);
    mom1graph[idir]->Draw("AC*");
    dmomcan[idir]->cd(3);
    mom2graph[idir]->Draw("AC*");
    dmomcan[idir]->cd(4);
    gapgraph[idir]->Draw("AC*");
    dmomcan[idir]->Draw();
    dmomcan[idir]->Write();
  }

  // now spatial derivatives: these are used to constrain continuity between traj pieces
  KTLine::PDER pder;
  ref_ktline.posDeriv(ttest,pder);
  Vec3 refpos, refdir;
  ref_ktline.position(ttest,refpos);
  ref_ktline.direction(ttest,refdir);
  double refpdot = refpos.Dot(refdir);

  TCanvas* dpdotcan = new TCanvas("dpdot","P dot D derivatives",1200,800);
  dpdotcan->Divide(3,2);
  for(size_t ipar=0;ipar< KTLine::NParams();ipar++){
    KTLine::ParamIndex ip = static_cast<KTLine::ParamIndex>(ipar);
    posgraph[ipar] = new TGraph(ndel);
    string title = KTLine::paramName(ip) + string(" #Delta #vec{P} #bullet #vec{D};exact;1st derivative");
    posgraph[ipar]->SetTitle(title.c_str());
    KTLine xktline(ref_ktline);
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id;
      xktline.params().parameters()[ipar] = ref_ktline.params().parameters()[ipar]*(1.0 + delta);
      Vec3 pos;
      xktline.position(ttest,pos);
      double dpdot = pos.Dot(refdir) - refpdot;
      // linear approximation
      double dirdpdot = pder[ipar]*delta*ref_ktline.params().parameters()[ipar];
      posgraph[ipar]->SetPoint(id,dpdot,dirdpdot);
    }
    dpdotcan->cd(ipar+1);
    posgraph[ipar]->Draw("AC*");
  }
  dpdotcan->Write();

  ktline_deriv.Write();
  ktline_deriv.Close();
  return 0;
}
