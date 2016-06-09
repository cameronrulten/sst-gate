/*
The following code includes a set of functions for producing pretty plots
Written by Cameron Rulten 2013 Observatoire de Paris, France.
 */
#ifndef SST_GATE_TEST_PERFORMANCE_H
#define SST_GATE_TEST_PERFORMANCE_H

#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <fstream>
#include <numeric>

#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TPad.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TPolyLine3D.h"
#include "TGaxis.h"
#include "TArc.h"

#include "sst_gate_globals.hpp"
#include "sst_gate_class.hpp"
#include "sst_gate_psf.hpp"
#include "sst_gate_data_output.hpp"
#include "sst_gate_plotter.hpp"

using namespace std;

double calculateRMS(vector<double> &vTimes);
double calculateMEAN(vector<double> &vTimes);


class SST_GATE_PERFORMANCE{
public:
  SST_GATE_PERFORMANCE();
  ~SST_GATE_PERFORMANCE();
  
  void TestPerformance(const SST_GATE& gate, string method, string ofile_name, int set_focal_plane); //set focal plane to choose camera=0, primary=1, secondary=2
  //std::pair<double,double> get_encircled_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn);
  //std::pair<double,double> get_ensquared_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn, double sqHalfLength);
  //double getMean(std::vector<double> inArray);
  void check_effective_area(double N_focused[], int tN, string method_in, double area_in[], double field_angles[], double total_photons, TGraph* effAreaAltMethod);
};

SST_GATE_PERFORMANCE::SST_GATE_PERFORMANCE()
{
}

SST_GATE_PERFORMANCE::~SST_GATE_PERFORMANCE()
{
}

void SST_GATE_PERFORMANCE::TestPerformance(const SST_GATE& gate, string method, string ofile_name, int set_focal_plane)
{
  const int kN = 10;         // 0 to 4.5 [deg]
  //const double kDegStep = 0.5;
  double kDegStep = 0.5; //0.5
  double fDegStep = 0.4; //step for final angle. This is the analysis limit angle. 
  TH2D* histSpot[kN];
  TH1D* histTime[kN];
  TProfile* hTime_xFP[kN];         // photon position in x direction and mean time taken from secondary mirror
  TProfile* hTime_yFP[kN];         // photon position in y direction and mean time taken from secondary mirror
  TProfile2D* hTime_xyFP_S[kN];    // photon position on focal plane and mean time taken from secondary mirror
  TProfile2D* hTime_xyS[kN];       // photon position on secondary mirror and mean time taken from primary mirror
  TProfile2D* hTime_xyP[kN];       // photon position on primary mirror and mean time taken from primary mirror
  TProfile2D* hTime_xyFP_Total[kN];// photon position on focal plane and mean time taken from primary mirror
  TProfile2D* hTime_xyFP_Delay[kN];// photon position on focal plane and mean delay between photon arrival times
  vector<double> delay_x;
  vector<double> delay_y;
  vector<double> delay_t;
  double delay_mean_fov[kN];
  double delay_meanErr_fov[kN];
  double mints[kN];
  double maxts[kN];

  //Added by CBR
  TH1D* spotY_cm[kN];
  TH1D* spotX_cm[kN];
  TH2D* spotH_cm[kN];
  TH2D* spotH_deg[kN];
  TH2D* spotH_arcmin[kN];
  TH2D* spotGround_cm[kN];
  vector<TH2D*> spot_array;
  vector<double> x_array;
  vector<double> y_array;
  vector<double> z_array;
  vector< vector<double> > x_array_all;
  vector< vector<double> > y_array_all;
  vector< vector<double> > z_array_all;
  vector<double> x_array_ps;
  vector<double> y_array_ps;
  vector< vector<double> > x_array_all_ps;
  vector< vector<double> > y_array_all_ps;
 
  vector<double>::iterator x_array_iter;
  vector<double>::iterator y_array_iter;
  vector<double> psf_array, angle_array;
  vector<double> fov_array;
  double psf[kN];
  double psf_error[kN];
  double fov_angle[kN];
  double psf_radii[kN];
  double psf_photons[kN];
  double effectiveArea[kN];
  double angularResSagittal[kN];
  double angularResTangential[kN];
  vector<double> x_coordinates;
  vector<double> y_coordinates;
  vector<double> z_coordinates;
  vector<double> t_coordinates;
  vector< vector<double> > x_coordinates_at_points;
  vector< vector<double> > y_coordinates_at_points;
  vector< vector<double> > z_coordinates_at_points;
  vector< vector<double> > t_coordinates_at_points;
  double x_check[kN];
  double N_focused_array[kN];
  double area_array[kN];
  double totalRaysNumber;

  TH2D* spot_progression_cm = new TH2D("spot_progression_cm","Spot progression on-axis (0 degrees) to off-axis (4.5 degrees);X (cm);Y (cm);Number of photons",1000,-1.,1.04*gate.fRf,1000,-1.,1.); //full camera
  // // spot_progression_cm->SetMinimum(1.0);
  // for(int i=1;i<=1000;i++)
  //    {
  //      for(int j=1;j<=1000;j++)
  //  	{
  //  	  spot_progression_cm->SetBinContent(i,j,-1000);
  //  	  //spot_progression_cm->Fill(TRandom::Uniform(-1,1.1*gate.fRf),TRandom::Uniform(-1.,1.),-1.);
  //  	}
  //    }

  vector<TPolyLine3D*> poly_array;

  // photon position on secondary mirror and mean time taken from primary mirror
  //TH2D* hist_xyP= new TH2D("hist_xyP", "#it{#theta} = 0 (deg); x (cm); y (cm); number",1000,-210,210,1000,-210,210);
  //TH2D* hist_xyS= new TH2D("hist_xyS", "#it{#theta} = 0 (deg); x (cm); y (cm); number",1000,-110,110,1000,-110,110);
  TProfile2D* hist_xyP[kN];
  TProfile2D* hist_xyS[kN];

  map<double,double> mA[kN];
  stringstream getNum, getAngle;
  string histNum, histName, angleName;
  double totalRays[kN];  

  double starting_angle=0.0;//starting angle 0 degrees
  double theta = starting_angle;

  for(int n=0; n!=kN; ++n)
    {
      std::cout << " n: " << n << " processing angle: " << theta << std::endl;
      getNum << n;
      histNum = "spot_cm_" + getNum.str();
      getAngle << theta;
      histName = "#it{#theta} = " + getAngle.str() + " (deg);X (cm);Y (cm)";
      //spotH_cm[n] = new TH2D(histNum.c_str(),histName.c_str(),2500,-1.11*gate.fRf,1.11*gate.fRf,2500,-1.11*gate.fRf,1.11*gate.fRf); //scaled by 1.11 just so image isn't butted left
      histName = "#it{#theta} = " + getAngle.str() + " (deg);X (deg);Y (deg)";
      histNum = "spot_deg_" + getNum.str();
      spotH_deg[n] = new TH2D(histNum.c_str(),histName.c_str(),500,-5.0,5.0,500,-5.0,5.0);
      histName = "#it{#theta} = " + getAngle.str() + " (deg);X (arcmin);Y (arcmin)";
      histNum = "spot_arcmin_" + getNum.str();
      spotH_arcmin[n] = new TH2D(histNum.c_str(),histName.c_str(),500,-1.11*270.0,1.11*270,500,-1.11*270,1.11*270); //scaled by 1.11 just so image isn't butted left
      histNum = "spotX_cm_" + getNum.str();
      histName = "#it{#theta} = " + getAngle.str() + " (deg);Photon_{X} (cm);Number";
      spotX_cm[n] = new TH1D(histNum.c_str(),histName.c_str(),2500,-18.0,18.0);
      histNum = "spotY_cm_" + getNum.str();
      histName = "#it{#theta} = " + getAngle.str() + " (deg);Photon_{Y} (cm);Number";
      spotY_cm[n] = new TH1D(histNum.c_str(),histName.c_str(),2500,-18.0,18.0);
      

      getNum.str("");
      getAngle.str("");
      if(n==8) theta += fDegStep;//need to increment before n==9
      else theta += kDegStep;
    }
  
  //reset angles to 0
  //theta=0;
  //kDegStep=0.5;
  // for(int i=0; i<kN;i++)
  //   {
  //     spotH_cm[i]->Fill(0.,0.,-1.);
  //     spotH_deg[i]->Fill(0.,0.,-1.);
  //     spotH_arcmin[i]->Fill(0.,0.,-1.);
  //   }



  TGraph* graEffArea = new TGraph;
  graEffArea->SetName("graEffArea");
  graEffArea->SetTitle("Evolution of the effective collecting area with field angle;field angle (degrees);effective area (m^{2})");
  
  TGraph* graSigma[2];
  TLegend* leg = new TLegend(0.15, 0.65, 0.5, 0.85);
  leg->SetFillStyle(0);
  leg->SetTextFont(132);
  for(int i = 0; i < 2; i++){
    graSigma[i] = new TGraph;
    graSigma[i]->SetName(Form("graSigma%d", i));
    graSigma[i]->SetTitle(";Field angle (deg);Angular Resolution (arcmin)");
    graSigma[i]->SetLineStyle(i + 1);
    leg->AddEntry(graSigma[i], i == 0 ? "2 #times #sigma_{Sagittal}" : "2 #times #sigma_{Tangential}" , "lp");
  } // i
  
  TGraph* graTime = new TGraph;
  graTime->SetTitle(";Field angle (deg);Photon propagation time spread (#sigma) (ns)");

  string hname_new, htitle_new;  

  //reset angles
  theta=starting_angle;

  // TH2D* test_spot = new TH2D("test_spot", "#it{#theta} = 4 (degrees);X (cm);Y (cm)",2000,14.,18.,2000,-1.,1.);  
  TH2D* test_spot = new TH2D("test_spot", "#it{#theta} = 4 (degrees);X (cm);Y (cm)",2000,-1.,1.,2000,-1.,1.);//spot centre 0,0  
  std::vector<double> test_spot_vector_x;
  std::vector<double> test_spot_vector_y;

  std::cout << std::endl;
  std::cout << " start ray tracing " << std::endl;
  std::cout << "------------------" << std::endl;

  for(int i=0; i!=kN; ++i)
    {
      std::cout << " i: " << i << " processing angle: " << theta << std::endl;
      TObject* obj;
          
      obj = gROOT->Get(Form("histSpot%d", i));
     
      if(obj){
       	delete obj;
       	obj = 0;
      } // if
      
      histSpot[i] = new TH2D(Form("histSpot%d", i), Form("#it{#theta} = %.1f (deg);X (arcmin);Y (arcmin)", theta), 1000, -10, 10, 1000, -10, 10);
      
      obj = gROOT->Get(Form("spot_cm_%d", i)); 
      if(obj){
      	delete obj;
      	obj = 0;
      } // if
      
      spotH_cm[i] = new TH2D(Form("spot_cm_%d", i),Form("#it{#theta} = %.1f (degrees);X (cm);Y (cm)", theta),1000,-1.11*gate.fRf,1.11*gate.fRf,1000,-1.11*gate.fRf,1.11*gate.fRf); //full camera
			//spotH_cm[i] = new TH2D(Form("spot_cm_%d", i),Form("#it{#theta} = %.1f (degrees);X (cm);Y (cm)",theta),1000,-0.65,0.65,1000,-0.65,0.65);//spot centre 0,0
      //spotH_cm[i] = new TH2D(Form("spot_cm_%d", i),Form("#it{#theta} = %.1f (degrees);X (cm);Y (cm)", theta),2000,-0.3,2.1,2000,-0.3,2.1);//spot centre 0,0 5deg phi, pixels
      //spotH_cm[i] = new TH2D(Form("spot_cm_%d", i),Form("#it{#theta} = %.1f (deg);X (cm);Y (cm)", theta),1000,-2.5,2.5,1000,-2.5,2.5);//MAPMT    
      std::cout << "before spot: " << i << " title: " << spotH_cm[i]->GetTitle() << std::endl;

      obj = gROOT->Get(Form("spotGround_cm_%d", i)); 
      if(obj){
      	delete obj;
      	obj = 0;
      } // if
      
      spotGround_cm[i] = new TH2D(Form("spotGround_cm_%d", i),Form("#it{#theta} = %.1f (degrees);X (cm);Y (cm)",theta),1000,-300.,300.,1000,-300.,300.);//spot centre 0,0 

      

      obj = gROOT->Get(Form("histTime%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      histTime[i]= new TH1D(Form("histTime%d",i), Form("#it{#theta} = %.1f (deg);Propagation delay (ns);Entries", theta), 120, -6, 6);
      
      //Time histograms CR

      obj = gROOT->Get(Form("hTime_xFP%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position in x direction and mean time taken from secondary mirror
      hTime_xFP[i]= new TProfile(Form("hTime_xFP%d",i), Form("#it{#theta} = %.1f (deg); x (cm); time_{mean} (ns)", theta), 1000, -0.7, 0.7);

      obj = gROOT->Get(Form("hTime_yFP%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position in y direction and mean time taken from secondary mirror      
      hTime_yFP[i]= new TProfile(Form("hTime_yFP%d",i), Form("#it{#theta} = %.1f (deg); y (cm); time_{mean} (ns)", theta), 1000, -0.7, 0.7);//1,3

      obj = gROOT->Get(Form("hTime_xyFP_S%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on focal plane and mean time taken from secondary mirror
      hTime_xyFP_S[i]= new TProfile2D(Form("hTime_xyFP_S%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); time_{mean} (ns)", theta),1000,-0.7,0.7,1000,-0.7,0.7);

      // obj = gROOT->Get(Form("hTime_xyFP_Total%d", i));
      // if(obj){
      // 	delete obj;
      // 	obj = 0;
      // } // if


      obj = gROOT->Get(Form("hTime_xyS%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on secondary mirror and mean time taken from primary mirror
      hTime_xyS[i]= new TProfile2D(Form("hTime_xyS%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); time_{mean} (ns)", theta),1000,-110,110,1000,-110,110);

      obj = gROOT->Get(Form("hTime_xyP%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on secondary mirror and mean time taken from primary mirror
      hTime_xyP[i]= new TProfile2D(Form("hTime_xyP%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); time_{mean} (ns)", theta),1000,-210,210,1000,-210,210);

      obj = gROOT->Get(Form("hTime_xyFP_Total%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on focal plane and mean time taken from primary mirror
      hTime_xyFP_Total[i]= new TProfile2D(Form("hTime_xyFP_Total%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); time_{mean} (ns)", theta),1000,-0.7,0.7,1000,-0.7,0.7);//centred on zero
      //hTime_xyFP_Total[i]= new TProfile2D(Form("hTime_xyFP_Total%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); time_{mean} (ns)", theta),1000,-18.0,+18.0,1000,-18.0,+18.0);//full camera

      obj = gROOT->Get(Form("hTime_xyFP_Delay%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on focal plane and mean delay between photon arrival times
      hTime_xyFP_Delay[i]= new TProfile2D(Form("hTime_xyFP_Delay%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); delay_{mean} (ns)", theta),1000,-0.7,0.7,1000,-0.7,0.7);

      obj = gROOT->Get(Form("hist_xyS%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on secondary mirror and mean time taken from primary mirror
      hist_xyS[i]= new TProfile2D(Form("hist_xyS%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); z (cm)", theta),1000,-110,110,1000,-110,110);

      obj = gROOT->Get(Form("hist_xyP%d", i));
      if(obj){
	delete obj;
	obj = 0;
      } // if

      // photon position on secondary mirror and mean time taken from primary mirror
      hist_xyP[i]= new TProfile2D(Form("hist_xyP%d",i), Form("#it{#theta} = %.1f (deg); x (cm); y (cm); z (cm)", theta),1000,-210,210,1000,-210,210);

      //-------------------------------
      // Ray tracing section
      //-------------------------------
      

      const double kZs = gate.fF/gate.fQ;                 // fF/fQ = separation distance between primary and secondary
      const double kZf = kZs - (1 - gate.fAlpha) * gate.fF; // (1 - fAlpha)*f = Fseparation distance between secondary & focal plane
      double lambda = 400*nm;
      TRandom3* ran3 = new TRandom3();
      ran3->SetSeed(static_cast<long unsigned int>(gSystem->Now()));
      //ran3->SetSeed(0);

      TVector3 dir;
      ARayArray* array;
      double showerAlt = 10.0;//[km]
      double scatAreaInc = 2.2; //increase the scattering area
      double rho = scatAreaInc * kZs;
      double phi = 0.0;
      //deg=theta
      if(method=="square")
	{
	  //theta only translation & rotation:
	  //TGeoTranslation raytr("raytr", -1.5*kZs*TMath::Sin(theta*TMath::DegToRad()), 0, 1.5*kZs*TMath::Cos(theta*TMath::DegToRad()));  //for x-direction
	  //TGeoTranslation raytr("raytr", 0, -1.5*kZs*TMath::Sin(theta*TMath::DegToRad()), 1.5*kZs*TMath::Cos(theta*TMath::DegToRad()));  //for y-direction
	  //TGeoTranslation raytr("raytr", 1.5*kZs*TMath::Sin(theta*TMath::DegToRad()), -1.5*kZs*TMath::Sin(theta*TMath::DegToRad()), 1.5*kZs*TMath::Cos(theta*TMath::DegToRad()));//for xy-direction
	  //dir.SetMagThetaPhi(1, TMath::Pi() - theta*TMath::DegToRad(), 0);//theta only //for +y dir: +trans, +theta and +90 in phi //for -y dir: -trans, -theta and +90 in phi

	  //theta and phi translation & rotation:
	  TGeoTranslation raytr("raytr", -rho*TMath::Sin(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad()), -rho*TMath::Sin(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad()), rho*TMath::Cos(theta*TMath::DegToRad())); //spherical
	  dir.SetMagThetaPhi(1, TMath::Pi() - theta*TMath::DegToRad(), phi*TMath::DegToRad());//spherical
	  
	  //define ray shooter
	  array = ARayShooter::Square(lambda, scatAreaInc*gate.fRpMax, 401, 0, &raytr, &dir);// 1 photon per 0.0025 m^2
	}
      else if(method=="cone")
	{
	  TGeoTranslation raytr("raytr", -showerAlt*km*TMath::Sin(theta*TMath::DegToRad()), 0, showerAlt*km*TMath::Cos(theta*TMath::DegToRad()));
	  //cout << " km: " << 10*km << endl;
	  TGeoRotation rotShoot("rotShoot", 90, 180 - theta, 0);
	  array = ARayShooter::RandomCone(lambda, 1.1*gate.fRpMax, showerAlt*km, 1.52053e5, &rotShoot, &raytr);    
	}
      else
	{
	  cout << "You have not provided a ray tracing method! e.g. square (for parallel rays) or random cone (for height specific)" << endl;
	  exit(1);
	}
      
      gate.GetManager()->TraceNonSequential(*array);

      TObjArray* focused = array->GetFocused();

      //new to get total photons for efective area:
      TObjArray* stopped = array->GetStopped();
      TObjArray* suspended = array->GetSuspended();
      TObjArray* exited = array->GetExited();
      TObjArray* focused2 = array->GetFocused();

      totalRays[i] = (stopped->GetEntries() + focused2->GetEntries() + exited->GetEntries() + suspended->GetEntries());// - (exited->GetEntries());
      totalRaysNumber = (stopped->GetEntries() + focused2->GetEntries() + exited->GetEntries() + suspended->GetEntries());// - (exited->GetEntries());

      cout << "field angle: " << theta << " totalRaysNumber: " << totalRaysNumber << " stopped: " << stopped->GetEntries() << " focused: " << focused2->GetEntries() << " exited: " << exited->GetEntries() << " suspended: " << suspended->GetEntries() << endl;

      // Check that rays have actually focused before contuining with the ray trace.
      // Really important when using the actual focal plane!!!
      if(focused->GetLast()<=0) continue;

      double Aeff = 0.;
      double Aeff_alt = 0.;
      double mypoints[focused->GetLast()];
      int countPoints=0;
      int countPoints2=0;
      double Aeff2 = 0.;

      //double N_focused=0.;

      N_focused_array[i] = focused->GetLast();
      if(method=="square") area_array[i] = (scatAreaInc*gate.fRpMax*1.e-2) * (scatAreaInc*gate.fRpMax*1.e-2); //convert from cm to m by multiplying by 10^-2
      else if(method=="cone") area_array[i] = TMath::Cos(theta*TMath::DegToRad())*TMath::Pi()*TMath::Power( (1.1*gate.fRpMax*1.e-2) ,2);//convert from cm to m by multiplying by 10^-2 //

      for(int j=0; j<=focused->GetLast(); j++)
	{

	  ARay* ray = (ARay*)(*focused)[j];

	  if(!ray) continue;

	  // Calculate the effective area from the number of focused photons
	  //if(method=="square") Aeff += scatAreaInc*gate.fRpMax*scatAreaInc*gate.fRpMax/400./400./m/m;
	  if(method=="square") Aeff += scatAreaInc*gate.fRpMax*scatAreaInc*gate.fRpMax/160000./m/m; //totalRaysNumber;
	  if(method=="square") Aeff2 += scatAreaInc*gate.fRpMax*scatAreaInc*gate.fRpMax/totalRaysNumber/m/m; //totalRaysNumber;
	  else if(method=="cone") Aeff2 += TMath::Cos(theta*TMath::DegToRad())*TMath::Pi()*TMath::Power(1.1*gate.fRpMax,2)/1.52053e5/m/m;
	  //N_focused++;TMath::Cos(theta*TMath::DegToRad())*

	  double *p;
	  if(method=="square")
	    {
	      switch(set_focal_plane)
		{
		case 0:
		  p = new double[4];
		  if(ray->GetNpoints()<4 || ray->GetNpoints()>4) cout << "Ray points (Camera as FP): " << ray->GetNpoints() << endl;
		  break;
		case 1:
		  p = new double[3];
		  if(ray->GetNpoints()<2 || ray->GetNpoints()>4) cout << "Ray points (Primary as FP): " << ray->GetNpoints() << endl;
		  break;
		case 2:
		  p = new double[3];
		  if(ray->GetNpoints()<3 || ray->GetNpoints()>4) cout << "Ray points (Secondary as FP): " << ray->GetNpoints() << endl;
		  break;
		}
	    }
	  else if(method=="cone")
	    {
	      p = new double[5];
	      if(ray->GetNpoints()<5 || ray->GetNpoints()>5) cout << "Ray points (Method cone): " << ray->GetNpoints() << endl;
	    }
	  
	  
	  //p = new double[ray->GetNpoints()];
	  ray->GetLastPoint(p);
	  ray->SetLineWidth(1);

	  //To draw polyline rays uncomment this section
	  if(theta == 0 && ran3->Uniform(1) < 0.001)//j==countPoints2)//ran3->Uniform(1) < 0.001) //gRandom->Uniform(1)
	    { 
	      TPolyLine3D* pol = ray->MakePolyLine3D();
	      pol->SetLineColor(4);
	      pol->SetLineStyle(2);
	      poly_array.push_back(pol);
	      //pol->Draw();
	      countPoints2+=1010;
	    } // if
	  
	  double x = theta*gate.deg2dist;
	  x_check[i] = x;
	  histSpot[i]->Fill((p[0] - x)/gate.deg2dist*60, p[1]/gate.deg2dist*60);
	  histTime[i]->Fill((p[3] - (1.5*kZs - kZf)/(TMath::C()*m))/1e-9); // ns

	  double px0,py0,pz0, px1,py1,pz1, px2,py2,pz2, px3,py3,pz3, px4,py4,pz4;
	  double pt0,pt1,pt2,pt3,pt4;

	  if(method=="square")
	    {
	      switch(set_focal_plane)
		{
		case 0: //camera
		  ray->GetPoint(0,px0,py0,pz0,pt0); // x-coord, y-coord, z-coord, t-coord at start point
		  ray->GetPoint(1,px1,py1,pz1,pt1); // x-coord, y-coord, z-coord, t-coord at primary mirror
		  ray->GetPoint(2,px2,py2,pz2,pt2); // x-coord, y-coord, z-coord, t-coord at secondary mirror
		  ray->GetPoint(3,px3,py3,pz3,pt3); // x-coord, y-coord, z-coord, t-coord at focal plane
		  
		  //Fill histograms for plotting etc.
		  countPoints++;
		  hTime_xFP[i]->Fill(px3-x,(pt3-pt2)*1.e9);           // photon position in x direction and mean time taken from secondary mirror
		  hTime_yFP[i]->Fill(py3,(pt3-pt2)*1.e9);             // photon position in y direction and mean time taken from secondary mirror
		  hTime_xyFP_S[i]->Fill(px3-x,py3,(pt3-pt2)*1.e9);    // photon position on focal plane and mean time taken from secondary mirror
		  hTime_xyP[i]->Fill(px1,py1,(pt1-pt0)*1.e9);         // photon position on secondary mirror and mean time taken from primary mirror
		  hTime_xyFP_Total[i]->Fill(px3-x,py3,(pt3-pt1)*1.e9);// centred on zero //photon position on focal plane and mean time taken from primary mirror	      
		  hTime_xyS[i]->Fill(px2,py2,(pt2-pt1)*1.e9);         // photon position on secondary mirror and mean time taken from primary mirror
		  //hTime_xyFP_Total[i]->Fill(px3,py3,(pt3-pt1)*1.e9); // full camera //photon position on focal plane and mean time taken from primary mirror
		  hist_xyP[i]->Fill(px1,py1,pz1);
		  hist_xyS[i]->Fill(px2,py2,pz2);
		  delay_x.push_back(px3-x);
		  delay_y.push_back(py3);
		  delay_t.push_back(1.e9*(pt3-pt1));//nanseconds
		  //if(i==1)
		  //cout << "pt0: " << pt0*1.e9 << " pt1: " << pt1*1.e9 << " pt2: " << pt2*1.e9 << " pt3: " << pt3*1.e9 << " M2-M1: " << TMath::C()*(pt2-pt1) << " m FP-M2: " << TMath::C()*(pt3-pt2) << " m FP-M1: " << TMath::C()*(pt3-pt1) << " m " << endl;//distances
		  //cout << i << "" << j << " pt0: " << pt0*1.e9 << " pt1: " << pt1*1.e9 << " pt2: " << pt2*1.e9 << " pt3: " << pt3*1.e9 << " M2-M1: " << 1.e9*(pt2-pt1) << " ns FP-M2: " << 1.e9*(pt3-pt2) << " ns FP-M1: " << 1.e9*(pt3-pt1) << " ns " << endl;//nanseconds
		  //cout << "t0: " << pt0/TMath::C()*m/1.e-9 << " t1: " << pt1/TMath::C()*m/1.e-9 << " t2: " << pt2/TMath::C()*m/1.e-9 << " t3: " << pt3/TMath::C()*m/1.e-9 << "t3-t1: " << (pt3/TMath::C()*m/1.e-9)-(pt1/TMath::C()*m/1.e-9) << endl;
		  break;
		case 1: //primary
		  ray->GetPoint(0,px0,py0,pz0,pt0); // x-coord, y-coord, z-coord, t-coord at start point
		  ray->GetPoint(1,px1,py1,pz1,pt1); // x-coord, y-coord, z-coord, t-coord at primary mirror
		  hist_xyP[i]->Fill(px1,py1,pz1);
		  hTime_xyP[i]->Fill(px1,py1,(pt1-pt0)*1.e9);
		  delay_x.push_back(px1-x);
		  delay_y.push_back(py1);
		  delay_t.push_back(1.e9*(pt1-pt0));//nanseconds
		  break;
		case 2: //secondary
		  ray->GetPoint(0,px0,py0,pz0,pt0); // x-coord, y-coord, z-coord, t-coord at start point
		  ray->GetPoint(1,px1,py1,pz1,pt1); // x-coord, y-coord, z-coord, t-coord at primary mirror
		  ray->GetPoint(2,px2,py2,pz2,pt2); // x-coord, y-coord, z-coord, t-coord at secondary mirror
		  hist_xyS[i]->Fill(px2,py2,pz2);
		  hTime_xyS[i]->Fill(px2,py2,(pt2-pt1)*1.e9);
		  delay_x.push_back(px2-x);
		  delay_y.push_back(py2);
		  delay_t.push_back(1.e9*(pt2-pt1));//nanseconds
		  break;
		}
	      
	    }
	  if(method=="cone")
	    {
	      ray->GetPoint(0,px0,py0,pz0,pt0); // x-coord, y-coord, z-coord, t-coord at start point
	      ray->GetPoint(1,px1,py1,pz1,pt1); // x-coord, y-coord, z-coord, t-coord at start point
	      ray->GetPoint(2,px2,py2,pz2,pt2); // x-coord, y-coord, z-coord, t-coord at primary mirror
	      ray->GetPoint(3,px3,py3,pz3,pt3); // x-coord, y-coord, z-coord, t-coord at secondary mirror
	      ray->GetPoint(4,px4,py4,pz4,pt4); // x-coord, y-coord, z-coord, t-coord at focal plane
	      
	      //fill histograms for plotting etc.
	      countPoints++;
	      //ray->GetPoint(4,px4,py4,pz4,pt4);
	      hTime_xFP[i]->Fill(px4-x,(pt4-pt3)*1.e9);           // photon position in x direction and mean time taken from secondary mirror
	      hTime_yFP[i]->Fill(py4,(pt4-pt3)*1.e9);             // photon position in y direction and mean time taken from secondary mirror
	      hTime_xyFP_S[i]->Fill(px4-x,py4,(pt4-pt3)*1.e9);    // photon position on focal plane and mean time taken from secondary mirror
	      hTime_xyP[i]->Fill(px2,py2,(pt2-pt0)*1.e9);  // photon position on primary mirror and mean time taken from primary mirror
	      hTime_xyS[i]->Fill(px3,py3,(pt3-pt2)*1.e9);  // photon position on secondary mirror and mean time taken from primary mirror
	      //hTime_xyFP_Total[i]->Fill(px4,py4,(pt4-pt2)*1.e9); // full camera //photon position on focal plane and mean time taken from primary mirror	   
	      //hTime_xyFP_Total[i]->Fill(px4-x,py4,(pt4-pt2)*1.e9);// centred on zero //photon position on focal plane and mean time taken from primary mirror

	      hist_xyP[i]->Fill(px2,py2,pz2);//);
	      hist_xyS[i]->Fill(px3,py3,pz3);
	      delay_x.push_back(px4-x);
	      delay_y.push_back(py4);
	      delay_t.push_back((pt4-pt2)*1.e9);
	    }
	  
	  // delay_x.push_back(px3-x);
	  // delay_y.push_back(py3);
	  // delay_t.push_back((pt3-pt1)*1.e9);
 
	  //if(method=="square") cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << ray->GetNpoints() << endl;
	  //if(method=="cone")  cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << " " << ray->GetNpoints() << endl;

	  // if(method=="square")
	  //   {
	  //     countPoints++;
	  //     hTime_xFP[i]->Fill(px3-x,(pt3-pt2)*1.e9);           // photon position in x direction and mean time taken from secondary mirror
	  //     hTime_yFP[i]->Fill(py3,(pt3-pt2)*1.e9);             // photon position in y direction and mean time taken from secondary mirror
	  //     hTime_xyFP_S[i]->Fill(px3-x,py3,(pt3-pt2)*1.e9);    // photon position on focal plane and mean time taken from secondary mirror
	  //     hTime_xyP[i]->Fill(px1,py1,(pt1-pt0)*1.e9);         // photon position on secondary mirror and mean time taken from primary mirror
	  //     hTime_xyFP_Total[i]->Fill(px3-x,py3,(pt3-pt1)*1.e9);// centred on zero //photon position on focal plane and mean time taken from primary mirror	      
	  //     hTime_xyS[i]->Fill(px2,py2,(pt2-pt1)*1.e9);         // photon position on secondary mirror and mean time taken from primary mirror
	  //     //hTime_xyFP_Total[i]->Fill(px3,py3,(pt3-pt1)*1.e9); // full camera //photon position on focal plane and mean time taken from primary mirror
	  //     hist_xyP[i]->Fill(px1,py1,pz1);
	  //     hist_xyS[i]->Fill(px2,py2,pz2);
	  //   }
	  // else if(method=="cone")
	  //   {
	  //     countPoints++;
	  //     //ray->GetPoint(4,px4,py4,pz4,pt4);
	  //     hTime_xFP[i]->Fill(px3-x,(pt3-pt2)*1.e9);           // photon position in x direction and mean time taken from secondary mirror
	  //     hTime_yFP[i]->Fill(py3,(pt3-pt2)*1.e9);             // photon position in y direction and mean time taken from secondary mirror
	  //     hTime_xyFP_S[i]->Fill(px3-x,py3,(pt3-pt2)*1.e9);    // photon position on focal plane and mean time taken from secondary mirror
	  //     hTime_xyP[i]->Fill(px2,py2,(pt2-pt0)*1.e9);  // photon position on primary mirror and mean time taken from primary mirror
	  //     hTime_xyS[i]->Fill(px3,py3,(pt3-pt2)*1.e9);  // photon position on secondary mirror and mean time taken from primary mirror
	  //     //hTime_xyFP_Total[i]->Fill(px4,py4,(pt4-pt2)*1.e9); // full camera //photon position on focal plane and mean time taken from primary mirror	   
	  //     //hTime_xyFP_Total[i]->Fill(px4-x,py4,(pt4-pt2)*1.e9);// centred on zero //photon position on focal plane and mean time taken from primary mirror

	  //     hist_xyP[i]->Fill(px2,py2,pz2);//);
	  //     hist_xyS[i]->Fill(px3,py3,pz3);   
	  //   }
	  //cout << "t0: " << p[0]/TMath::C()*m/1.e-9 << " t1: " << p[1]/TMath::C()*m/1.e-9 << " t2: " << p[2]/TMath::C()*m/1.e-9 << " t3: " << p[3]/TMath::C()*m/1.e-9 << endl;

	  double fillX = p[0];  // (p[0] - (x/gate.deg2dist)); //cm

	  //cout << theta << " deg " << p[0] << " cm " << p[0]/gate.deg2dist << " deg " << x << " cm " << x/gate.deg2dist << " deg " << p[0]-x << " cm " <<  endl;
	  //cout << " p0: " << pt0*1.e9 << " p1: " << pt1*1.e9 << " p2: " << pt2*1.e9 << " p3: " << pt3*1.e9 << "pt2-pt1: " << (pt2-pt1)*1.e9 << " pt3-pt2: " << (pt3-pt2)*1.e9 << endl;

	  double fillY = p[1]; //cm
	  spotX_cm[i]->Fill(p[0]);
	  spotY_cm[i]->Fill(p[1]);
	  spotH_cm[i]->Fill(p[0],p[1]);//spot centre 0,0 except 4.5 degrees.
	  spotH_deg[i]->Fill(p[0]/gate.deg2dist,p[1]/gate.deg2dist);
	  spotH_arcmin[i]->Fill(p[0]/gate.deg2dist*60.0,p[1]/gate.deg2dist*60.0);
	  spot_progression_cm->Fill(p[0],p[1]);
	  spotGround_cm[i]->Fill(p[0],p[1]);
	  if(i==7)
	    {
	      test_spot->Fill(p[0],p[1]);
	      test_spot_vector_x.push_back(p[0]);
	      test_spot_vector_y.push_back(p[1]);
	    }
	  x_array.push_back(p[0]);
	  y_array.push_back(p[1]);
	  z_array.push_back(p[2]);
	  x_array_ps.push_back(p[0]-x);
	  y_array_ps.push_back(p[1]);
	} // j

      cout << "totalRaysNumber: " << totalRaysNumber << " Aeff: " << Aeff << " Aeff2: " << Aeff2 << endl;

      vector<double> dta;
      
      double minT = (*min_element(delay_t.begin(),delay_t.end()));
      double maxT = (*max_element(delay_t.begin(),delay_t.end()));
      mints[i] = minT;
      maxts[i] = maxT;

      if(minT<=0) cout << i << " minT: " << minT << " maxT: " << maxT << endl;

      for(int j = 0; j <= focused->GetLast(); j++)
      	{
      	  hTime_xyFP_Delay[i]->Fill(delay_x[j], delay_y[j], delay_t[j]-minT);
	  dta.push_back(delay_t[j] - minT);
      	}
      
      // if(i==1)
      //  	{
      //  	  for(int j=0; j<=focused->GetLast(); j++)
      //  	    {
      //  	      cout << "time delay: " << delay_t[j] << " " << minT << " " << delay_t[j]-minT << " " <<  maxT << endl;
      //  	    }
      // 	}



      delay_mean_fov[i] = TMath::Mean(focused->GetLast(), &dta[0], 0);
      delay_meanErr_fov[i] =TMath::RMS(dta.begin(), dta.end());// TMath::RMS(focused->GetLast(), &dta[0]); //
      //delay_meanErr_fov[i] = calculateRMS(dta);
      
      cout << " TIME mean: " << TMath::Mean(focused->GetLast(), &dta[0]) << " rms: " << TMath::RMS(dta.begin(), dta.end()) << " rms2: " << TMath::RMS(focused->GetLast(), &dta[0]) << " rms3: " << calculateRMS(dta) << " mean2: " << calculateMEAN(dta) << endl;

      delay_x.clear();
      delay_y.clear();
      delay_t.clear();


      effectiveArea[i] = Aeff2;
      graEffArea->SetPoint(graEffArea->GetN(), theta, Aeff2);//original
      
      double rmsx = histSpot[i]->GetRMS(1);
      double rmsy = histSpot[i]->GetRMS(2);
      //cout << "RMSX: " << rmsx << " RMSXerr: " << TMath::Sqrt(histSpot[i]->GetStdDev(1))/histSpot[i]->GetEntries() << " RMSY: " << rmsy << " RMSYerr: " << TMath::Sqrt(histSpot[i]->GetStdDev(2))/histSpot[i]->GetEntries() << endl;
      graSigma[0]->SetPoint(graSigma[0]->GetN(), theta, 2*rmsx);
      graSigma[1]->SetPoint(graSigma[1]->GetN(), theta, 2*rmsy);
      
      angularResSagittal[i] = 2.*rmsx;
      angularResTangential[i] = 2.*rmsy;

      graTime->SetPoint(graTime->GetN(), theta, histTime[i]->GetRMS());
      std::cerr << "angle number: " << i << " theta: " << theta <<  " x-size: " << x_array.size() << " y-size: " << y_array.size() << " z-size: " << z_array.size() << std::endl;
      std::cout << "------------------" << std::endl;
      fov_array.push_back(theta); //only need to store the fov angles once
      fov_angle[i] = theta;
      x_array_all.push_back(x_array);
      y_array_all.push_back(y_array);
      z_array_all.push_back(z_array);
      x_array_all_ps.push_back(x_array_ps);
      y_array_all_ps.push_back(y_array_ps);

      x_array.clear();
      y_array.clear();
      z_array.clear();
      x_array_ps.clear();
      y_array_ps.clear();
      if(i==8) theta+=fDegStep;//need to increment befor i=9
      else theta+=kDegStep;
      //theta+=0.5;

      delete array;
           
    } // i
  
  //----------------------------------------------------------------------------------------------
  //----END OF RAY TRACING SECTION----------------------------------------------------------------
  //----------------------------------------------------------------------------------------------

  //Test effective area
  check_effective_area(N_focused_array, kN, method, area_array, fov_angle, totalRaysNumber, graEffArea);


  //Centre spot data and fill spot histograms
  //-----centre all spot images-----
  int countV=0;
  std::vector< std::vector<double> > x_coord_all_centred;
  std::vector< std::vector<double> > y_coord_all_centred;
  std::vector<double> x_coord_centred;
  std::vector<double> y_coord_centred;

  //The next 2 lines are just a test of the function in the PSF class which produces the same
  //mean and psf results derived below. This is safe to use for PSF comparisons of sim_telarray data.
  //SST_GATE_PSF* getPSF = new SST_GATE_PSF();
  //getPSF->centre_coordinates(x_array_all, y_array_all, x_coord_all_centred, y_coord_all_centred, 0, kN);

  
  TGraph* gr_psf[kN];

  cout << " x-size: " << x_array_all.size() << " y-size: " << y_array_all.size() << endl;
  std::vector< vector<double> >::iterator x_array_all_iter;
  std::vector< vector<double> >::iterator y_array_all_iter=y_array_all.begin();

  for(x_array_all_iter=x_array_all.begin();x_array_all_iter!=x_array_all.end();++x_array_all_iter)
    {

      cout << " x-s: " << (*x_array_all_iter).size() << " y-s: " << (*y_array_all_iter).size() << endl;


      double sumX = std::accumulate((*x_array_all_iter).begin(),(*x_array_all_iter).end(),0.0);
      double meanX = sumX / (*x_array_all_iter).size();

      double sumY = std::accumulate((*y_array_all_iter).begin(),(*y_array_all_iter).end(),0.0);
      double meanY = sumY / (*y_array_all_iter).size(); 

      vector<double>::iterator x_coord_array_iter;
      vector<double>::iterator y_coord_array_iter = (*y_array_all_iter).begin();
      for( x_coord_array_iter=(*x_array_all_iter).begin(); x_coord_array_iter!=(*x_array_all_iter).end(); ++x_coord_array_iter)
	{
	  double getX = (*x_coord_array_iter);
	  double getY = (*y_coord_array_iter);
	  double getX_psc = getX-x_check[countV]; //platescale_centred
	  double getX_mc = getX-meanX;            //mean_centred

	  //mean centred spot
	  spotH_cm[countV]->Fill(getX-meanX, getY);
	  x_coord_centred.push_back(getX_mc);
	  y_coord_centred.push_back(getY);
	  //platescale centred spot
	  //spotH_cm[countV]->Fill(getX_psc, getY);
	  //x_coord_centred.push_back(getX_psc);
	  //y_coord_centred.push_back(getY);

	  //if(countV==9) cout << "num: " << countV << " x: " << getX << " y: " << getY << " x-<x>: " << getX-meanX << " y-<y>: " << getY-meanY << " meanX: " << meanX << " ps: " << (fov_angle[countV]*gate.deg2dist) << " x-check: " << x_check[countV] << " getX-ps: " << getX-(fov_angle[countV]*gate.deg2dist) << " fov: " << fov_angle[countV] << endl;
	  y_coord_array_iter++;
	}

      gr_psf[countV] = new TGraph(x_coord_centred.size(), &x_coord_centred[0], &y_coord_centred[0]);
      x_coord_all_centred.push_back(x_coord_centred);
      y_coord_all_centred.push_back(y_coord_centred);
      x_coord_centred.clear();
      y_coord_centred.clear();

      countV++;
      y_array_all_iter++;
    }
  
  cout << " x-all-centred-size: " << x_coord_all_centred.size() << " y-all-centred-size: " << y_coord_all_centred.size() << endl;
  
  for(int i=0;i!=kN;++i)
    {
      cout << "fov: " << fov_angle[i] << " mean delay: " << delay_mean_fov[i] << " spread: " << delay_meanErr_fov[i] << " min: " << mints[i] << " max: " << maxts[i] << " maxdiff: " << maxts[i]-mints[i] << endl;

      cout << "Total Rays: Field Angle: " << endl;
      cout << totalRays[i] << " " << fov_angle[i] << endl;

    }
  
  //----------------------------------------------------------------------------------------------
  // Check that the vectors are filled with the correct values.
  //----------------------------------------------------------------------------------------------

  // std::vector<double>::iterator x_coord_centred_iter;
  // for(int i=9; i!=10; ++i)
  //   {
  //     std::vector<double>::iterator y_coord_centred_iter = y_coord_all_centred[i].begin();
  //     std::vector<double>::iterator x_coord_ps_iter = x_array_all_ps[i].begin();
  //     std::vector<double>::iterator y_coord_ps_iter = y_array_all_ps[i].begin();
  //     for(x_coord_centred_iter=x_coord_all_centred[i].begin();x_coord_centred_iter!=x_coord_all_centred[i].end();++x_coord_centred_iter)
  // 	{
  // 	  cout << " x-c: " << (*x_coord_centred_iter) << " x-ps: " << (*x_coord_ps_iter) << " y-c: " << (*y_coord_centred_iter) << " y-ps: " << (*y_coord_ps_iter) << endl;
  // 	  x_coord_ps_iter++;
  // 	  y_coord_ps_iter++;
  // 	  y_coord_centred_iter++;
  // 	}
  //   }
  
  
  //----------------------------------------------------------------------------------------------
  // Calculate PSF using centred spots
  //----------------------------------------------------------------------------------------------

  
  SST_GATE_PSF* getPSF = new SST_GATE_PSF();

  bool derive_psf=true;
  int centred_method=0; //mean-centred method = 0 or platescale-centred method = 1. Any other integer value will default to the mean-centred method.

  vector< vector<double> > psf_encircled_all;
  vector< vector<double> > psf_ensquared_all;
  vector< vector<double> > psf_encircled_all_ps;
  vector< vector<double> > psf_ensquared_all_ps;
  
  SST_GATE_DATA_OUTPUT* output_data2file = new SST_GATE_DATA_OUTPUT();

  if(set_focal_plane==0 and derive_psf==true)
    {
      //must provide centred spot coordinates in order to calculate the PSF.      
      getPSF->calculate_psf(gate.deg2dist, x_coord_all_centred, y_coord_all_centred, fov_array, method, ofile_name, psf_encircled_all, psf_ensquared_all, centred_method);
      
      cout << "encircled array size: " << psf_encircled_all.size() << " ensquared array size: " << psf_ensquared_all.size() << endl;
      
      for(int i=0;i!=psf_encircled_all.size();++i)
	{
	  cout << " psf encircled array sizes: " << psf_encircled_all[i].size() << " array#: " << i << endl;
	}
      for(int i=0;i!=psf_ensquared_all.size();++i)
	{
	  cout << " psf ensquared array sizes: " << psf_ensquared_all[i].size() << " array#: " << i << endl;
	}
      
      // //-----test between methods-----
       // getPSF->calculate_psf(gate.deg2dist, x_array_all_ps, y_array_all_ps, fov_array, method, ofile_name, psf_encircled_all_ps, psf_ensquared_all_ps, platescale_method);

      // for(int i=0;i!=10;++i)
      //  	{
      //  	  cout << "encircled ps: " << psf_encircled_all_ps[0][i] << " centred: " << psf_encircled_all[0][i] << endl;
      //  	}

      
      output_data2file->write_data(psf_encircled_all, psf_ensquared_all, fov_array, spotH_cm, effectiveArea, angularResSagittal, angularResTangential, gate.fRpMax, method, ofile_name);
 
      //<<<REMOVED SECTION 1 March 2014>>>
      
      //---------------------------------------------------------------------------------------
      // New plotting methods
      //---------------------------------------------------------------------------------------
      
      // //-----Plot centred spots-----
      // TCanvas* spot_canvas33 = new TCanvas("spot_canvas33", "spot_canvas33", 1440, 500);
      // // gStyle->SetOptStat(0);
      // gStyle->SetTitleFontSize(0.06);
      // //TCanvas* spot_canvas33 = new TCanvas("spot_canvas33", "spot_canvas33", 3508, 2480); //optimised for A4 landscape 300dpi
      // spot_canvas33->Divide(5, 2);
      // TPad *spot_pads[kN];
      // std::stringstream get_padNum;
      // std::string padNum;
      // TPaveStats *mystats[kN];
      // TPaveStats *newStats;
      // TH2D* spotH_cm_clone2[kN];
      
      // for(int i=0;i!=kN;++i)
      // 	{
      // 	  spot_pads[i] = (TPad*)spot_canvas33->GetPad(i+1);
      // 	  spotH_cm_clone2[i] = (TH2D*)spotH_cm[i]->Clone("spotH_cm_clone2");
      // 	}
      
      // for(int i=0;i!=kN;++i)
      // 	{
      // 	  spot_pads[i]->SetRightMargin(0.15);
      // 	  spot_pads[i]->Draw();
      // 	}
      
      // spot_canvas33->Update();
      
      // //theta=0.0; //reset angle
      // //photonsFile << "# FOVangle(deg): \t photons: \t photons: \t psf(cm): \t effArea(m2): \t angResS(arcmin): \t angResT(arcmin): \t M1radius(cm): \t method:" << std::endl;
      // for(int i=0; i!=kN; ++i)
      // 	{
      // 	  //change pad
      // 	  spot_pads[i]->cd();
      // 	  //plot unbinned spots
      // 	  // gr_psf[i]->GetXaxis()->SetLimits(-0.6,0.6);
      // 	  // gr_psf[i]->GetYaxis()->SetRangeUser(-0.6,0.6);
      // 	  // gr_psf[i]->SetTitle("; xcoord (cm); ycoord (cm)");
      // 	  // gr_psf[i]->GetXaxis()->SetLabelSize(0.035);
      // 	  // gr_psf[i]->GetYaxis()->SetLabelSize(0.035);
      // 	  // gr_psf[i]->GetYaxis()->SetTitleSize(0.03);
      // 	  // gr_psf[i]->GetXaxis()->SetTitleSize(0.03);
      // 	  // gr_psf[i]->GetXaxis()->SetTitleOffset(1.3);
      // 	  // gr_psf[i]->GetYaxis()->SetTitleOffset(1.6);
      // 	  // gr_psf[i]->SetLineColor(1);
      // 	  // gr_psf[i]->SetMarkerStyle(6);
      // 	  // gr_psf[i]->SetMarkerColor(4);
      // 	  // gr_psf[i]->Draw("ap");
      // 	  //plot binned spots
      // 	  spot_pads[i]->SetLogz();
      // 	  //spotH_cm_clone2[i]->GetXaxis()->SetLimits(-0.,1);
      // 	  //spotH_cm_clone2[i]->GetYaxis()->SetRangeUser(-1,1);
      // 	  spotH_cm_clone2[i]->GetZaxis()->SetRangeUser(1,100000);
      // 	  spotH_cm_clone2[i]->GetXaxis()->SetLabelSize(0.035);
      // 	  spotH_cm_clone2[i]->GetYaxis()->SetLabelSize(0.035);
      // 	  spotH_cm_clone2[i]->GetYaxis()->SetTitleSize(0.03);
      // 	  spotH_cm_clone2[i]->GetXaxis()->SetTitleSize(0.03);
      // 	  spotH_cm_clone2[i]->GetXaxis()->SetTitleOffset(1.3);
      // 	  spotH_cm_clone2[i]->GetYaxis()->SetTitleOffset(1.6);
      // 	  spotH_cm_clone2[i]->SetLineColor(1);
      // 	  spotH_cm_clone2[i]->SetMarkerStyle(6);
      // 	  spotH_cm_clone2[i]->SetMarkerColor(4);
      // 	  spotH_cm_clone2[i]->Draw("colz");
      // 	}
      
      // spot_canvas33->Modified();
      // spot_canvas33->Update();
      

      //additional plots
      double xD[kN], yD[kN];
      TEllipse* psfEllipse[kN];
      TBox* pixelBox[kN];
      double psf_encircled_cm[kN];
      double psf_encircled_photons[kN];
      for(int i=0;i!=kN;++i)
	{
	  psf_encircled_cm[i] = psf_encircled_all[0][i];
	  psf_encircled_photons[i] = psf_encircled_all[3][i];
	}

      create_psf_ellipse(psfEllipse, spotH_cm, kN, psf_encircled_cm);//psf_radii
      create_pixel_box(pixelBox, spotH_cm, kN, 3.0);
      //plot_spot(spotH_cm, pixelBox, psfEllipse, kN, psf_encircled_cm, method, ofile_name, xD, yD, psf_encircled_photons, psf_encircled_cm, effectiveArea, angularResSagittal, angularResTangential, gate.fRpMax, kDegStep, fDegStep);//psf_radii

      
      
      //spot progression plot
			TH2D* cameraFrame = new TH2D("cameraFrame","Spot progression on-axis (0 degrees) to off-axis (4.5 degrees);X (cm);Y (cm)",1000,-1.,1.04*gate.fRf,1000,-1.,1.); //full camera
			TArc *cameraEdge = new TArc(0.0,0.0,gate.fRf,0.0, 6.0);
      TCanvas* spot_prog = new TCanvas("spot_prog", "spot_prog", 1500, 300);
      gPad->SetLogz();
			//spot_prog->SetMargin(0.15, 0.85, 0.15, 0.85);
			spot_prog->SetRightMargin(0.15);
			spot_progression_cm->SetStats(0);
			spot_progression_cm->GetZaxis()->SetTitle("number of photons");
			spot_progression_cm->SetTitleSize(0.08);
			spot_progression_cm->GetXaxis()->SetTitleOffset(0.7);
			spot_progression_cm->GetYaxis()->SetTitleOffset(0.7);
			spot_progression_cm->GetZaxis()->SetTitleOffset(0.4);
			spot_progression_cm->GetXaxis()->SetTitleSize(0.07);
			spot_progression_cm->GetYaxis()->SetTitleSize(0.07);
			spot_progression_cm->GetZaxis()->SetTitleSize(0.07);
			spot_progression_cm->GetXaxis()->SetLabelSize(0.05);
			spot_progression_cm->GetYaxis()->SetLabelSize(0.05);
			spot_progression_cm->GetZaxis()->SetLabelSize(0.05);
			//spot_progression_cm->GetZaxis()->SetLabelOffset(0);
			
			spot_progression_cm->Draw("colz");
			//spot_progression_cm->SetStats(kFALSE);
			//spot_progression_cm->Draw("colz");
			cameraEdge->SetLineColor(3);
			cameraEdge->SetFillColor(0);
			cameraEdge->SetFillStyle(0);
			cameraEdge->SetTheta(357.0);
			cameraEdge->Draw("same");
			//spot_progression_cm->Draw("same");
			spot_prog->Update();
			string canvas_name = method + "_spot_progression_" + ofile_name + ".png";
			spot_prog->SaveAs(canvas_name.c_str());
      
    }//end of psf section



  //plot_spot_ground(spotGround_cm, kN, method, ofile_name, kDegStep, fDegStep);
  //void plot_spot_ground(TH2D* spotGround_cm[], int tN, std::string method, std::string ofile_name, double kDegStep, double fDegStep);


  
  //----- PLOT Spot ON vs OFF and Effective area --------------------------------------


  // //spot ON vs OFF plot
  //     TPaveStats* spot_stat;
  //     TPaveStats* spot_stat2;
  
  //     TH2D* spotH_cm_clone[kN];
  //     stringstream clone_num;
  //     string clone_name;
  //     for(int i=0;i<kN;i++)
  // 	{
  // 	  clone_num << i+1;
  // 	  clone_name = "spotH_cm_clone_" + clone_num.str();
  // 	  spotH_cm_clone[i] = (TH2D*)spotH_cm[i]->Clone(clone_name.c_str());
  // 	  clone_num.str("");
  // 	}

  //     TEllipse* fovEllipse = new TEllipse(0.0,0.0,0.4,2.0,-0.2);// -17.6977,0.0,17.6977, 0.0, 2.0,-0.2);//,180);
  //     fovEllipse->SetLineStyle(2);
  //     fovEllipse->SetLineColor(3);
  //     fovEllipse->SetLineWidth(2);
  //     fovEllipse->SetFillStyle(0);
      
  //     TCanvas* spot_on_v_off = new TCanvas("spot_on_v_off", "spot_on_v_off", 1200, 600);
  //     gStyle->SetOptStat(0);
  //     gStyle->SetTitleFontSize(0.04);
  //     spot_on_v_off->Divide(2,1);
  //     TPad *pads[2];

  //     for(int i=0;i!=2;++i)
  // 	{
  // 	  pads[i] = (TPad*)spot_on_v_off->GetPad(i+1);
  // 	}
      
  //     for(int i=0;i!=2;++i)
  // 	{
  // 	  pads[i]->SetRightMargin(0.15);
  // 	  pads[i]->Draw();
  // 	}
  //     spot_on_v_off->Update();
      
  //     spot_on_v_off->cd(1);
  //     pads[0]->SetLogz();
  //     //spot_stat = (TPaveStats*)spotH_cm_clone[0]->FindObject("stats_1");
  //     //spot_stat->SetOptStat(0);
  //     spotH_cm_clone[0]->SetTitle("#theta_{field angle} = 0 (degrees)");
  //     spotH_cm_clone[0]->SetAxisRange(-0.5,0.5,"X");
  //     spotH_cm_clone[0]->SetAxisRange(-0.5,0.5,"Y");
  //     spotH_cm_clone[0]->GetZaxis()->SetRangeUser(1,100);//spotH_cm[0]->GetMaximum());
  //     spotH_cm_clone[0]->SetLabelFont(42,"X;Y;Z");
  //     spotH_cm_clone[0]->SetLabelSize(0.035,"X;Y;Z");
  //     spotH_cm_clone[0]->SetTitleSize(0.04,"X;Y;Z");
  //     spotH_cm_clone[0]->SetTitleFont(42,"X;Y;Z");
  //     spotH_cm_clone[0]->SetTitleOffset(1.2,"Y");
  //     spotH_cm_clone[0]->Draw("colz");
  //     spot_on_v_off->cd(2);

  //     pads[1]->SetLogz();
  //     //spot_stat2 = (TPaveStats*)spotH_cm_clone[8]->FindObject("stats_9");
  //     //spot_stat2->SetOptStat(0);
  //     spotH_cm_clone[kN-2]->SetTitle("#theta_{field angle} = 4.4 (degrees)");
  //     spotH_cm_clone[kN-2]->SetAxisRange(-0.5,0.5,"X");
  //     spotH_cm_clone[kN-2]->SetAxisRange(-0.5,0.5,"Y");
  //     spotH_cm_clone[kN-2]->GetZaxis()->SetRangeUser(1,100);//spotH_cm[kN-1]->GetMaximum());
  //     spotH_cm_clone[kN-2]->SetLabelFont(42,"X;Y;Z");
  //     spotH_cm_clone[kN-2]->SetLabelSize(0.035,"X;Y;Z");
  //     spotH_cm_clone[kN-2]->SetTitleSize(0.04,"X;Y;Z");
  //     spotH_cm_clone[kN-2]->SetTitleFont(42,"X;Y;Z");
  //     spotH_cm_clone[kN-2]->SetTitleOffset(1.2,"Y");
  //     spotH_cm_clone[kN-2]->Draw("colz");
  //     fovEllipse->Draw("only same");
  //     string output_name = method + "/spot_ON_versus_OFF_" + ofile_name + ".eps";
  //     spot_on_v_off->SaveAs(output_name.c_str());
  //     output_name = method + "/spot_ON_versus_OFF_" + ofile_name + ".png";
  //     spot_on_v_off->SaveAs(output_name.c_str());

  //     TCanvas* effective_area_plot = new TCanvas("effective_area_plot", "effective_area_plot", 900, 700);
  //     gStyle->SetTitleFontSize(0.04);
  //     graEffArea->SetMarkerStyle(20);
  //     graEffArea->SetMarkerColor(4);
  //     graEffArea->SetMarkerSize(1);
  //     graEffArea->GetXaxis()->SetLabelFont(42);
  //     graEffArea->GetYaxis()->SetLabelFont(42);
  //     graEffArea->GetXaxis()->SetLabelSize(0.028);
  //     graEffArea->GetYaxis()->SetLabelSize(0.028);
  //     graEffArea->GetXaxis()->SetTitleSize(0.03);
  //     graEffArea->GetYaxis()->SetTitleSize(0.03);
  //     graEffArea->GetXaxis()->SetTitleFont(42);
  //     graEffArea->GetYaxis()->SetTitleFont(42);
  //     graEffArea->SetMaximum(15.0);
  //     graEffArea->SetMinimum(1.0);
  //     graEffArea->Draw("APL");
  //     output_name = method + "/effective_area_" + ofile_name + ".eps";
  //     effective_area_plot->SaveAs(output_name.c_str());
  //     output_name = method + "/effective_area_" + ofile_name + ".png";
  //     effective_area_plot->SaveAs(output_name.c_str());


  //------- Additional plotting options -----------------------------------------------------------------
  

      // TCanvas* test2 = new TCanvas("test2","test2",800,600);
      // graTime->SetMarkerColor(4);
      // graTime->Draw("APL");
      
      //plot_focal_plane(x_array_all, y_array_all, z_array_all);
      //plot_rays(gate,poly_array);
      //plot_time_xFP(hTime_xFP, kN);
      //plot_time_yFP(hTime_yFP, kN);
      //plot_time_xyS(hTime_xyS, kN, method, ofile_name);
      //plot_time_xyP(hTime_xyP, kN, method, ofile_name);
      //plot_time_xyFP_S(hTime_xyFP_S, kN, method, ofile_name);
      //plot_time_xyFP_Total(hTime_xyFP_Total, kN, method, ofile_name);
      //plot_time_xyFP_Delay(hTime_xyFP_Delay, kN, method, ofile_name);
      //plot_mean_time_delay(kN, fov_angle, delay_mean_fov, delay_meanErr_fov, method, ofile_name);
      //plot_onAxis(hTime_xyFP_Total, method, ofile_name);
      //plot_onAxis(hist_xyP, method, ofile_name);
      //plot_onAxis(hist_xyS, method, ofile_name);
      //}
  
      cout << "inputs pre: "  << kDegStep << " " << fDegStep << endl;
      switch(set_focal_plane)
      	{
      	case 0: //camera
      	  //a higher quality on-off plot is done above
      	  break;
      	case 1: //primary
      	  plot_on_v_off(hist_xyP, method, ofile_name, 0, 9, 0, kDegStep, fDegStep); //plot ON vs OFF
      	  plot_time_xyP(hTime_xyP, kN, method, ofile_name); //plot all field (theta) angles
      	  //plot_on_v_off(TProfile2D* hist_xy[], std::string method, std::string ofile_name, int on_axis, int off_axis, int mirror_index)
      	  break;
      	case 2: //secondary
      	  plot_on_v_off(hist_xyS, method, ofile_name, 0, 9, 1, kDegStep, fDegStep);//plot ON vs OFF
      	  //plot_time_xyS(hTime_xyS, kN, method, ofile_name); //plot all field (theta) angles
      	  break;
      	}
      
      //plot_effArea_angRes(graEffArea, graSigma,method, ofile_name, leg);

      //append effective data to file for comparison

      ofstream effArea_outfile;
  effArea_outfile.open("square/effArea_data_file.txt", ios::app);
  //output effective area values to file for checking
  for(int i=0; i!=kN; i++)
    {
	  effArea_outfile.width(5);
	  effArea_outfile.precision(1);
	  effArea_outfile << fov_angle[i] << " ";
	  effArea_outfile.width(3);
	  effArea_outfile.precision(1);
	  effArea_outfile << "0.0" << " "; //phi angle
	  effArea_outfile.width(10);
	  effArea_outfile << fixed;
	  effArea_outfile.precision(5);
	  effArea_outfile << effectiveArea[i] << " ";
	  effArea_outfile.width(10);
	  effArea_outfile << effectiveArea[i] << " ";
	  //specify whether with or wihout shadowing components
	  //effArea_outfile << "noComponents" << " ";
	  effArea_outfile << "withComponents_v2" << " "; //withComponents or withComponents_v2
	  //specify whether spot analysis or shadowing analysis
	  effArea_outfile << "spot" << endl;
	  //effArea_outfile << "shadowing" << endl;
    }
  
  effArea_outfile.close();





      
      //---------------------------------------------------------------------------------------
}

void SST_GATE_PERFORMANCE::check_effective_area(double N_focused[], int tN, string method_in, double area_in[], double field_angles[], double total_photons, TGraph* effAreaAltMethod)
{

  double N_simulated[tN];
  double effective_area[tN];

  for(int i=0;i!=tN;++i)
    {
      if(method_in=="square") N_simulated[i] = total_photons;//400. * 400.;  //photon number
      else if(method_in=="cone") N_simulated[i] = 1.52053e5; //photon number
    }

  for(int i=0;i!=tN;++i)
    {
      effective_area[i] = ( N_focused[i] / N_simulated[i] ) * area_in[i];
    }

  TGraph* gr_eEffA = new TGraph(tN, field_angles, effective_area);
  TCanvas* cEffA = new TCanvas("cEffA", "cEffA", 800,600);
  cEffA->SetGrid(1,1);

  TMultiGraph* mg1 = new TMultiGraph();
  // gr_eEffA->SetMaximum(12.0);
  // gr_eEffA->SetMinimum(0.0);
  gr_eEffA->SetMarkerStyle(24);  
  gr_eEffA->SetMarkerColor(4);
  gr_eEffA->SetLineColor(4);
  // gr_eEffA->Draw("ALP");
  mg1->Add(gr_eEffA);
  effAreaAltMethod->SetMarkerStyle(26);
  effAreaAltMethod->SetMarkerColor(1);
  effAreaAltMethod->SetLineColor(1);
  //effAreaAltMethod->Draw("sameLP");
  mg1->Add(effAreaAltMethod);
  mg1->SetMaximum(9.5); //withComponents = 8.2 noComponents = 9.5
  mg1->SetMinimum(5.5);  //withComponents = 5.5 noComponents = 7.5
  
  mg1->SetTitle("Effective area (calculation cross-check); #theta_{field angle} (deg.); A_{effective} (m^{2})");
  mg1->Draw("ALP");
  
  //TLegend* leg1 = new TLegend(0.15, 0.65, 0.5, 0.85);//noComponents
  TLegend* leg1 = new TLegend(0.136935, 0.198953, 0.487437, 0.399651); //withComponents
  leg1->SetFillStyle(1001);
  leg1->SetFillColor(0);  
  leg1->AddEntry(gr_eEffA, "A_{eff} = #frac{N_{focused}}{N_{simulated}} #times A_{scattering}" , "lp");
  leg1->AddEntry(effAreaAltMethod, "A_{eff} = dA_{single photon} #times N_{simulated}" , "lp");
  leg1->Draw();
  
  cEffA->Update();
  cEffA->Modified();

  cEffA->SaveAs("effective_area_cross_check.pdf");
  cEffA->SaveAs("effective_area_cross_check.png");

}

double calculateRMS(vector<double> &vTimes)
{

  const int kN = vTimes.size();
  double sumMeans;
  
  for(auto itTime = vTimes.begin(); itTime != vTimes.end(); ++itTime)
    {
      double getVal = (*itTime);
      sumMeans += (getVal*getVal);
    }

  double mean;
  if(kN!=0) mean = (sumMeans/kN);
  else mean=0.0;

  double RMS = sqrt(mean);

  return RMS;
}

double calculateMEAN(vector<double> &vTimes)
{

  const int kN = vTimes.size();
  double sumMeans;
  
  for(auto itTime = vTimes.begin(); itTime != vTimes.end(); ++itTime)
    {
      sumMeans += (*itTime);
    }

  double mean;
  if(kN!=0) mean = (sumMeans/kN);
  else mean=0.0;

  return mean;
}

#endif


// //---------  CACULATE PSF SECTION ---------//

// std::pair<double,double> SST_GATE_PERFORMANCE::get_encircled_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn)
// {

//   std::vector<double>::iterator x_iter;
//   std::vector<double>::iterator y_iter;
//   double totalPhotonsX = x_arrayIn.size(); 
//   double totalPhotonsY = y_arrayIn.size();

//   //cout << "NEW: " << totalPhotonsX << " " << totalPhotonsY << endl;

//   double mX = getMean(x_arrayIn);
//   double mY = getMean(y_arrayIn);
//   //cout << "NEW mX: " << mX << " mY: " << mY << endl;
//   double radius = mX+stepIn;
//   double ratio = 0.0;
//   int photonsInside;

//   while(ratio<80)
//     {
//       double photon_radius=0.0;  
//       photonsInside=0;      

//       y_iter=y_arrayIn.begin();
//       for(x_iter=x_arrayIn.begin();x_iter!=x_arrayIn.end();++x_iter)
// 	{
// 	  double x = (*x_iter);
// 	  double y = (*y_iter);
// 	  photon_radius = sqrt(x*x + y*y);
// 	  if(photon_radius<=radius)
// 	    {
// 	      //cout << " x: " << x << " y: " << y << " photon_radius: " << photon_radius << endl;
// 	      photonsInside++;
// 	    }
	  
// 	  ++y_iter;
// 	}
      
//       ratio = static_cast<int>(photonsInside/totalPhotonsX*100.0);
//       //cout << " photonsInside: " << photonsInside << " photonsInside/totalPhotons :" << ratio << endl;
//       radius+=stepIn;
//     }

//   cout << " final radius = " << radius << " photons: " << photonsInside << endl;
//   cout << endl;

//   return std::make_pair(radius,photonsInside);
// }

// std::pair<double,double> SST_GATE_PERFORMANCE::get_ensquared_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn, double sqHalfLength)
// {

//   std::vector<double>::iterator x_iter;
//   std::vector<double>::iterator y_iter;
//   double totalPhotonsX = x_arrayIn.size(); 
//   double totalPhotonsY = y_arrayIn.size();

//   //cout << "NEW: " << totalPhotonsX << " " << totalPhotonsY << endl;

//   double mX = getMean(x_arrayIn);
//   double mY = getMean(y_arrayIn);
//   //cout << "NEW mX: " << mX << " mY: " << mY << endl;

//   double squareXminus = mX-stepIn;
//   double squareXplus = mX+stepIn;
//   double squareYminus = mY-stepIn;
//   double squareYplus = mY+stepIn;

//   //double radius = mX+stepIn;
//   double ratio = 0.0;
//   int photonsInside;

//   while(ratio<80)
//     {
//       double photon_radius=0.0;  
//       photonsInside=0;      

//       y_iter=y_arrayIn.begin();
//       for(x_iter=x_arrayIn.begin();x_iter!=x_arrayIn.end();++x_iter)
// 	{
// 	  double x = (*x_iter);
// 	  double y = (*y_iter);
// 	  if(x>=(mX-sqHalfLength) and x<=(mX+sqHalfLength) and y>=(mY-sqHalfLength) and y<=(mY+sqHalfLength))
// 	    {
// 	      photonsInside++;
// 	    }
// 	  // photon_radius = sqrt(x*x + y*y);
// 	  // if(photon_radius<=radius)
// 	  //   {
// 	  //     //cout << " x: " << x << " y: " << y << " photon_radius: " << photon_radius << endl;
// 	  //     photonsInside++;
// 	  //   }
	  
// 	  ++y_iter;
// 	}
      
//       ratio = static_cast<int>(photonsInside/totalPhotonsX*100.0);
//       //cout << " photonsInside: " << photonsInside << " photonsInside/totalPhotons :" << ratio << endl;
//       sqHalfLength+=stepIn;
//       //radius+=stepIn;
//     }

//   cout << " final square half length = " << sqHalfLength << " photons: " << photonsInside << endl;
//   //cout << " final radius = " << radius << " photons: " << photonsInside << endl;
//   cout << endl;

//   return std::make_pair(sqHalfLength,photonsInside);
// }


// double SST_GATE_PERFORMANCE::getMean(std::vector<double> inArray)
// {

//   std::vector<double>::iterator iter;
//   double sum=0.0;
//   double mean=0.0;
//   for(iter=inArray.begin();iter!=inArray.end();++iter)
//     {
//       sum += (*iter);
//     }
//   if(inArray.size()!=0)
//     {
//       mean = sum/inArray.size();
//     }
//   else
//     {
//       cout << "Your array is empty! You cannot perform this calculation!" << endl;
//       exit(1);
//     }
//   return mean;
// }

//#endif


// //Old code to print X,Y 1 dimensional coordinate distributions.
//   double Xmin[kN],Xmax[kN],Ymin[kN],Ymax[kN];
  
//   cout<< endl;
//   for(int i=0;i<kN;i++)
//     {
//       //cout << "max: " << spotX_cm[i]->FindFirstBinAbove(0,1) << " max: " << spotX_cm[i]->FindLastBinAbove(0,1) << endl;
//       //cout << "first: " << spotX_cm[i]->GetBinLowEdge(spotX_cm[i]->FindFirstBinAbove(0,1)) << " last: " << spotX_cm[i]->GetBinLowEdge(spotX_cm[i]->FindLastBinAbove(0,1)) << endl;
//       cout<< endl;
      
//       Xmin[i] = spotX_cm[i]->GetBinLowEdge(spotX_cm[i]->FindFirstBinAbove(0,1));
//       Xmax[i] = spotX_cm[i]->GetBinLowEdge(spotX_cm[i]->FindLastBinAbove(0,1));
//       Ymin[i] = spotY_cm[i]->GetBinLowEdge(spotY_cm[i]->FindFirstBinAbove(0,1));
//       Ymax[i] = spotY_cm[i]->GetBinLowEdge(spotY_cm[i]->FindLastBinAbove(0,1)) + spotY_cm[i]->GetBinWidth(i);
//     }
  
//   TF1* fX1[kN];
  
//   for(int i=0;i<kN;i++)
//     {
// 	fX1[i] = new TF1("fX1","gaus(0)",Xmin[i],Xmax[i]);
// 	//fX1[i]->SetParameter(2,0.8);
//     }
  
//   // for(int i=0;i<kN;i++)
//   //   {
//   //     spotX_cm[i]->Fit(fX1[i]);//RB0
//   //   }
  

//   //To plot the X,Y coordinate distributions uncomment the following section.

//   TCanvas* cSpotX = new TCanvas("cSpotX","cSpotX",1000,700);
//   cSpotX->Divide(5,2);
//   for(int i=0;i<kN;i++)
//     {
//       cSpotX->cd(i+1);
//       //gPad->SetLogy();
//       spotX_cm[i]->SetAxisRange(Xmin[i],Xmax[i],"X");
//       spotX_cm[i]->Draw();
//     }

//   TCanvas* cSpotY = new TCanvas("cSpotY","cSpotY",1000,700);
//   cSpotY->Divide(5,2);
//   for(int i=0;i<kN;i++)
//     {
//       cSpotY->cd(i+1);
//       spotY_cm[i]->SetAxisRange(Ymin[i],Ymax[i],"X");
//       spotY_cm[i]->Draw();
//     }
  
//REMOVED FROM MAIN CODE 01 MARCH 2014

//void write_data(std::vector< std::vector<double> > psf_encircled_in, std::vector< std::vector<double> > psf_ensquared_in, std::vector<double> fov_angles_in, TH2D* spotH_cm[], double effectiveArea[], double angularResSagittal[], double angularResTangential[], double maxRprimary, std::string method, std::string ofile_name);


  // if(set_focal_plane==0)
  //   {  
  //     std::vector< vector<double> >::iterator xIter;
  //     std::vector< vector<double> >::iterator yIter;
  //     double psf_encircled_cm[kN], psf_encircled_deg[kN], psf_encircled_arcmin[kN];
  //     double psf_encircled_error[kN], psf_encircled_photons[kN];
  //     int countB=0;
  //     double step=0.0001;
      
  //     yIter = y_array_all.begin();
      
  //     for(xIter=x_array_all.begin();xIter!=x_array_all.end();++xIter)
  // 	{
  // 	  //cout << "d2size: " << x_array_all.size() << " " << y_array_all.size() << " " << (*xIter).size() << " " << (*yIter).size() << endl;
  // 	  std::pair<double,double> psf_vals = get_encircled_PSF((*xIter),(*yIter),step); //[cm] [first=psf,second=photons]
  // 	  //std::pair<double,double> psf_vals = get_ensquared_PSF((*xIter),(*yIter),step,0.01);//start with a 1mm pixel i.e. 0.1cm
  // 	  psf_encircled_cm[countB] = psf_vals.first;                  //[cm]
  // 	  psf_encircled_deg[countB] = psf_vals.first/gate.deg2dist;        //[deg]
  // 	  psf_encircled_arcmin[countB] = psf_vals.first/gate.deg2dist*60.0;//[arcmin]
  // 	  psf_encircled_photons[countB] = psf_vals.second;
  // 	  psf_encircled_error[countB] = step/gate.deg2dist;                //[deg] for [arcmin] * 60.0
  // 	  countB++;
  // 	  ++yIter;
  // 	}
      
  //     //setup output file
  //     ofstream photonsFile;
  //     string output_name = method + "/data_" + ofile_name + ".txt";
  //     photonsFile.open(output_name.c_str(),ios::app);

  //     //print file headers
  //     photonsFile << "#FOV(deg): Tot.Phot.: PSF.Phot.: PSF(cm): PSF(deg): PSF(arcmin): EffArea(m2): AngResS(arcmin): AngResT(arcmin): M1radius(cm): Method:" << endl;

  //     for(int i=0;i<10;i++)
  // 	{
  // 	  //output values to file
  // 	  photonsFile << fov_angle[i] << "\t" << spotH_cm[i]->GetEntries() << "\t" << psf_encircled_photons[i] << "\t" << psf_encircled_cm[i] << "\t" << psf_encircled_deg[i] << "\t" << psf_encircled_arcmin[i] << "\t" << effectiveArea[i] << "\t" << angularResSagittal[i] << "\t" << angularResTangential[i] << "\t" << gate.fRpMax << "\t" << method << endl;

  // 	  //output values to screen
  // 	  cout << "\tpsf: \t\tfov: \tphotons:" << endl;
  // 	  cout << "\t" << psf_encircled_cm[i] << "\t" << fov_angle[i] << "\t" << psf_encircled_photons[i] << "\t [cm, deg, number]" << endl;
  // 	  cout << "\t" << psf_encircled_deg[i] << "\t" << fov_angle[i] << "\t" << psf_encircled_photons[i] << "\t [deg, deg, number]" <<endl;
  // 	  cout << "\t" << psf_encircled_arcmin[i] << "\t" << fov_angle[i] << "\t" << psf_encircled_photons[i] << "\t [arcmin, deg, number]" <<endl;
  // 	  cout << endl;
  // 	}
  //     photonsFile.close();
  
      
  //     //TGraph* psf_graph = new TGraph(kN,fov_angle,psf);
  //     TGraphErrors* psf_graph = new TGraphErrors(kN,fov_angle,psf_encircled_deg,0,psf_encircled_error);//psf_error
  //     TCanvas* plotPSF = new TCanvas("plotPSF", "plotPSF", 900, 700);
  //     gStyle->SetTitleFontSize(0.04);
  //     psf_graph->SetTitle("Evolution of the 80% containment radius with field angle; field angle (degrees);PSF_{80% containment radius} (degrees)");
  //     psf_graph->GetXaxis()->SetRangeUser(0.0,4.9);  
  //     psf_graph->GetYaxis()->SetRangeUser(0.02,0.07);
  //     psf_graph->GetXaxis()->SetTitleSize(0.03);
  //     psf_graph->GetYaxis()->SetTitleSize(0.03);
  //     psf_graph->GetYaxis()->SetTitleOffset(1.5);
  //     psf_graph->GetXaxis()->SetLabelSize(0.028);
  //     psf_graph->GetYaxis()->SetLabelSize(0.028);
  //     psf_graph->SetMarkerStyle(20);
  //     psf_graph->SetMarkerColor(4);
  //     psf_graph->Draw("ALP");
  //     TGaxis* arcminAxis = new TGaxis(4.9,0.02,4.9,0.07,1.2,4.2,510,"+L");
  //     arcminAxis->SetName("arcminAxis");
  //     arcminAxis->SetTitle("PSF_{80% containment radius} (arcmin)");
  //     arcminAxis->SetTitleOffset(0.8);
  //     arcminAxis->SetTitleSize(0.03);
  //     arcminAxis->SetLabelSize(0.028);
  //     arcminAxis->SetLabelFont(42);
  //     arcminAxis->SetTitleFont(42);
  //     arcminAxis->Draw();
      
  //     output_name = method + "/psf_NEW_deg_full_" + ofile_name + ".png";    //  string output_name;
  //     plotPSF->SaveAs(output_name.c_str());
  //     output_name = method + "/psf_NEW_deg_full_" + ofile_name + ".eps";
  //     plotPSF->SaveAs(output_name.c_str());
      
  //     //  plotPSF->Divide(5, 2);
  //     // for(int i=0;i<kN;i++)
  //     // {
  //     //  plotPSF->cd(i+1);
  //     //  psf_graphs[i]->Draw("A*");
  //     // }
