#ifndef SST_GATE_H
#define SST_GATE_H

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

#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
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
#include "TSystem.h"

#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoCompositeShape.h"
#include "TGeoPara.h"
#include "TGeoPgon.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoTrd2.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoTorus.h"
#include "TGeoTrd1.h"
#include "TGeoXtru.h"

#include "AOpticsManager.h"
#include "AOpticalComponent.h"
#include "AGeoAsphericDisk.h"
#include "AMirror.h"
#include "AFocalSurface.h"
#include "AObscuration.h"
#include "ARay.h"
#include "ARayShooter.h"
#include "ALens.h"
#include "ABorderSurfaceCondition.h"
#include "AGlassCatalog.h"

#include "segmented_mirrors.hpp"
#include "sector_segmented_mirrors.hpp"
#include "sst_gate_globals.hpp"

//using namespace std;

class SST_GATE: public clonable{
private:
  AOpticsManager* fManager;
  AGeoAsphericDisk* fPrimaryV;
  AGeoAsphericDisk* fSecondaryV;
  
public:
  SST_GATE();
  virtual ~SST_GATE();
  virtual SST_GATE* clone() const { return new SST_GATE( *this ); }
  
  void StartBuildingTelescope();
  void AddIdealFocalPlane();
  void AddIdealGroundFocalPlane();
  void AddCameraBody();
  void AddCameraLid(double deltaPos, double depth);
  void AddTelescopeFrame();
  void AddTelescopeFrame_version2();
  void Add_CHEC_MAPMT_FocalPlane(double deltaZ);
  void AddMAPMTFocalPlane(double deltaZ);
  void AddPrimaryMask();
  void AddPrimaryMirror(const char* name, SegmentedMirror* mirror);
  void AddSecondaryMirror(const char* name, SegmentedMirror* mirror);
  
  void AddSecondaryObscuration();
  void CloseGeometry() {fManager->CloseGeometry();}
  AOpticsManager* GetManager() const {return fManager;}
  void InitParameters();
  void setNewPrimaryMaxRadius(double new_radius);
  void AddPrimaryMirrorAsFP(const char* name, SegmentedMirror* mirror);
  void AddSecondaryMirrorAsFP(const char* name, SegmentedMirror* mirror);
  
  // Define the telescope parameters
  double fF;       // Focal length
  double fAlpha;   // \alpha
  double fQ;       // q
  double fRpMax;   // Primary radius max
  double fRpMin;   // Primary radius min
  double fRsMax;   // Secondary radius max
  double fRsMin;   // Secondary radius min
  double fKappa1;  // Focal plane sag constant
  double fKappa2;  // Focal plane sag constant
  double fZs;      // position of secondary
  double fZf;      // position of focal plane
  double fPlate;   // plate scale in mm/arcsec
  double deg2dist; //convert plate scale from mm/arcmin to cm.
  
  int fNp;         // Number of coefficients for the primary
  int fNs;         // Number of coefficients for the secondary
  double* fP;      // Polynomial coefficients (p0, p1 ...)
  double* fS;      // Polynomial coefficients (s0, s1 ...)
    
  // Focal plane
  double fRf;      // radius of focal plane
};

SST_GATE::SST_GATE()
{
  fManager = new AOpticsManager("manager", "the optics manager for the SST-GATE prototype");
  fManager->SetVisLevel(5);
  fManager->SetNsegments(50);
  
  // Make dummy material
  TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
  mat->SetTransparency(70); // needed in OpenGL view
  new TGeoMedium("med", 1, mat);
  
  // Make the world
  TGeoBBox* worldbox = new TGeoBBox("worldbox", 30*m, 30*m, 30*m);
  AOpticalComponent* world = new AOpticalComponent("world", worldbox);
  fManager->SetTopVolume(world);
}

SST_GATE::~SST_GATE()
{
  // Note that all the geometry components are automatically deleted by ROOT
  delete fManager;
  delete [] fP;
  delete [] fS;
}

void SST_GATE::StartBuildingTelescope()
{
  InitParameters();
  const double kZp = 0.*m;
  const double kZs = fF/fQ;  // fF/fQ = separation distance between primary and secondary
  
  // Make the ideal volume of the primary mirror
  fPrimaryV = new AGeoAsphericDisk("primaryV", kZp + fP[0] - 1*um, 0, kZp + fP[0] , 0, fRpMax, fRpMin);
  fPrimaryV->SetPolynomials(fNp - 1, &fP[1], fNp - 1, &fP[1]);
  std::cout << "CHECK!" << std::endl;
  // Make the ideal volume of the secondary mirror
  fSecondaryV = new AGeoAsphericDisk("secondaryV",kZs + fS[0], 0, kZs + fS[0] + 1.*cm, 0, fRsMax, fRsMin);// + 1*um
  fSecondaryV->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);
}

void SST_GATE::AddPrimaryMirrorAsFP(const char* name, SegmentedMirror* mirror)
{
  AFocalSurface* primaryMirrorAsFocalP = new AFocalSurface("primaryMirrorAsFocalP", fPrimaryV);
  //AMirror* mir=mirror->BuildMirror(name, fPrimaryV, kTRUE);
  TGeoCombiTrans* combi=mirror->BuildMirrorCombiTrans(fPrimaryV, kTRUE);
  
  ABorderSurfaceCondition * condition = new ABorderSurfaceCondition((AOpticalComponent*)fManager->GetTopVolume(), primaryMirrorAsFocalP);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  //  fManager->GetTopVolume()->AddNode(mir, 1, combi);
  fManager->GetTopVolume()->AddNode(primaryMirrorAsFocalP, 1);
}


void SST_GATE::AddPrimaryMirror(const char* name, SegmentedMirror* mirror)
{
  AMirror* mir=mirror->BuildMirror(name, fPrimaryV, kTRUE);
  TGeoCombiTrans* combi=mirror->BuildMirrorCombiTrans(fPrimaryV, kTRUE);
  
  ABorderSurfaceCondition * condition = new ABorderSurfaceCondition((AOpticalComponent*)fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  fManager->GetTopVolume()->AddNode(mir, 1, combi);
}

void SST_GATE::AddPrimaryMask()
{
  const double kZp = 0.*m;
  
  AGeoAsphericDisk* disk=new AGeoAsphericDisk("primaryObsV", kZp + fP[0] + 1*um, 0, kZp + fP[0] , 0, fRpMax, fRpMin);
  disk->SetPolynomials(fNp - 1, &fP[1], fNp - 1, &fP[1]);
  
  //TGeoMedium* med=fManager->GetMedium("med");
  AObscuration* primaryObs=new AObscuration("primaryObs", disk);//, med);
  //fManager->GetTopVolume()->AddNode(primaryObs, 1);

  TGeoTrd1 * primarySupportFrame = new TGeoTrd1("primarySupportFrame", 0.*mm, 620.*mm * TMath::Tan(30*TMath::DegToRad()), 1.*um, 310.*mm);
  AObscuration* primaryObs2 = new AObscuration("primaryObs2", primarySupportFrame);
  for(int i=0;i<6;i++)
    {
      switch(i)
	{
	case 0:
	  fManager->GetTopVolume()->AddNodeOverlap(primaryObs2, i+1, new TGeoCombiTrans(TGeoTranslation(2.*m-310.*mm,0.,17.*cm), TGeoRotation("", 90, 80, 0)));
	  break;
	case 1:
	  fManager->GetTopVolume()->AddNodeOverlap(primaryObs2, i+1, new TGeoCombiTrans(TGeoTranslation(-2.*m+310.*mm,0.,17.*cm), TGeoRotation("", 90, -80, 0)));
	  break;
	case 2:
	  fManager->GetTopVolume()->AddNodeOverlap(primaryObs2, i+1, new TGeoCombiTrans(TGeoTranslation((2.*m-310.*mm)*TMath::Cos(60*TMath::DegToRad()),(2.*m-310.*mm)*TMath::Sin(60*TMath::DegToRad()),17.*cm), TGeoRotation("", 60+90, 80, 0)));
	  break;
	case 3:
	  fManager->GetTopVolume()->AddNodeOverlap(primaryObs2, i+1, new TGeoCombiTrans(TGeoTranslation((2.*m-310.*mm)*TMath::Cos(120*TMath::DegToRad()),(2.*m-310.*mm)*TMath::Sin(120*TMath::DegToRad()),17.*cm), TGeoRotation("", 120+90, 80, 0)));
	  break;
	case 4:
	  fManager->GetTopVolume()->AddNodeOverlap(primaryObs2, i+1, new TGeoCombiTrans(TGeoTranslation((2.*m-310.*mm)*TMath::Cos(240*TMath::DegToRad()),(2.*m-310.*mm)*TMath::Sin(240*TMath::DegToRad()),17.*cm), TGeoRotation("", 240+90, 80, 0)));
	  break;
	case 5:
	  fManager->GetTopVolume()->AddNodeOverlap(primaryObs2, i+1, new TGeoCombiTrans(TGeoTranslation((2.*m-310.*mm)*TMath::Cos(300*TMath::DegToRad()),(2.*m-310.*mm)*TMath::Sin(300*TMath::DegToRad()),17.*cm), TGeoRotation("", 300+90, 80, 0)));
	  break;
	}
    }
}

void SST_GATE::AddSecondaryMirrorAsFP(const char* name, SegmentedMirror* mirror)
{
  AFocalSurface* secondaryMirrorAsFocalP = new AFocalSurface("secondaryMirrorAsFocalP", fSecondaryV);
  //AMirror* mir=mirror->BuildMirror(name, fSecondaryV, kFALSE);
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fSecondaryV, kFALSE);
  
  ABorderSurfaceCondition * condition=new ABorderSurfaceCondition((AOpticalComponent*)fManager->GetTopVolume(), secondaryMirrorAsFocalP);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  
  //fManager->GetTopVolume()->AddNode(mir, 1, combi);
  fManager->GetTopVolume()->AddNode(secondaryMirrorAsFocalP, 1);
}


void SST_GATE::AddSecondaryMirror(const char* name, SegmentedMirror* mirror)
{
  AMirror* mir=mirror->BuildMirror(name, fSecondaryV, kFALSE);
  TGeoCombiTrans* combi = mirror->BuildMirrorCombiTrans(fSecondaryV, kFALSE);
  
  ABorderSurfaceCondition * condition=new ABorderSurfaceCondition((AOpticalComponent*)fManager->GetTopVolume(), mir);
  condition->SetGaussianRoughness(mirror->GetRoughness()*TMath::DegToRad());
  
  fManager->GetTopVolume()->AddNode(mir, 1, combi);
}

void SST_GATE::AddSecondaryObscuration()
{
  const double kZs = fF/fQ; // fF/fQ = separation distance between primary and secondary
  
  AGeoAsphericDisk* disk=new AGeoAsphericDisk("secondaryObsV", kZs + fS[0], 0, kZs + fS[0] + 1.*cm, 0, fRsMax, 0);
  disk->SetPolynomials(fNs - 1, &fS[1], fNs - 1, &fS[1]);

  //TGeoMedium* med=fManager->GetMedium("med");
  AObscuration* secondaryObs=new AObscuration("secondaryObs", disk);
  secondaryObs->SetLineColor(24);
  //fManager->GetTopVolume()->AddNode(secondaryObs, 1);
  fManager->GetTopVolume()->AddNode(secondaryObs, 1, new TGeoCombiTrans(TGeoTranslation(0.,0.,1.*cm), TGeoRotation("", 0, 0, 0)));

  //TGeoBBox * primarySupportFrame = new TGeoBBox("primarySupportFrame", 25.0*mm, 125.0*mm, v4_PS.Mag());
  // TGeoPgon * secondaryPolySupportFrame = new TGeoPgon("secondaryPolySupportFrame", 360., 360.,8,2);
  // AObscuration* secondaryObs2 = new AObscuration("secondaryObs2", secondaryPolySupportFrame);

  // fManager->GetTopVolume()->AddNodeOverlap(secondaryObs2, 1, new TGeoCombiTrans(TGeoTranslation(0.,0.,kZs), TGeoRotation("", 0, 0, 0)));

}

void SST_GATE::setNewPrimaryMaxRadius(double new_radius)
{
fRpMax  = new_radius*m;
}

void SST_GATE::InitParameters()
{

  //initialise sst_gate parameters
  fF      = 2.283*m;       //404*m;//5km //2.283*m;     // Focal length @infinity   2.28352@10km 2.28404@5km
  fAlpha  = 0.776303109;   // (1 - fAlpha)*fF = separation distance between secondary & focal plane 
  fQ      = 0.641094044;   // fF/fQ = separation distance between primary and secondary
  //fRpMax  = newM1Radius*m;         // Primary radius max D=4.0m
  //fRpMax  = 2.0*m;         // Primary radius max D=4.0m
  fRpMin  = 0.65*m;        // Primary radius min D=1.3m //0.65*m;
  fRsMax  = 1.*m;          // Secondary radius max. D=2.002m
  fRsMin  = 0.0*m;         // Secondary radius min
  fKappa1 = 0.0;           // Focal plane sag constant
  fKappa2 = 16.517*mm;     // Focal plane sag constant
  fZs = 3.5611*m;          // Secondary position //Frederic=3561.075 //original=3.5611*m;
  fZf = (fZs-0.5107)*m;    // Focal plane position //Frederic=510.694 fZf should equal 3050.381 thus fZs must be 3561.075 //original=0.5107
  fPlate = 0.67037037;     //plate scale [mm/arcmin]
  deg2dist = fPlate*mm*60.;//convert plate scale from mm/arcmin to cm.

  fRf = 18.1*cm; //36.2*cm * 0.5;       //radius of the Focal plane
  
  //initialise polynomial coefficients these are takem from Jurgen Schmoll v.17
  fNp = 9;
  // Olivier & Andreas = {5.4280255e-4, 3.3912879e-10, -1.3451359e-13, 1.2900035e-17, -6.8508142e-22, 2.0059722e-26, -3.0563336e-31, 1.8853301e-36};
  // Jurgen = {5.4280255e-5, 3.3912879e-13, -1.3451359e-18, 1.2900035e-24, -6.8508142e-31, 2.0059722e-37, -3.0563336e-44, 1.8853301e-51}; //signs flipped
  // Denis = {5.4280255e-005, 3.3912879e-013, -1.3451359e-018, 1.2900035e-024, -6.8508142e-031, 2.0059722e-037, -3.0563336e-044, 1.8853301e-051}
  // Denis separation distance: M1->M2: 3561.1mm, M2->Camera: 510.7mm, Camera Diam: 362mm, M2 Diam:2002mm m1 Diam: 4000mm

  fP = new double[fNp];
  fP[0] = 0.0;
  fP[1] = 5.4280255e-4;//4
  fP[2] = 3.3912879e-10;//10
  fP[3] = -1.3451359e-13;//13
  fP[4] = 1.2900035e-17;//17
  fP[5] = -6.8508142e-22;//22
  fP[6] = 2.0059722e-26;//26
  fP[7] = -3.0563336e-31;//31
  fP[8] = 1.8853301e-36;//36

  fNs = 9;
  //polynomial: a2, a4, a6, a8, a10, a12, a14, a16
  // Olivier & Andreas = {-2.435033e-3, -3.8794144e-8, 1.3111154e-11, -2.8830939e-15, 3.9781971e-19, -3.3371289e-23, 1.542123e-27, -2.987865e-32};
  // Jurgen = {-2.435033e-4, -3.8794144e-11, 1.3111154e-16, -2.8830939e-22, 3.9781971e-28, -3.3371289e-34, 1.542123e-40, -2.987865e-47};//signs flipped
  // Denis = {0.0002435033, 3.8794144e-011, -1.3111154e-016, 2.8830939e-022, -3.9781971e-028, 3.3371289e-034, -1.542123e-040, 2.987865e-047}
  fS = new double[fNs]; 
  fS[0] = 0.0;
  fS[1] = -2.435033e-3;//3
  fS[2] = -3.8794144e-8;//8
  fS[3] = 1.3111154e-11;//11
  fS[4] = -2.8830939e-15;//15
  fS[5] = 3.9781971e-19;//19
  fS[6] = -3.3371289e-23;//23
  fS[7] = 1.542123e-27;//27
  fS[8] = -2.987865e-32;//32

  // fSv.push_back(0.0);
  // fSv.push_back(-2.435033e-3);
  // fSv.push_back(-3.8794144e-8);
  // fSv.push_back(1.3111154e-11);
  // fSv.push_back(-2.8830939e-15);
  // fSv.push_back(3.9781971e-19);
  // fSv.push_back(-3.3371289e-23);
  // fSv.push_back(1.542123e-27);
  // fSv.push_back(-2.987865e-32);


}

//void AddPrimaryF(SST_GATE* sct, double margin = 10.0*mm)
void AddPrimaryF(SST_GATE& sct, double marginIn, bool primaryAsFP)
{
  int count = 0;
  
  // P1 mirrors
  for(int i = 0; i < 6; i++){
    double rmin = (1.3*m * 0.5) - marginIn/TMath::Cos(11.25/2.*TMath::DegToRad());
    double rmax = 4.0*m * 0.5;
    
    double phimin = 60.0*i;
    double phimax = 60.0*(i + 1);
    SectorSegmentedMirrors mirror(rmin, rmax, phimin, phimax);
    //To use marginIn to increase the inner radius uncomment the following line
    //SectorSegmentedMirrors mirror(rmin + marginIn, rmax, phimin, phimax);
    mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
    mirror.SetRotationErrors(0., 0., 0.);
    mirror.SetRoughness(0.0);
    mirror.SetMargin(marginIn);
    if(primaryAsFP==true) sct.AddPrimaryMirrorAsFP(Form("primary_mirror%d", count), &mirror);
    else sct.AddPrimaryMirror(Form("primary_mirror%d", count), &mirror);
    //sct.AddPrimaryMirrorAsFP(Form("primary_mirror%d", count), &mirror); //To add primary mirror as a focal plane uncomment this line.
    count++;
  } // i
}

//void AddSecondaryJ(SST_GATE* sct, double margin = 0.0*mm)
void AddSecondaryJ(SST_GATE& sct, double marginIn, bool setMonolithic, bool secondaryAsFP)
{
  int count = 0;
  
  if(setMonolithic)
    {
      double rmin = 0.*m;//InitParam has not been run yet so we have to declare these manually. If changing diameters remember this!
      double rmax = 1.*m;
      double phimin = 0.0;
      double phimax = 360.0;
      SectorSegmentedMirrors mirror(rmin, rmax, phimin, phimax);
      mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
      mirror.SetRotationErrors(0., 0., 0.);
      mirror.SetRoughness(0.0);
      //mirror.SetMargin(marginIn); //If monolithic, then can have no margin as no petals.
      if(secondaryAsFP==true) sct.AddSecondaryMirrorAsFP(Form("secondary_mirror%d", 1), &mirror);//To add secondary mirror as a focal plane uncomment this line.
      else sct.AddSecondaryMirror(Form("secondary_mirror%d", 1), &mirror);
    }
  else
    {
      // S1 mirrors
      for(int i=0; i<2; i++)
	{
	  double rmin = 0.*m;//15.24*cm;//0.0381*m;
	  double rmax = 2.002*m * 0.5;
	  
	  double phimin = 180.*i;
	  double phimax = 180.*(i + 1);
	  SectorSegmentedMirrors mirror(rmin, rmax, phimin, phimax);
	  //To use marginIn to increase the inner radius uncomment the following line
	  //SectorSegmentedMirrors mirror(rmin + marginIn, rmax, phimin, phimax);
	  mirror.SetPositionErrors(0*mm, 0*mm, 0*mm);
	  mirror.SetRotationErrors(0., 0., 0.);
	  mirror.SetRoughness(0.0);
	  mirror.SetMargin(marginIn);
	  if(secondaryAsFP==true) sct.AddSecondaryMirrorAsFP(Form("secondary_mirror%d", count), &mirror);//To add secondary mirror as a focal plane uncomment this line.
	  else sct.AddSecondaryMirror(Form("secondary_mirror%d", count), &mirror);
	  count++;
	} // i
    }
}

void SST_GATE::AddIdealFocalPlane()
{
  const double kZs = fF/fQ;                 // fF/fQ = separation distance between primary and secondary
  const double kZf = kZs - (1 - fAlpha)*fF; // (1 - fAlpha)*fF = separation distance between secondary & focal plane
  //const double kZf = fZf;

  AGeoAsphericDisk* idealCameraV = new AGeoAsphericDisk("idealCameraV", kZf - 1*um, 0, kZf, 0, fRf, 0);
  //double sagPar[2] = {fKappa1*TMath::Power(fF, -1),
  //                    fKappa2*TMath::Power(fF, -3)};
  //idealCameraV->SetPolynomials(2, sagPar, 2, sagPar);

  double sagPar[5] = {-5.0e-3, -1.25e-7, -6.25e-12, -3.90625e-16, -2.734375e-20}; // Yi focal plane surface
  idealCameraV->SetPolynomials(5, sagPar, 5, sagPar);

  AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV", idealCameraV->CalcF1(fRf) - 1*cm, 0, idealCameraV->CalcF1(fRf), 0, fRf, 0.);
  
  AFocalSurface* idealCamera = new AFocalSurface("idealCamera", idealCameraV);
  AObscuration* idealCameraObs = new AObscuration("idealCameraObs", idealCameraV);
  //AObscuration* idealCameraObs = new AObscuration("idealCameraObs", focalObsV);
  fManager->GetTopVolume()->AddNode(idealCamera, 1);
  fManager->GetTopVolume()->AddNode(idealCameraObs, 1, new TGeoTranslation(0, 0, -5*um));


  //try joing components
  // TGeoCone* cameraLEDcone = new TGeoCone("cameraLEDcone", (3.77 + 8.02)*cm*0.5, 2.54*cm, 2.54*cm, 2.54*cm, 2.54*cm);
  // TGeoTranslation* myTR1 = new TGeoTranslation("myTR1",(18.78) * TMath::Cos(45.*TMath::DegToRad()) * cm, (18.78) * TMath::Sin(45.*TMath::DegToRad()) * cm, kZf-(3.77 + 8.02)*cm*0.5);
  // TGeoTranslation* myTR2 = new TGeoTranslation("myTR2",0,0,0);
  // TGeoCompositeShape *myCS1 = new TGeoCompositeShape("myCS1", Form("%s:%s", cameraLEDcone->GetName(), myTR1->GetName()));
  // AFocalSurface* idealCamera = new AFocalSurface("idealCamera", myCS1);

  // //AOpticalComponent* cameraLEDobs = new AOpticalComponent("cameraLEDobs", cameraLEDcone);
  // fManager->GetTopVolume()->AddNode(idealCamera, 1);
  //fManager->GetTopVolume()->AddNode(myVol, 2, new TGeoTranslation((18.78) * TMath::Cos(45.*TMath::DegToRad()) * cm, (18.78) * TMath::Sin(45.*TMath::DegToRad()) * cm, kZf-(3.77 + 8.02)*cm*0.5));
  //fManager->GetTopVolume()->AddNode(idealCameraObs, 1, new TGeoTranslation(0, 0, -1*um));

  std::cout << "kZs: " << kZs << " kZf: " << kZf << std::endl;

}

void SST_GATE::AddIdealGroundFocalPlane()
{
  TGeoBBox *idealGroundV = new TGeoBBox("idealGroundV",250*cm,250*cm,0.1*cm);
  
  AFocalSurface* idealGroundFocalArea = new AFocalSurface("idealGroundFocalArea", idealGroundV);
  AObscuration* idealGroundFocalAreaObs = new AObscuration("idealGroundFocalAreaObs", idealGroundV);
  fManager->GetTopVolume()->AddNode(idealGroundFocalArea, 1,new TGeoTranslation(0, 0, 0));
  fManager->GetTopVolume()->AddNode(idealGroundFocalAreaObs, 1, new TGeoTranslation(0, 0, -5*um));
}

//---------  ADD CAMERA BODY ---------//

void SST_GATE::AddCameraBody() 
{
  //This is not the final camera body dimensions rather just a dummy body

  const double kZs = fF/fQ; //position of secondary [cm]
  const double kZf = kZs - (1 - fAlpha)*fF; // where (1 - fAlpha)*fF = separation distance between secondary & focal plane

  // add camera housing obscuration
  //TGeoCone* cameraFrame1V = new TGeoCone("cameraFrame1V", 48.935*cm*0.5, 1.2*fRf*cm, 1.2*fRf*cm, fRf*cm, fRf*cm);
  TGeoTrd2* cameraFrame2 = new TGeoTrd2("cameraFrame2", 1.2*fRf*cm, fRf*cm, 1.2*fRf*cm, fRf*cm, 48.935*cm*0.5);
  //TGeoTube* cameraFrame1V = new TGeoTube("cameraFrame1V", fRf*cm, fRf*cm, 48.935*cm*0.5);
  //AObscuration* cameraFrame1Obs = new AObscuration("cameraFrame1Obs", cameraFrame1V);
  AObscuration* cameraFrame1Obs = new AObscuration("cameraFrame1Obs", cameraFrame2);
  //fManager->GetTopVolume()->AddNodeOverlap(cameraFrame1Obs, 1, new TGeoTranslation(0., 0., kZf*cm - 48.935*cm*0.5 - 1.8*cm));

//Camera section 1
  double xCS1 = 0.5 * 345.*mm;
  double yCS1 = 0.5 * 345.*mm;
  double zCS1 = 0.5 * 64.*mm;
  TGeoBBox* CS1 = new TGeoBBox("CS1", xCS1, yCS1, zCS1);

  double rminLEDhole = 0.*cm;
  double rmaxLEDhole = 0.5 * 5.08*cm;
  double zLEDhole = 0.5 * 1.*cm;
  double ledDeltaZ = zCS1 - 0.5*cm + 0.1*cm;

  TGeoTube* ledHole = new TGeoTube("ledHole", rminLEDhole, rmaxLEDhole, zLEDhole);
  TGeoTranslation* ledHoleT1 = new TGeoTranslation("ledHoleT1", 187.8*mm*TMath::Cos(45*TMath::DegToRad()), 187.8*mm*TMath::Sin(45*TMath::DegToRad()),ledDeltaZ);// kZf - 3.77*cm - 0.5*cm);
  ledHoleT1->RegisterYourself();
  TGeoTranslation* ledHoleT2 = new TGeoTranslation("ledHoleT2", 187.8*mm*TMath::Cos(135*TMath::DegToRad()), 187.8*mm*TMath::Sin(135*TMath::DegToRad()),ledDeltaZ);// kZf - 3.77*cm - 0.5*cm);
  ledHoleT2->RegisterYourself();
  TGeoTranslation* ledHoleT3 = new TGeoTranslation("ledHoleT3", 187.8*mm*TMath::Cos(225*TMath::DegToRad()), 187.8*mm*TMath::Sin(225*TMath::DegToRad()),ledDeltaZ);// kZf - 3.77*cm - 0.5*cm);
  ledHoleT3->RegisterYourself();
  TGeoTranslation* ledHoleT4 = new TGeoTranslation("ledHoleT4", 187.8*mm*TMath::Cos(315*TMath::DegToRad()), 187.8*mm*TMath::Sin(315*TMath::DegToRad()),ledDeltaZ);// kZf - 3.77*cm - 0.5*cm);
  ledHoleT4->RegisterYourself();
  TGeoTranslation* CS1initialT = new TGeoTranslation("CS1initialT", 0.,0.,0.);
  CS1initialT->RegisterYourself();
  //combine elements
  TGeoCompositeShape* combinedCS1 = new TGeoCompositeShape("combinedCS1", Form("((%s:%s)-(%s:%s)-(%s:%s)-(%s:%s)-(%s:%s))", CS1->GetName(), CS1initialT->GetName(), ledHole->GetName(), ledHoleT1->GetName(), ledHole->GetName(), ledHoleT2->GetName(),ledHole->GetName(), ledHoleT3->GetName(),ledHole->GetName(), ledHoleT4->GetName()));

  //AObscuration* cameraCS1combined = new AObscuration("cameraCS1combined",combinedCS1);
  //fManager->GetTopVolume()->AddNodeOverlap(cameraCS1combined, 1, new TGeoTranslation(0.,0.,0.));

  TGeoTranslation* CS1T = new TGeoTranslation("CS1T", 0.0, 0., kZf - 31.*mm - zCS1);
  CS1T->RegisterYourself();

  //Camera section 2
  double zCS2 = 0.5*83.5*mm;
  TGeoArb8 *CS2 = new TGeoArb8("CS2",zCS2);
  CS2->SetVertex(0,220*mm,210.*mm);
  CS2->SetVertex(1,220*mm,-210.*mm);
  CS2->SetVertex(2,-200.*mm,-210.*mm);
  CS2->SetVertex(3,-200.*mm,210.*mm);
  CS2->SetVertex(4,185.*mm,185.*mm);
  CS2->SetVertex(5,185.*mm,-185.*mm);
  CS2->SetVertex(6,-185.*mm,-185.*mm);
  CS2->SetVertex(7,-185.*mm,185.*mm);
  TGeoTranslation* CS2T = new TGeoTranslation("CS2T", 0., 0., kZf - 92.2*mm - zCS2);
  CS2T->RegisterYourself();


  //Camera section 3
  //angle box front
  double zCS3 = 0.5 * 100.5 * mm;

  TGeoArb8 *CS3 = new TGeoArb8("CS3",zCS3);
  CS3->SetVertex(0, 220*mm, 210.*mm);
  CS3->SetVertex(1, 220*mm, -210.*mm);
  CS3->SetVertex(2, -220.*mm, -210.*mm);
  CS3->SetVertex(3, -220*mm, 210.*mm);
  CS3->SetVertex(4, 220.*mm, 210.*mm);
  CS3->SetVertex(5, 220.*mm, -210.*mm);
  CS3->SetVertex(6, -200.*mm, -210.*mm);
  CS3->SetVertex(7, -200.*mm, 210.*mm);
  TGeoTranslation* CS3T = new TGeoTranslation("CS3T", 0., 0., kZf - 175.2*mm - zCS3);
  CS3T->RegisterYourself();
  
  //Camera section 4
  double xCS4 = 0.5*440.*mm;
  double yCS4 = 0.5*420.*mm;
  double zCS4 = 0.5*183.5*mm;//289.*mm;
  double zCS4_Position = kZf - 175.2*mm - zCS4;
  TGeoBBox* CS4 = new TGeoBBox("CS4", xCS4, yCS4, zCS4);
  TGeoTranslation* CS4T = new TGeoTranslation("CS4T", 0., 0., kZf - (175.2*mm+100.5*mm) - zCS4);
  CS4T->RegisterYourself();
  //AObscuration* cameraHousing = new AObscuration("cameraHousing",cameraBox);
  
  //combine elements
  TGeoCompositeShape* cameraBackC1 = new TGeoCompositeShape("cameraBackC1", Form("((%s:%s)+(%s:%s)+(%s:%s)+(%s:%s))", combinedCS1->GetName(), CS1T->GetName(), CS2->GetName(), CS2T->GetName(), CS3->GetName(), CS3T->GetName(),CS4->GetName(), CS4T->GetName()));
  TGeoTranslation* cameraBackC1T = new TGeoTranslation("cameraBackC1T",0.,0.,0.);
  cameraBackC1T->RegisterYourself();



  // CS3->GetName(), CS3T->GetName(), box1->GetName(), box1T->GetName(),CS2->GetName(),CS2T->GetName(),CS1->GetName(), CS1T->GetName(),CS2->GetName()));

  //to add camera body without cooler unit uncomment these 2 lines.
  //AObscuration* cameraHousing = new AObscuration("cameraHousing",cameraBackC1);
  //fManager->GetTopVolume()->AddNodeOverlap(cameraHousing, 1, new TGeoCombiTrans( TGeoTranslation(0., 0., 0.), TGeoRotation("", 90., 0., 0.)));

  //create and combine cooling unit
  //Camera section 7
  double xCS7 = 0.5*39.*mm;
  double yCS7 = 0.5*370.*mm;
  double zCS7 = 0.5*260.*mm;
  double zCooler_Position = kZf - 175.2*mm - zCS7;
  TGeoBBox* coolerbox = new TGeoBBox("coolerbox", xCS7, yCS7, zCS7);
  TGeoTranslation* coolerboxT1 = new TGeoTranslation("coolerboxT1", xCS4+xCS7, 0., zCooler_Position);
  coolerboxT1->RegisterYourself();
  TGeoTranslation* arb3T = new TGeoTranslation("arb3T", 0., 0., 0.);
  arb3T->RegisterYourself();
  //Camera section 8
  double zCS8 = 0.5*260.*mm;
  TGeoArb8 *CS8 = new TGeoArb8("CS8",zCS8);
  CS8->SetVertex(0,10.5*mm, 164.*mm);
  CS8->SetVertex(1,10.5*mm, -164.*mm);
  CS8->SetVertex(2,-10.5*mm, -185.*mm);
  CS8->SetVertex(3,-10.5*mm, 185.*mm);
  CS8->SetVertex(4,10.5*mm, 164.*mm);
  CS8->SetVertex(5,10.5*mm,-164.*mm);
  CS8->SetVertex(6,-10.5*mm,-185.*mm);
  CS8->SetVertex(7,-10.5*mm,185.*mm);
  TGeoTranslation* CS8T = new TGeoTranslation("CS8T", 220.*mm + 39.*mm + 10.5*mm, 0., kZf - 175.2*mm - zCS8);
  CS8T->RegisterYourself();
  TGeoCompositeShape* coolerUnitComb = new TGeoCompositeShape("coolerUnitComb", Form("((%s:%s)+(%s:%s))",coolerbox->GetName(),coolerboxT1->GetName(),CS8->GetName(),CS8T->GetName()));
  TGeoTranslation* coolerUnitCombT = new TGeoTranslation("coolerUnitCombT", 0., 0., 0.);
  coolerUnitCombT->RegisterYourself();

  TGeoCompositeShape* cameraBackC2 = new TGeoCompositeShape("cameraBackC2", Form("((%s:%s)+(%s:%s))",cameraBackC1->GetName(),cameraBackC1T->GetName(),coolerUnitComb->GetName(),coolerUnitCombT->GetName()));

  //to add camera body with cooler unit uncomment these 2 lines.  
  AObscuration* cameraHousing = new AObscuration("cameraHousing",cameraBackC2);
  cameraHousing->RegisterYourself();
  fManager->GetTopVolume()->AddNode(cameraHousing, 1, new TGeoCombiTrans( TGeoTranslation(0., 0., 0.), TGeoRotation("", 90., 0., 0.)));
  
  //fManager->GetTopVolume()->AddNodeOverlap(cameraHousing, 1, new TGeoTranslation(0,0,zbox_Position));
  //  ew TGeoCombiTrans(TGeoTranslation(0, v1.Y()+v4.Y(), v1.Z() + v4.Z()), TGeoRotation("", 0, 180 - theta, 0))
 
  //Camera section 5
  double xCS5 = 0.5 * 450.*mm;
  double yCS5 = 0.5 * 450.*mm;
  double zCS5 = 0.5 * 35.*mm;
  TGeoBBox* CS5 = new TGeoBBox("CS5", xCS5, yCS5, zCS5);
  AObscuration* cameraFixingPlate1 = new AObscuration("cameraFixingPlate1",CS5);
  cameraFixingPlate1->RegisterYourself();
  TGeoTranslation* CS5T = new TGeoTranslation("CS5T", 0.0, 0., kZf - 31.*mm - zCS5);
  CS5T->RegisterYourself();
  fManager->GetTopVolume()->AddNode(cameraFixingPlate1, 1, new TGeoTranslation(0.,0.,kZf - 500.2*mm + zCS5));

  //Camera section 5
  double xCS6 = 0.5 * 480.*mm;
  double yCS6 = 0.5 * 480.*mm;
  double zCS6 = 0.5 * 5.*mm;
  TGeoBBox* CS6 = new TGeoBBox("CS6", xCS6, yCS6, zCS6);
  AObscuration* cameraFixingPlate2 = new AObscuration("cameraFixingPlate2",CS6);
  cameraFixingPlate2->RegisterYourself();
  TGeoTranslation* CS6T = new TGeoTranslation("CS6T", 0.0, 0., kZf - 500.2*mm + 35.*mm + zCS6);
  CS6T->RegisterYourself();
  fManager->GetTopVolume()->AddNode(cameraFixingPlate2, 1, new TGeoTranslation(0.,0.,kZf - 500.2*mm + 35.*mm + zCS6));

  

  //Focal plane test sheet
  TGeoBBox* fpSurf = new TGeoBBox("fpSurf", 0.5*420.*mm, 0.5*420.*mm, 0.1*mm);
  AObscuration* fpSurfObs = new AObscuration("fpSurfObs",fpSurf);
  //fManager->GetTopVolume()->AddNodeOverlap(fpSurfObs, 1, new TGeoTranslation(0,0,kZf));

}

void SST_GATE::AddCameraLid(double deltaPos, double depth)
{
  const double kZs = fF/fQ;                 // fF/fQ = separation distance between primary and secondary
  const double kZf = kZs - (1 - fAlpha)*fF; // (1 - fAlpha)*fF = separation distance between secondary & focal plane
  // //const double kZf = fZf;

  // AGeoAsphericDisk* camera_lid_surface = new AGeoAsphericDisk("camera_lid_surface", kZf + deltaPos, 0, kZf + deltaPos + depth, 0, fRf, 0);
  // //AGeoAsphericDisk* camera_lid_surface = new AGeoAsphericDisk("camera_lid_surface", kZf + 1.*mm, 0, kZf + 1.*um, 0, fRf, 0);
  // //double sagPar[2] = {fKappa1*TMath::Power(fF, -1),
  // //                    fKappa2*TMath::Power(fF, -3)};
  // //idealCameraV->SetPolynomials(2, sagPar, 2, sagPar);

  // double sagPar[5] = {-5.0e-3, -1.25e-7, -6.25e-12, -3.90625e-16, -2.734375e-20}; // Yi focal plane surface
  // camera_lid_surface->SetPolynomials(5, sagPar, 5, sagPar);

  // AFocalSurface* camera_lid_focus = new AFocalSurface("idealCamera", camera_lid_surface);
  // fManager->GetTopVolume()->AddNodeOverlap(camera_lid_focus, 1);
  
  // //AGeoAsphericDisk* focalObsV = new AGeoAsphericDisk("focalObsV", camera_lid_surface->CalcF1(fRf) - 1*cm, 0, camera_lid_surface->CalcF1(fRf), 0, fRf, 0.);  
  // //AObscuration* idealCameraObs = new AObscuration("idealCameraObs", idealCameraV);
  // //AObscuration* idealCameraObs = new AObscuration("idealCameraObs", focalObsV);
  // //fManager->GetTopVolume()->AddNode(idealCameraObs, 1, new TGeoTranslation(0, 0, -1*um));

  //NEW LID
 
  double xLid = 0.5*345.*mm;//365.*mm;
  double yLid = 0.5*53.*mm;//0.5*413.22*mm;
  double zLid = 0.5*413.22*mm;//0.5*53.*mm;
  double zLid_Position = kZf - 31.*mm + zLid;
  TGeoBBox* cameraLid = new TGeoBBox("cameraLid", xLid, yLid, zLid);
  AObscuration* cameraLidObs = new AObscuration("cameraLidObs",cameraLid);
  cameraLidObs->RegisterYourself();
  cameraLidObs->SetFillColor(4);
  //lid open
  fManager->GetTopVolume()->AddNodeOverlap(cameraLidObs, 1,new TGeoCombiTrans( TGeoTranslation(0.,-208.*mm - (0.5*413.22*mm*TMath::Sin(11.1*TMath::DegToRad())),kZf - (0.5*413.22*mm*TMath::Cos(11.1*TMath::DegToRad()))), TGeoRotation("", 0.,-11.1,0.)));//90., -90 + 11.1, 0.)));
  //TGeoCombiTrans( TGeoTranslation( -208.*mm - (0.5*413.22*mm*TMath::Sin(11.1*TMath::DegToRad())),0.,zLid_Position - (0.5*413.22*mm*TMath::Cos(11.1*TMath::DegToRad()))), TGeoRotation("", 90., -90 + 11.1, 0.)));
  // new TGeoTranslation(-208.*mm - (0.5*413.22*mm),0,zLid_Position));
  //lid closed
  //fManager->GetTopVolume()->AddNodeOverlap(cameraLidObs, 1,new TGeoTranslation(0.,0.,zLid_Position));
  double xLidFlange = 0.5*345.*mm;//0.5*25.*mm;//365.*mm;
  double yLidFlange = 0.5*25.*mm;//0.5*345.*mm;
  double zLidFlange = 0.5*5.*mm;//0.5*5.*mm;
  TGeoBBox* cameraLidFlange = new TGeoBBox("cameraLidFlange", xLidFlange, yLidFlange, zLidFlange);
  AObscuration* cameraLidFlangeObs = new AObscuration("cameraLidFlangeObs",cameraLidFlange);
  cameraLidFlangeObs->RegisterYourself();
  fManager->GetTopVolume()->AddNodeOverlap(cameraLidFlangeObs, 1,new TGeoTranslation(0.,-185.*mm+2.*mm,kZf - 28.2*mm));//-185.*mm+1.*mm,0.,kZf - 28.2*mm));





}


//---------  ADD TELESCOPE FRAME ---------//

void SST_GATE::AddTelescopeFrame()
{

  const double kZs = fF/fQ;                 // position of secondary [cm]
  const double kZf = kZs - (1 - fAlpha)*fF; // where (1 - fAlpha)*fF = separation distance between secondary & focal plane
  
  std::cout << " kZs: " << kZs << " kZf: " << kZf << std::endl;

  // Camera support frame
  TVector3 v1(0*m, 290.*mm , kZf*cm - 549.35*mm);    // Point at back of camera edge
  TVector3 v2(0*m, 1068.88*mm, kZf*cm + 300.65*mm);  // Point on secondary edge
  
  TVector3 v3 = v2 - v1;
  TVector3 v4 = v3*0.5;                              // Midpoint of vector v2v1

  TVector3 vFoot1(0*m, 0.*m,  kZf*cm - 549.35*mm);   // Support foot vector
  TVector3 vFoot2(0*m, 290.*mm,  kZf*cm - 549.35*mm);// Support foot vector
  TVector3 vFoot3 = vFoot2 - vFoot1; vFoot3 *= 0.5;
  
  double theta = v4.Theta()*TMath::RadToDeg();
  double phi   = v4.Phi()*TMath::RadToDeg();

  // Primary support frame
  // Note: the structure under the mirror is not as per the actual structure
  // We are only interested in the fixing points as the steel beneath the primary does not affect optical ray tracing.

  TVector3 v1_PS(0*m, 0.*m, 0. - 166.2*mm);   // vector for Secondary Support tube number 1
  TVector3 v2_PS(0*m, -2216.72*mm, 0. - 166.2*mm);
  TVector3 v3_PS = v2_PS - v1_PS;
  TVector3 v4_PS = 0.5 * v3_PS;
  
  // Secondary support frame
  // Note: the structure above the secondary is not complete. However it is not expected to greatly effect the shadowing.

  // v4_PS.Mag()*sPhi3, v4_PS.Mag()*cPhi3
  //TMath::Cos(TMath::DegToRad()*anglePhi2);
  //double sPhi2 = TMath::Sin(TMath::DegToRad()*anglePhi2);
  TVector3 v1_SS1(1958.43*mm, 1177.06*mm, -166.2*mm);
  TVector3 v2_SS1(0*m, 1068.88*mm, kZs+11.967*cm);
  TVector3 v3_SS1 = v2_SS1 - v1_SS1;
  TVector3 v4_SS1 = 0.5 * v3_SS1;
  
  double theta_SS1 = v4_SS1.Theta()*TMath::RadToDeg();
  
  TVector3 v1_SS2(0*m,-v1_SS1.Mag(),-166.2*mm);
  
  TVector3 v2_SS2(1068.88*mm*TMath::Cos(TMath::DegToRad()*30), -1068.88*mm*TMath::Sin(TMath::DegToRad()*30), kZs+11.967*cm);
  TVector3 v3_SS2 = v2_SS2 - v1_SS2;
  TVector3 v4_SS2 = 0.5 * v3_SS2;
  
  double theta_SS2 = v4_SS2.Theta()*TMath::RadToDeg();
  double phi_SS2 = v4_SS2.Phi()*TMath::RadToDeg();
  
  TVector3 v1_SS3(1919.73*mm, 1108.36*mm, -166.2*mm);
  TVector3 v2_SS3(1068.88*mm*TMath::Cos(TMath::DegToRad()*-30), 1068.88*mm*TMath::Sin(TMath::DegToRad()*-30), kZs+11.967*cm);
  TVector3 v3_SS3 = v2_SS3 - v1_SS3;
  TVector3 v4_SS3 = 0.5 * v3_SS3;
  
  double theta_SS3 = v4_SS3.Theta()*TMath::RadToDeg();
  double phi_SS3 = v4_SS3.Phi()*TMath::RadToDeg();
  
  TVector3 v1_SS4(-1919.73*mm, 1108.36*mm, -166.2*mm);
  TVector3 v2_SS4(-1068.88*mm*TMath::Cos(TMath::DegToRad()*-30), 1068.88*mm*TMath::Sin(TMath::DegToRad()*-30), kZs+11.967*cm);
  TVector3 v3_SS4 = v2_SS4 - v1_SS4;
  TVector3 v4_SS4 = 0.5 * v3_SS4;
  
  double theta_SS4 = v4_SS4.Theta()*TMath::RadToDeg();
  double phi_SS4 = v4_SS4.Phi()*TMath::RadToDeg();
  
  
  
  std::cout << std::endl; 
  std::cout << " v1_SS1: x: " << v1_SS1.X() << " y: " << v1_SS1.Y() << " z: " << v1_SS1.Z() << " mag: " << v1_SS1.Mag() << " theta: " << v1_SS1.Theta()*TMath::RadToDeg() << " phi: " << v1_SS1.Phi()*TMath::RadToDeg() << std::endl;
  std::cout << " v2_SS1: x: " << v2_SS1.X() << " y: " << v2_SS1.Y() << " z: " << v2_SS1.Z() << " mag: " << v2_SS1.Mag() << " theta: " << v2_SS1.Theta()*TMath::RadToDeg() << " phi: " << v2_SS1.Phi()*TMath::RadToDeg() << std::endl;
  std::cout << " v3_SS1: x: " << v3_SS1.X() << " y: " << v3_SS1.Y() << " z: " << v3_SS1.Z() << " mag: " << v3_SS1.Mag() << " theta: " << v3_SS1.Theta()*TMath::RadToDeg() << " phi: " << v3_SS1.Phi()*TMath::RadToDeg() << std::endl;
  std::cout << " v4_SS1: x: " << v4_SS1.X() << " y: " << v4_SS1.Y() << " z: " << v4_SS1.Z() << " mag: " << v4_SS1.Mag() << " theta: " << v4_SS1.Theta()*TMath::RadToDeg() << " phi: " << v4_SS1.Phi()*TMath::RadToDeg() << std::endl;

  std::cout << 0. - (v1.Y()+v4.Y())*TMath::Cos(TMath::DegToRad()*30) << " " << 0. - (v1.Y()+v4.Y())*TMath::Sin(TMath::DegToRad()*30) << std::endl;
  std::cout << std::endl;  
  
  //TGeoTube* cameraSupportFrame = new TGeoTube("cameraSupportFrame", 0*mm, 50.*mm, v4.Mag());
  TGeoBBox * cameraSupportFrame = new TGeoBBox("cameraSupportFrame", 25.0*mm, 50.0*mm, v4.Mag());
  AObscuration* cameraSupportFrameObs = new AObscuration("cameraSupportFrameObs", cameraSupportFrame);
  cameraSupportFrameObs->RegisterYourself();  

  TGeoBBox * cameraSupportFrameFoot = new TGeoBBox("cameraSupportFrame", 25.0*mm, 50.0*mm, vFoot3.Mag());
  AObscuration* cameraSupportFrameObs2 = new AObscuration("cameraSupportFrameObs2", cameraSupportFrameFoot);
  cameraSupportFrameObs2->RegisterYourself();
  
  TGeoBBox * primarySupportFrame = new TGeoBBox("primarySupportFrame", 25.0*mm, 125.0*mm, v4_PS.Mag());
  AObscuration* primarySupportFrameObs = new AObscuration("primarySupportFrameObs", primarySupportFrame);
  
  
  //TGeoTorus* secondarySupport = new TGeoTorus(fF/fQ,0,70.*mm,0,360);
  //TGeoBBox * secondarySupport = new TGeoBBox("secondarySupport", 35.*mm, 35.*mm, v4_SS1.Mag());
  TGeoTube* secondarySupport = new TGeoTube("secondarySupport", 0., 35.*mm, v4_SS1.Mag());
  AObscuration* secondarySupportObs = new AObscuration("secondarySupportObs", secondarySupport);
  secondarySupportObs->RegisterYourself();  


  cameraSupportFrameObs->SetLineColor(45);
  cameraSupportFrameObs2->SetLineColor(45);
  primarySupportFrameObs->SetLineColor(45);
  secondarySupportObs->SetLineColor(45);

  // fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs, 1, new TGeoCombiTrans(TGeoTranslation(0, v1.Y()+v4.Y(), v1.Z() + v4.Z()), TGeoRotation("", 0, 180 - theta, 0)));

  for(int i=0; i<3; i++)
    {
      double anglePhi=30;
      double cPhi = TMath::Cos(TMath::DegToRad()*anglePhi);
      double sPhi = TMath::Sin(TMath::DegToRad()*anglePhi);
      
      double anglePhi2=60;
      double cPhi2 = TMath::Cos(TMath::DegToRad()*anglePhi2);
      double sPhi2 = TMath::Sin(TMath::DegToRad()*anglePhi2);
      
      double anglePhi3=60;
      double cPhi3 = TMath::Cos(TMath::DegToRad()*anglePhi2);
      double sPhi3 = TMath::Sin(TMath::DegToRad()*anglePhi2);
      
      TGeoRotation rot1("", -60., theta, 0.);
      //TGeoRotation rot2("", 0, theta, 0);
      
      TGeoRotation rot1_PS("", -30, 0, 0);
      TGeoRotation rot2_PS("", 0, 90, 0);
      
      switch(i)
	{
	case 0:
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs, i + 1, new TGeoCombiTrans(TGeoTranslation(0, v1.Y()+v4.Y(), v1.Z() + v4.Z()), TGeoRotation("", 0, -theta, 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs2, i + 1, new TGeoCombiTrans(TGeoTranslation(0, 0. + vFoot3.Y(), vFoot1.Z()), TGeoRotation("", 0, 90., 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1, new TGeoCombiTrans(TGeoTranslation(0, v4_PS.Y(), v1_PS.Z() + v4_PS.Z()), TGeoRotation("", 0, 90., 0)));
	  break;
	case 1:
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs, i + 1, new TGeoCombiTrans(TGeoTranslation( 0. - (v1.Y()+v4.Y())*cPhi, 0. - (v1.Y()+v4.Y())*sPhi, v1.Z() + v4.Z()), TGeoRotation("", -60., theta, 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs2, i + 1, new TGeoCombiTrans(TGeoTranslation(0. - vFoot3.Y()*cPhi, 0. - vFoot3.Y()*sPhi, vFoot1.Z()), TGeoRotation("", -60., 90., 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1, new TGeoCombiTrans(TGeoTranslation(v4_PS.Mag()*sPhi3, v4_PS.Mag()*cPhi3, v1_PS.Z() + v4_PS.Z()), TGeoRotation("", -60., 90., 0)));
	  break;
	case 2:
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs, i + 1, new TGeoCombiTrans(TGeoTranslation( 0. + (v1.Y()+v4.Y())*cPhi, 0. - (v1.Y()+v4.Y())*sPhi, v1.Z() + v4.Z()), TGeoRotation("", 60., theta, 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs2, i + 1, new TGeoCombiTrans(TGeoTranslation(0. + vFoot3.Y()*cPhi, 0. - vFoot3.Y()*sPhi, vFoot1.Z()), TGeoRotation("", 60., 90., 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1, new TGeoCombiTrans(TGeoTranslation(0. - v4_PS.Mag()*sPhi3, v4_PS.Mag()*cPhi3, v1_PS.Z() + v4_PS.Z()), TGeoRotation("", 60., 90., 0)));
	  break;
	  }
    }
  
  for(int i=0;i<6;i++)
    {
      switch(i)
	{
	case 0:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1, new TGeoCombiTrans(TGeoTranslation(v2_SS1.X() - v4_SS1.X(),v2_SS1.Y() - v4_SS1.Y(),v1_SS1.Z() + v4_SS1.Z()), TGeoRotation("", 90, -theta_SS1, 0)));
	  break;
	case 1:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1, new TGeoCombiTrans(TGeoTranslation(v2_SS1.X() + v4_SS1.X(),v2_SS1.Y() - v4_SS1.Y(),v1_SS1.Z() + v4_SS1.Z()), TGeoRotation("", 90, theta_SS1, 0)));
	  break;
	case 2:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1, new TGeoCombiTrans(TGeoTranslation(v2_SS2.X() - v4_SS2.X(),v1_SS2.Y() + v4_SS2.Y() ,v1_SS2.Z() + v4_SS2.Z()), TGeoRotation("", 90 + phi_SS2, theta_SS2, 0)));//20.863 //68.6223 //27.2994 phi: 62.2111
	  break;
	case 3:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1, new TGeoCombiTrans(TGeoTranslation(-v4_SS2.X(),v1_SS2.Y() + v4_SS2.Y() ,v1_SS2.Z() + v4_SS2.Z()), TGeoRotation("", 90 - phi_SS2 , -theta_SS2, 0)));//27.2994 - 62.2111
	  break;
	case 4:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1, new TGeoCombiTrans(TGeoTranslation(v1_SS3.X() + v4_SS3.X(), v1_SS3.Y() + v4_SS3.Y() ,v1_SS3.Z() + v4_SS3.Z()), TGeoRotation("", phi_SS3 + 90, theta_SS3, 0)));
	  break;
	case 5:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1, new TGeoCombiTrans(TGeoTranslation(v1_SS4.X() + v4_SS4.X(), v1_SS4.Y() + v4_SS4.Y() ,v1_SS4.Z() + v4_SS4.Z()), TGeoRotation("", phi_SS4 + 90, theta_SS4, 0)));
	  break;
	}
    }
}

void SST_GATE::AddTelescopeFrame_version2()
{

  const double kZs = fF/fQ;                 // position of secondary [cm]
  const double kZf = kZs - (1 - fAlpha)*fF; // where (1 - fAlpha)*fF = separation distance between secondary & focal plane
  
  std::cout << " kZs: " << kZs << " kZf: " << kZf << std::endl;

  //-------------------------------
  // Define camera support trusses
  //-------------------------------

  //camera back side
  //TVector3 v1_CS1(0*m*TMath::Cos(90.*TMath::DegToRad()), 380.*mm*TMath::Sin(90.*TMath::DegToRad()) , kZf*cm - 549.35*mm);    // Point at back of camera edge
  //TVector3 v2_CS1(0*m*TMath::Cos(90.*TMath::DegToRad()), 1068.88*mm*TMath::Sin(90. * TMath::DegToRad()), kZf*cm + 300.65*mm);  // Point on secondary edge
  TVector3 v1_CS1(0*m*TMath::Cos(90.*TMath::DegToRad()), 575.66*mm*TMath::Sin(90.*TMath::DegToRad()) , kZf*cm - 643.85*mm);    // Point at back of camera edge
  TVector3 v2_CS1(0*m*TMath::Cos(90.*TMath::DegToRad()), 1469.86*mm*TMath::Sin(90. * TMath::DegToRad()), kZf*cm + 421.81*mm);  // Point on secondary edge
  TVector3 v3_CS1 = v2_CS1 - v1_CS1;
  TVector3 v4_CS1 = v3_CS1 * 0.5;                              // Midpoint of vector v2v1
  double theta_CS1 = v4_CS1.Theta()*TMath::RadToDeg();
  double phi_CS1   = v4_CS1.Phi()*TMath::RadToDeg();

  //camera lid side
  //TVector3 v1_CS2(0*m*TMath::Cos(270.*TMath::DegToRad()), 380.*mm*TMath::Sin(270. * TMath::DegToRad()) , kZf*cm - 549.35*mm);    // Point at back of camera edge
  //TVector3 v2_CS2(0*m*TMath::Cos(270.*TMath::DegToRad()), 1068.88*mm*TMath::Sin(270. * TMath::DegToRad()), kZf*cm + 300.65*mm);  // Point on secondary edge
  TVector3 v1_CS2(0*m*TMath::Cos(270.*TMath::DegToRad()), 575.66*mm*TMath::Sin(270. * TMath::DegToRad()) , kZf*cm - 643.85*mm);    // Point at back of camera edge
  TVector3 v2_CS2(0*m*TMath::Cos(270.*TMath::DegToRad()), 1208.47*mm*TMath::Sin(270. * TMath::DegToRad()), kZf*cm + 110.3*mm);  // Point on secondary edge
  TVector3 v3_CS2 = v2_CS2 - v1_CS2;
  TVector3 v4_CS2 = v3_CS2 * 0.5;                              // Midpoint of vector v2v1
  double theta_CS2 = v4_CS2.Theta()*TMath::RadToDeg();
  double phi_CS2   = v4_CS2.Phi()*TMath::RadToDeg();

  //TGeoTube* cameraSupportFrame = new TGeoTube("cameraSupportFrame", 0*mm, 50.*mm, v4.Mag());
  TGeoBBox * cameraSupportFrame = new TGeoBBox("cameraSupportFrame", 25.0*mm, 50.0*mm, v4_CS1.Mag());//old dims: 25mm x 50mm
  AObscuration* cameraSupportFrameObs = new AObscuration("cameraSupportFrameObs", cameraSupportFrame);
  cameraSupportFrameObs->RegisterYourself();

  TGeoBBox * cameraSupportFrame2 = new TGeoBBox("cameraSupportFrame2", 25.0*mm, 50.0*mm, v4_CS2.Mag());//old dims: 25mm x 50mm 42x75
  AObscuration* cameraSupportFrameObs2 = new AObscuration("cameraSupportFrameObs2", cameraSupportFrame2);
  cameraSupportFrameObs2->RegisterYourself();  

  std::cout << "camera trusses = theta: " << theta_CS1 << "; " << theta_CS2 << "; phi: " << phi_CS1 << "; " << phi_CS2 << "; mag: " <<  v4_CS1.Mag() << "; " << v4_CS2.Mag() << std::endl;
  std::cout << "camera trusses = v1_CS1 X: " << v1_CS1.X() << "; Y: " << v1_CS1.Y() << "; Z: " << v1_CS1.Z() << "; v4_CS1 X: " << v4_CS1.X() << "; Y: " << v4_CS1.Y() << "; Z: " << v4_CS1.Z() << "; v4_CS1 length: " <<  2*v4_CS1.Mag() << std::endl;
  std::cout << "camera trusses = v1_CS2 X: " << v1_CS2.X() << "; Y: " << v1_CS2.Y() << "; Z: " << v1_CS2.Z() << "; v4_CS2 X: " << v4_CS2.X() << "; Y: " << v4_CS2.Y() << "; Z: " << v4_CS2.Z() << "; v4_CS2 length: " <<  2*v4_CS2.Mag() << std::endl;
  
  //-------------------------------
  // Define camera support footings
  //-------------------------------
  
  //TVector3 v1_Foot1(0*m*TMath::Cos(90.*TMath::DegToRad()), 0.*m*TMath::Sin(90.*TMath::DegToRad()),  kZf*cm - 549.35*mm);   // Support foot vector
  //TVector3 v2_Foot1(0*m*TMath::Cos(90.*TMath::DegToRad()), 380.*mm*TMath::Sin(90.*TMath::DegToRad()),  kZf*cm - 549.35*mm);// Support foot vector
  TVector3 v1_Foot1(0*m*TMath::Cos(90.*TMath::DegToRad()), 280.*mm*TMath::Sin(90.*TMath::DegToRad()),  kZf*cm - (611.22 + 25.0)*mm);   // Support foot vector
  TVector3 v2_Foot1(0*m*TMath::Cos(90.*TMath::DegToRad()), 596.63*mm*TMath::Sin(90.*TMath::DegToRad()),  kZf*cm - (611.22 + 25.0)*mm);// Support foot vector
  TVector3 v3_Foot1 = v2_Foot1 - v1_Foot1;
  TVector3 v4_Foot1 = v3_Foot1 * 0.5;
  double theta_Foot1 = v4_Foot1.Theta()*TMath::RadToDeg();
  double phi_Foot1   = v4_Foot1.Phi()*TMath::RadToDeg();
    
  //TVector3 v1_Foot2(0*m*TMath::Cos(270.*TMath::DegToRad()), 0.*m*TMath::Sin(270.*TMath::DegToRad()),  kZf*cm - 549.35*mm);   // Support foot vector
  //TVector3 v2_Foot2(0*m*TMath::Cos(270.*TMath::DegToRad()), 380.*mm*TMath::Sin(270.*TMath::DegToRad()),  kZf*cm - 549.35*mm);// Support foot vector
  TVector3 v1_Foot2(0*m*TMath::Cos(270.*TMath::DegToRad()), 280.*mm*TMath::Sin(270.*TMath::DegToRad()),  kZf*cm - (611.22 + 25.0)*mm);   // Support foot vector
  TVector3 v2_Foot2(0*m*TMath::Cos(270.*TMath::DegToRad()), 596.63*mm*TMath::Sin(270.*TMath::DegToRad()),  kZf*cm - (611.22 + 25.0)*mm);// Support foot vector
  TVector3 v3_Foot2 = v2_Foot2 - v1_Foot2;
  TVector3 v4_Foot2 = v3_Foot2 * 0.5;
  double theta_Foot2 = v4_Foot2.Theta()*TMath::RadToDeg();
  double phi_Foot2   = v4_Foot2.Phi()*TMath::RadToDeg();

  TGeoBBox * cameraSupportFrameFoot = new TGeoBBox("cameraSupportFrame", 25.0*mm, 50.0*mm, v4_Foot1.Mag());//old dims: 25mm x 50mm
  AObscuration* cameraSupportFrameObs3 = new AObscuration("cameraSupportFrameObs3", cameraSupportFrameFoot);
  cameraSupportFrameObs3->RegisterYourself();

  std::cout << "camera footings = theta: " << theta_Foot1 << "; " << theta_Foot2 << "; phi: " << phi_Foot1 << "; " << phi_Foot2 << "; mag: " <<  v4_Foot1.Mag() << "; " << v4_Foot2.Mag() << std::endl;

  //-------------------------------
  // Define camera support bars
  //-------------------------------

  TVector3 v1_CS3(225.*mm, 0.*mm, kZf*cm - 500.2*mm + (0.5*35*mm)); //kZf*cm - 549.35*mm);    // Point on fixing plate
  TVector3 v2_CS3(1068.88*mm, 0.*mm, kZf*cm + 300.65*mm);  // Point on secondary edge
  TVector3 v3_CS3 = v2_CS3 - v1_CS3;
  TVector3 v4_CS3 = v3_CS3 * 0.5;                              // Midpoint of vector v2v1
  double theta_CS3 = v4_CS3.Theta()*TMath::RadToDeg();
  double phi_CS3   = v4_CS3.Phi()*TMath::RadToDeg();

  TVector3 v1_CS4(-225.*mm, 0.*mm, kZf*cm - 500.2*mm + (0.5*35*mm)); // kZf*cm - 549.35*mm);    // Point on fixing plate
  TVector3 v2_CS4(-1068.88*mm, 0.*mm, kZf*cm + 300.65*mm);  // Point on secondary edge
  TVector3 v3_CS4 = v2_CS4 - v1_CS4;
  TVector3 v4_CS4 = v3_CS4 * 0.5;                              // Midpoint of vector v2v1
  double theta_CS4 = v4_CS4.Theta()*TMath::RadToDeg();
  double phi_CS4   = v4_CS4.Phi()*TMath::RadToDeg();
  
  TGeoTube* cameraSupportBar_1 = new TGeoTube("cameraSupportBar_1", 0., 6.*mm, v4_CS3.Mag());
  AObscuration* cameraSupportBarObs_1 = new AObscuration("cameraSupportBarObs_1", cameraSupportBar_1);
  cameraSupportBarObs_1->RegisterYourself();

  TGeoTube* cameraSupportBar_2 = new TGeoTube("cameraSupportBar_2", 0., 6.*mm, v4_CS4.Mag());
  AObscuration* cameraSupportBarObs_2 = new AObscuration("cameraSupportBarObs_2", cameraSupportBar_2);
  cameraSupportBarObs_2->RegisterYourself();

  std::cout << "camera bars = theta: " << theta_CS3 << "; " << theta_CS4 << "; phi: " << phi_CS3 << "; " << phi_CS4 << "; mag: " <<  v4_CS3.Mag() << "; " << v4_CS4.Mag() << std::endl;
  
  //-------------------------------------
  // Define primary mirror support frame
  //-------------------------------------
  
  // Note: the structure under the mirror is not as per the actual structure
  // We are only interested in the fixing points as the steel beneath the primary does not affect optical ray tracing.

  TVector3 v1_PS1(0*m, 0.*m, 0. - 166.2*mm);   // vector for Secondary Support tube number 1
  TVector3 v2_PS1(2216.72*mm*TMath::Cos(45. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(45. * TMath::DegToRad()), -166.2*mm);
  TVector3 v3_PS1 = v2_PS1 - v1_PS1;
  TVector3 v4_PS1 = 0.5 * v3_PS1;
  double phi_PS1 = v4_PS1.Phi()*TMath::RadToDeg();
  
  TVector3 v1_PS2(0*m, 0.*m, 0. - 166.2*mm);   // vector for Secondary Support tube number 1
  TVector3 v2_PS2(2216.72*mm*TMath::Cos(135. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(135. * TMath::DegToRad()), -166.2*mm);
  TVector3 v3_PS2 = v2_PS2 - v1_PS2;
  TVector3 v4_PS2 = 0.5 * v3_PS2;
  double phi_PS2 = v4_PS2.Phi()*TMath::RadToDeg();

  TVector3 v1_PS3(0*m, 0.*m, 0. - 166.2*mm);   // vector for Secondary Support tube number 1
  TVector3 v2_PS3(2216.72*mm*TMath::Cos(225. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(225. * TMath::DegToRad()), -166.2*mm);
  TVector3 v3_PS3 = v2_PS3 - v1_PS3;
  TVector3 v4_PS3 = 0.5 * v3_PS3;
  double phi_PS3 = v4_PS3.Phi()*TMath::RadToDeg();

  TVector3 v1_PS4(0*m, 0.*m, 0. - 166.2*mm);   // vector for Secondary Support tube number 1
  TVector3 v2_PS4(2216.72*mm*TMath::Cos(315. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(315. * TMath::DegToRad()), -166.2*mm);
  TVector3 v3_PS4 = v2_PS4 - v1_PS4;
  TVector3 v4_PS4 = 0.5 * v3_PS4;
  double phi_PS4 = v4_PS4.Phi()*TMath::RadToDeg();

  TGeoBBox * primarySupportFrame = new TGeoBBox("primarySupportFrame", 25.0*mm, 125.0*mm, v4_PS1.Mag());
  AObscuration* primarySupportFrameObs = new AObscuration("primarySupportFrameObs", primarySupportFrame);

  //---------------------------------------
  // Define secondary mirror support masts
  //---------------------------------------
  
  // Note: the structure above the secondary is not complete. However it is not expected to greatly effect the shadowing.

  TVector3 v1_SS1(2216.72*mm*TMath::Cos(45. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(45. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS1(1068.88*mm*TMath::Cos(0. * TMath::DegToRad()),1068.88*mm*TMath::Sin(0. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS1 = v2_SS1 - v1_SS1;
  TVector3 v4_SS1 = 0.5 * v3_SS1;
  
  double theta_SS1 = v4_SS1.Theta()*TMath::RadToDeg();
  double phi_SS1 = v4_SS1.Phi()*TMath::RadToDeg();
 
  TVector3 v1_SS2(2216.72*mm*TMath::Cos(45. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(45. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS2(1068.88*mm*TMath::Cos(90. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(90. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS2 = v2_SS2 - v1_SS2;
  TVector3 v4_SS2 = 0.5 * v3_SS2;
  
  double theta_SS2 = v4_SS2.Theta()*TMath::RadToDeg();
  double phi_SS2 = v4_SS2.Phi()*TMath::RadToDeg();
  
  TVector3 v1_SS3(2216.72*mm*TMath::Cos(135. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(135. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS3(1068.88*mm*TMath::Cos(90. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(90. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS3 = v2_SS3 - v1_SS3;
  TVector3 v4_SS3 = 0.5 * v3_SS3;
  
  double theta_SS3 = v4_SS3.Theta()*TMath::RadToDeg();
  double phi_SS3 = v4_SS3.Phi()*TMath::RadToDeg();
  
  TVector3 v1_SS4(2216.72*mm*TMath::Cos(135. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(135. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS4(1068.88*mm*TMath::Cos(180. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(180. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS4 = v2_SS4 - v1_SS4;
  TVector3 v4_SS4 = 0.5 * v3_SS4;
  
  double theta_SS4 = v4_SS4.Theta()*TMath::RadToDeg();
  double phi_SS4 = v4_SS4.Phi()*TMath::RadToDeg();

  TVector3 v1_SS5(2216.72*mm*TMath::Cos(225. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(225. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS5(1068.88*mm*TMath::Cos(180. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(180. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS5 = v2_SS5 - v1_SS5;
  TVector3 v4_SS5 = 0.5 * v3_SS5;

  double theta_SS5 = v4_SS5.Theta()*TMath::RadToDeg();
  double phi_SS5 = v4_SS5.Phi()*TMath::RadToDeg();

  TVector3 v1_SS6(2216.72*mm*TMath::Cos(225. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(225. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS6(1068.88*mm*TMath::Cos(270. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(270. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS6 = v2_SS6 - v1_SS6;
  TVector3 v4_SS6 = 0.5 * v3_SS6;

  double theta_SS6 = v4_SS6.Theta()*TMath::RadToDeg();
  double phi_SS6 = v4_SS6.Phi()*TMath::RadToDeg();
  
  TVector3 v1_SS7(2216.72*mm*TMath::Cos(315. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(315. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS7(1068.88*mm*TMath::Cos(270. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(270. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS7 = v2_SS7 - v1_SS7;
  TVector3 v4_SS7 = 0.5 * v3_SS7;

  double theta_SS7 = v4_SS7.Theta()*TMath::RadToDeg();
  double phi_SS7 = v4_SS7.Phi()*TMath::RadToDeg();
  
  TVector3 v1_SS8(2216.72*mm*TMath::Cos(315. * TMath::DegToRad()), 2216.72*mm*TMath::Sin(315. * TMath::DegToRad()), -166.2*mm);
  TVector3 v2_SS8(1068.88*mm*TMath::Cos(0. * TMath::DegToRad()), 1068.88*mm*TMath::Sin(0. * TMath::DegToRad()), kZs+11.967*cm);
  TVector3 v3_SS8 = v2_SS8 - v1_SS8;
  TVector3 v4_SS8 = 0.5 * v3_SS8;

  double theta_SS8 = v4_SS8.Theta()*TMath::RadToDeg();
  double phi_SS8 = v4_SS8.Phi()*TMath::RadToDeg();

  //TGeoTorus* secondarySupport = new TGeoTorus(fF/fQ,0,70.*mm,0,360);
  //TGeoBBox * secondarySupport = new TGeoBBox("secondarySupport", 35.*mm, 35.*mm, v4_SS1.Mag());
  TGeoTube* secondarySupport = new TGeoTube("secondarySupport", 0., 35.*mm, v4_SS1.Mag());
  AObscuration* secondarySupportObs = new AObscuration("secondarySupportObs", secondarySupport);
  secondarySupportObs->RegisterYourself();  

  std::cout << "SS1 = theta: " << theta_SS1 << "; " << theta_SS2 << "; " << theta_SS3 << "; " << theta_SS4 << " ; " << theta_SS5 << "; " << theta_SS6 << "; " << theta_SS6 << " ; " << theta_SS8 << "; phi: " << phi_SS1 << "; " << phi_SS2 << "; " << phi_SS3 << "; " << phi_SS4 <<  " ; " << phi_SS5 << " ; "  << phi_SS6 << "; " << phi_SS7 << "; " << phi_SS8 << " ; mag: " <<  v4_SS1.Mag() << "; " << v4_SS2.Mag() << std::endl;
  
  //---------------------------------------
  // Define mast bars
  //---------------------------------------

  TVector3 vBarOne_1((1068.88*mm + 235.88*mm)*TMath::Cos(phi_SS2 * TMath::DegToRad()), (1068.88*mm+235.88*mm)*TMath::Sin(phi_SS2 * TMath::DegToRad()), kZs+11.967*cm-55.16872*cm);
  TVector3 vBarOne_2((1068.88*mm+570.05*mm)*TMath::Cos(phi_SS6 * TMath::DegToRad()), (1068.88*mm+570.05*mm)*TMath::Sin(phi_SS6 * TMath::DegToRad()), kZs+11.967*cm-133.32442*cm);
  TVector3 vBarOne_3 = vBarOne_2 - vBarOne_1;
  TVector3 vBarOne_4 = 0.5 * vBarOne_3;

  double theta_BarOne_1 = vBarOne_1.Theta()*TMath::RadToDeg();
  double phi_BarOne_1 = vBarOne_1.Phi()*TMath::RadToDeg();
  
  TGeoTube* secondarySupportBar_1 = new TGeoTube("secondarySupportBar_1", 0., 35.*mm, vBarOne_4.Mag());
  AObscuration* secondarySupportBarObs_1 = new AObscuration("secondarySupportBarObs_1", secondarySupportBar_1);
  secondarySupportBarObs_1->RegisterYourself();  

  std::cout << "BarOne_1 = theta: " << theta_BarOne_1 << " phi: " << phi_BarOne_1 << " mag: " << vBarOne_1.Mag() << std::endl;
  std::cout << "BarOne_3 = x: " << vBarOne_4.X() << " y: " << vBarOne_4.Y() << " z: " << vBarOne_4.Z() << " mag: " << vBarOne_4.Mag() << " ; " << vBarOne_4.Mag2() << std::endl;
  
  //-------------------------------
  //  Add geometry objects to model
  //-------------------------------
  
  cameraSupportFrameObs->SetLineColor(45);
  cameraSupportFrameObs2->SetLineColor(45);
  cameraSupportFrameObs3->SetLineColor(45);
  cameraSupportBarObs_1->SetLineColor(45);
  cameraSupportBarObs_2->SetLineColor(45);
  primarySupportFrameObs->SetLineColor(45);
  secondarySupportObs->SetLineColor(45);
 
  //----------------------------------------
  //  Add camera support trusses & footings
  //----------------------------------------
  
  for(int i=0; i!=2; i++)
    {
      switch(i)
	{
	case 0:
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_CS1.X() + v4_CS1.X(), v1_CS1.Y() + v4_CS1.Y(), v1_CS1.Z() + v4_CS1.Z()),
								      TGeoRotation("",0,180-theta_CS1, 0)));
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs3, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_Foot1.X() + v4_Foot1.X(), v1_Foot1.Y() + v4_Foot1.Y(), v1_Foot1.Z() + v4_Foot1.Z()),
								      TGeoRotation("", 0., 90., 0.)));
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportBarObs_1, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_CS3.X() + v4_CS3.X(), v1_CS3.Y() + v4_CS3.Y(), v1_CS3.Z() + v4_CS3.Z()),
								      TGeoRotation("",90,theta_CS3, 0)));
	  
	  break;
	case 1:
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs2, i + 1,
	  					   new TGeoCombiTrans(TGeoTranslation(v1_CS2.X() + v4_CS2.X(), v1_CS2.Y() + v4_CS2.Y(), v1_CS2.Z() + v4_CS2.Z()),
	  							      TGeoRotation("", 0, theta_CS2, 0)));//90-theta_CS2
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportFrameObs3, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_Foot2.X() + v4_Foot2.X(), v1_Foot2.Y() + v4_Foot2.Y(), v1_Foot2.Z() + v4_Foot2.Z()),
								      TGeoRotation("", 0., 90., 0.)));
	  fManager->GetTopVolume()->AddNodeOverlap(cameraSupportBarObs_2, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_CS4.X() + v4_CS4.X(), v1_CS4.Y() + v4_CS4.Y(), v1_CS4.Z() + v4_CS4.Z()),
								      TGeoRotation("",-90,theta_CS4, 0)));
	  break;
	  }
    }

  //----------------------------------------
  //  Add primary mirror support structure
  //----------------------------------------

  for(int i=0;i!=4;i++)
    {
      switch(i)
	{
	case 0:
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_PS1.X() + v4_PS1.X(), v1_PS1.Y() + v4_PS1.Y(), v1_PS1.Z() + v4_PS1.Z()),
								      TGeoRotation("", phi_PS1 + 90., 90., 0)));
	  break;
	case 1:
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_PS2.X() + v4_PS2.X(), v1_PS2.Y() + v4_PS2.Y(), v1_PS2.Z() + v4_PS2.Z()),
								      TGeoRotation("", phi_PS2 + 90., 90., 0)));
	  break;
	case 2:
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_PS3.X() + v4_PS3.X(), v1_PS3.Y() + v4_PS3.Y(), v1_PS3.Z() + v4_PS3.Z()),
								      TGeoRotation("", phi_PS3 + 90., 90., 0)));
	  break;
	case 3:
	  fManager->GetTopVolume()->AddNodeOverlap(primarySupportFrameObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_PS4.X() + v4_PS4.X(), v1_PS4.Y() + v4_PS4.Y(), v1_PS4.Z() + v4_PS4.Z()),
								      TGeoRotation("",phi_PS4 + 90., 90., 0)));
	  break;
	}
    }
  
  
  //-------------------------------
  //  Add secondary support masts
  //-------------------------------
  
  for(int i=0;i!=8;i++)
    {
      switch(i)
	{
	case 0:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
	   					   new TGeoCombiTrans(TGeoTranslation(v1_SS1.X() + v4_SS1.X(), v1_SS1.Y() + v4_SS1.Y(), v1_SS1.Z() + v4_SS1.Z()),
	   							      TGeoRotation("", phi_SS1 + 90, theta_SS1, 0)));
	  //fManager->GetTopVolume()->AddNodeOverlap(secondarySupportBarObs_1, i + 1,
	  // 					   new TGeoCombiTrans(TGeoTranslation(vBarOne_4.X(), vBarOne_4.Y(), vBarOne_1.Z() + vBarOne_4.Z()),
	  // 							      TGeoRotation("", 0, 0, 0)));
	  break;
	case 1:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS2.X() + v4_SS2.X(), v1_SS2.Y() + v4_SS2.Y(), v1_SS2.Z() + v4_SS2.Z()),
								      TGeoRotation("", phi_SS2 + 90, theta_SS2, 0)));
	  break;
	case 2:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS3.X() + v4_SS3.X(), v1_SS3.Y() + v4_SS3.Y() ,v1_SS3.Z() + v4_SS3.Z()),
								      TGeoRotation("", phi_SS3 + 90, theta_SS3, 0)));
	  break;
	case 3:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS4.X() + v4_SS4.X(), v1_SS4.Y() + v4_SS4.Y() ,v1_SS4.Z() + v4_SS4.Z()),
								      TGeoRotation("", phi_SS4 + 90 , theta_SS4, 0)));
	  break;
	case 4:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS5.X() + v4_SS5.X(), v1_SS5.Y() + v4_SS5.Y() ,v1_SS5.Z() + v4_SS5.Z()),
								      TGeoRotation("", phi_SS5 + 90, theta_SS5, 0)));
	  break;
	case 5:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS6.X() + v4_SS6.X(), v1_SS6.Y() + v4_SS6.Y() ,v1_SS6.Z() + v4_SS6.Z()),
								      TGeoRotation("", phi_SS6 + 90, theta_SS6, 0)));
	  break;
	  case 6:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS7.X() + v4_SS7.X(), v1_SS7.Y() + v4_SS7.Y() ,v1_SS7.Z() + v4_SS7.Z()),
								      TGeoRotation("", phi_SS7 + 90, theta_SS7, 0)));
	  break;
	  case 7:
	  fManager->GetTopVolume()->AddNodeOverlap(secondarySupportObs, i + 1,
						   new TGeoCombiTrans(TGeoTranslation(v1_SS8.X() + v4_SS8.X(), v1_SS8.Y() + v4_SS8.Y() ,v1_SS8.Z() + v4_SS8.Z()),
								      TGeoRotation("", phi_SS8 + 90, theta_SS8, 0)));
	  break;
	}
    }
}





// void SST_GATE::AddMAPMTFocalPlane()
// {

//   const double kZs = fF/fQ;                 // fF/fQ = separation distance between primary and secondary
//   const double kZf = kZs - (1 - fAlpha)*fF; // (1 - fAlpha)*fF = separation distance between secondary & focal plane
//   //const double kZf = fZf;
  
//   double pixel_WX = 49. * mm;
//   double pixel_WY = 49. * mm;
//   double pixel_DZ = 0.5 * 25.8 * mm;
//   double box_WX = 52. * mm;
//   double box_WY = 52. * mm;
//   double box_DZ = 0.5 * 25.8 * mm;

//   double pix_X = 52.0 * mm; //effective area 49mm
//   double pix_Y = 52.0 * mm; //effective area 49mm
  
//   double pix_diagonal = sqrt(pix_X*pix_X + pix_Y*pix_Y);
//   double pix_start_radius = (2*pix_diagonal) + (0.5*pix_diagonal);
//   double startX = -1.0 * pix_start_radius * TMath::Cos(45*TMath::DegToRad());
//   double startY = pix_start_radius * TMath::Sin(45*TMath::DegToRad()); //a^2 + b^2 = c^2 (find y when x=0)
  
//   std::cout << "pix_diagonal: " << pix_diagonal << " pix_start_radius: " << pix_start_radius << " startX: " << startX << " startY: " << startY << std::endl;
  
  
//   double pix_halfX = 0.5 * pix_X;
//   double pix_halfY = 0.5 * pix_Y;
//   int row = 6;
//   int col = 6;
  
//   std::cout << " pixel_DZ: " << pixel_DZ << std::endl;
  
//   TGeoBBox* pixel_effArea[row][col];
//   TGeoBBox* pixelBox[row][col];
//   TGeoBBox* pixelHousing[row][col];
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  pixel_effArea[i][j] = new TGeoBBox(Form("pixel_effArea_R%d_C%d",i,j), pixel_WX, pixel_WY, pixel_DZ);
// 	  pixelBox[i][j] = new TGeoBBox(Form("pixelBox_R%d_C%d",i,j), pixel_WX, pixel_WY, pixel_DZ);
// 	  pixelHousing[i][j] = new TGeoBBox(Form("pixelHousing_R%d_C%d",i,j), box_WX, box_WY, box_DZ);
// 	}
//     }
  
//   TGeoCompositeShape* spacer[row][col];
//   TGeoTranslation* dummyTranslation = new TGeoTranslation("dummyTranslation", 0., 0., 0.);
//   dummyTranslation->RegisterYourself();
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  spacer[i][j] = new TGeoCompositeShape(Form("spacer_R%d_C%d", i, j), Form("((%s:%s)-(%s:%s))", pixelBox[i][j]->GetName(), dummyTranslation->GetName(), pixelHousing[i][j]->GetName(), dummyTranslation->GetName()));
// 	}
//     }
  
  
//   AObscuration* deadspace[row][col];
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  deadspace[i][j] = new AObscuration(Form("deadspace_R%d_C%d", i, j),spacer[0][0]);
// 	  deadspace[i][j]->SetFillColor(1);
// 	}
//     }
  
//   TGeoTranslation* pixTR1[row][col];
  
//   //double startX = -pix_halfX - (2 * pix_X);
//   //double startY = pix_halfY + (2 * pix_Y);
//   double pixSeparationX = startX;
//   double pixSeparationY = startY;
//   double pixel_Xarray[row][col];
//   double pixel_Yarray[row][col];
  
//   std::cout << "startX: " << startX << " startY: " << startY << " pix_X: " << pix_X << " pix_Y: " << pix_Y << " pix_halfX: " << pix_halfX << " pix_halfY: " << pix_halfY << std::endl;
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  pixTR1[i][j] = new TGeoTranslation(Form("pixTR_R%d_C%d", i, j), pixSeparationX, pixSeparationY, kZf);
// 	  pixTR1[i][j]->RegisterYourself();
// 	  pixel_Xarray[i][j] =  pixSeparationX;
// 	  pixel_Yarray[i][j] = pixSeparationY;
// 	  std::cout << "pixSeparationX: " << pixSeparationX << " pixSeparationY: " << pixSeparationY << std::endl;
// 	  pixSeparationX = pixSeparationX + pix_X;	  
// 	}
//       pixSeparationY = pixSeparationY - pix_Y;
//       pixSeparationX = startX;
//     }
  
//   AGeoAsphericDisk* idealCamera = new AGeoAsphericDisk("idealCamera", kZf - 1*um, 0, kZf, 0, fRf, 0);
//   double sagPar[5] = {-5.0e-3, -1.25e-7, -6.25e-12, -3.90625e-16, -2.734375e-20}; // Yi focal plane surface
//   idealCamera->SetPolynomials(5, sagPar, 5, sagPar);
  
  
//   TVector3 optical_axis;
//   optical_axis.SetMagThetaPhi(kZf,0.0,TMath::Pi()/2);
  
//   TVector3* pixel_positions[row][col];
//   TVector3* pixel_spherical[row][col];
//   double radius[row][col];
//   double rho[row][col];
//   double phi[row][col];
//   double theta[row][col];
  
//   AFocalSurface* pixel[row][col];
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  pixel[i][j] = new AFocalSurface(Form("pixel_R%d_C%d", i, j), pixel_effArea[i][j]);
// 	  pixel[i][j]->RegisterYourself();
// 	  radius[i][j] = sqrt( (pixel_Xarray[i][j]*pixel_Xarray[i][j]) + (pixel_Yarray[i][j]*pixel_Yarray[i][j]) );
// 	  pixel_positions[i][j] = new TVector3(pixel_Xarray[i][j], pixel_Yarray[i][j], idealCamera->CalcF1(radius[i][j]));
// 	}
//     }
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  rho[i][j] = pixel_positions[i][j]->Mag();
// 	  theta[i][j] = TMath::ATan(pixel_positions[i][j]->Y()/pixel_positions[i][j]->X());
// 	  phi[i][j] = (pixel_positions[i][j]->Z()/pixel_positions[i][j]->Mag());
// 	}
//     }
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  pixel_spherical[i][j] = new TVector3(rho[i][j],theta[i][j],phi[i][j]);
// 	}
//     }
  
//   double angle_cosTheta[row][col];
//   double angle_Theta[row][col];
//   double angle_Phi[row][col];
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  angle_cosTheta[i][j] = TMath::ACos((pixel_spherical[i][j]->Dot(optical_axis) / (pixel_spherical[i][j]->Mag()*optical_axis.Mag())))*TMath::RadToDeg();
// 	  angle_Theta[i][j] = pixel_positions[i][j]->Theta()*TMath::RadToDeg();
// 	  angle_Phi[i][j] = pixel_positions[i][j]->Phi()*TMath::RadToDeg();
// 	}
//     }
  
//   std::cout << " angle, r0c1: " << pixel_spherical[0][1]->Angle(optical_axis)*TMath::RadToDeg() << " r3c2: " << pixel_spherical[1][5]->Angle(optical_axis)*TMath::RadToDeg() << std::endl;
//   std::cout << " theta, r0c1: " << pixel_positions[0][1]->Theta()*TMath::RadToDeg() << " " << theta[0][1]*TMath::RadToDeg() << " r3c2: " << pixel_positions[3][2]->Theta()*TMath::RadToDeg() << " " << theta[3][2]*TMath::RadToDeg() << std::endl;
//   std::cout << " phi, r0c1: " << pixel_positions[0][1]->Phi()*TMath::RadToDeg() << " " << phi[0][1]*TMath::RadToDeg() << " r3c2: " << pixel_positions[3][2]->Phi()*TMath::RadToDeg() << " " << phi[3][2]*TMath::RadToDeg() << std::endl;
  
//   TGeoRotation* pixelRot[row][col];
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  pixel[i][j]->RegisterYourself();
// 	  if(i<=2 && j<=2) pixelRot[i][j] = new TGeoRotation("pixelRot",angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
// 	  if(i<=2 && j>=3) pixelRot[i][j] = new TGeoRotation("pixelRot", angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
// 	  if(i>=3 && j<=2) pixelRot[i][j] = new TGeoRotation("pixelRot", angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
// 	  if(i>=3 && j>=3) pixelRot[i][j] = new TGeoRotation("pixelRot", angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
	  
// 	}
//     }
  
//   TGeoCombiTrans* pixelCTR[row][col];
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{
// 	  pixelCTR[i][j] = new TGeoCombiTrans(Form("pixelCTR_R%d_C%d", i, j), pixel_Xarray[i][j], pixel_Yarray[i][j], idealCamera->CalcF1(radius[i][j] - pixel_DZ), pixelRot[i][j]);
// 	  //if(i<=2 && j>=3)pixelCTR[i][j] = new TGeoCombiTrans(Form("pixelCTR_R%d_C%d", i, j), pixel_Xarray[i][j], pixel_Yarray[i][j], idealCamera->CalcF1(radius[i][j]), pixelRot[i][j]);
// 	  //else pixelCTR[i][j] = new TGeoCombiTrans(Form("pixelCTR_R%d_C%d", i, j), pixel_Xarray[i][j], pixel_Yarray[i][j], idealCamera->CalcF1(radius[i][j]), 0);
// 	}
//     }
  
  
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{ 
// 	  if(i==0 && j==0) std::cout << "------- ";// cout << "test 1: " << i << " " << j << endl;
// 	  else if(i==0 && j==5) std::cout << "-------";// cout << "test 2: " << i << " " << j << endl;
// 	  else if(i==5 && j==0) std::cout << "------- ";// cout << "test 3: " << i << " " << j << endl;
// 	  else if(i==5 && j==5) std::cout << "-------";// cout << "test 4: " << i << " " << j << endl;
// 	  else
// 	    {
// 	      fManager->GetTopVolume()->AddNode(pixel[i][j], 1, pixelCTR[i][j]);
// 	      std::cout << pixel_Xarray[i][j] << "," << pixel_Yarray[i][j] <<  " ";
// 	    }
// 	}
//       std::cout << std::endl;
//     }
  
//   std::cout << std::endl;
//   std::cout << std::endl;
//   std::cout << std::endl;  
  
//   for(int i=0; i<row; i++)
//     {
//       for(int j=0; j<col; j++)
// 	{ 
// 	  if(i==0 && j==0) std::cout << "------- ";// cout << "test 1: " << i << " " << j << endl;
// 	  else if(i==0 && j==5) std::cout << "-------";// cout << "test 2: " << i << " " << j << endl;
// 	  else if(i==5 && j==0) std::cout << "------- ";// cout << "test 3: " << i << " " << j << endl;
// 	  else if(i==5 && j==5) std::cout << "-------";// cout << "test 4: " << i << " " << j << endl;
// 	  else
// 	    {
// 	      std::cout << i << "," << j << "(" << angle_cosTheta[i][j] << "," << angle_Theta[i][j] << "," << angle_Phi[i][j] << ") ";
// 	    }
// 	}
//       std::cout << std::endl;
//     }
  
//   std::cout << "kZs: " << kZs << " kZf: " << kZf << std::endl;
  
// }

void SST_GATE::Add_CHEC_MAPMT_FocalPlane(double deltaZ)
{
  const double kZs = fF/fQ;                 // fF/fQ = separation distance between primary and secondary
  const double kZf = kZs - (1 - fAlpha)*fF; // (1 - fAlpha)*fF = separation distance between secondary & focal plane
  double radius_of_curvature = 1000.*mm;    // define focal plane radius of curvature
  double kRp = kZf - radius_of_curvature;               //define the rotation point for the radius of curvature.

  //const double kZf = fZf;

  std::cout << "kZf: " << kZf << " kZf-1000mm: " << kZf - 1000.*mm << std::endl;

//define pixel size (square pixels)
  double pixel_half_length = 5.8*mm * 0.5;
  double pixel_half_height = 1.0*mm * 0.5;
  double pixel_deadspace_inner_half_length = 5.8*mm * 0.5;
  double pixel_deadspace_inner_half_height = 1.2*mm * 0.5;
  double pixel_deadspace_outer_half_length = (5.8*mm + 0.325*mm) * 0.5;
  double pixel_deadspace_outer_half_height = 1.0*mm * 0.5;
  double module_deadspace_inner_half_length = 49.*mm * 0.5;
  double module_deadspace_outer_half_length = 52.*mm * 0.5;
  double module_deadspace_inner_half_height = 2.0*mm * 0.5;//26.8*mm * 0.5;
  double module_deadspace_outer_half_height = 25.8*mm * 0.5;

  TGeoBBox* pixel_volume = new TGeoBBox("pixel_volume", pixel_half_length, pixel_half_length, pixel_half_height);
  TGeoBBox* pixel_deadspace_inner_volume = new TGeoBBox("pixel_deadspace_inner_volume", pixel_deadspace_inner_half_length, pixel_deadspace_inner_half_length, pixel_deadspace_inner_half_height);
  TGeoBBox* pixel_deadspace_outer_volume = new TGeoBBox("pixel_deadspace_outer_volume", pixel_deadspace_outer_half_length, pixel_deadspace_outer_half_length, pixel_deadspace_outer_half_height);
  TGeoBBox* module_deadspace_inner_volume = new TGeoBBox("module_deadspace_inner_volume", module_deadspace_inner_half_length, module_deadspace_inner_half_length, module_deadspace_inner_half_height);
  TGeoBBox* module_deadspace_outer_volume = new TGeoBBox("module_deadspace_outer_volume", module_deadspace_outer_half_length, module_deadspace_outer_half_length, module_deadspace_outer_half_height);

  TGeoTranslation* dummyTranslation = new TGeoTranslation("dummyTranslation", 0., 0., 0.);
  dummyTranslation->RegisterYourself();

  TGeoTranslation* moduleInnerTranslation = new TGeoTranslation("moduleInnerTranslation", 0., 0., module_deadspace_outer_half_height);
  moduleInnerTranslation->RegisterYourself();

  TGeoCompositeShape* pixel_spacer = new TGeoCompositeShape("pixel_spacer", Form("((%s:%s)-(%s:%s))", pixel_deadspace_outer_volume->GetName(), dummyTranslation->GetName(), pixel_deadspace_inner_volume->GetName(), dummyTranslation->GetName()));

  TGeoCompositeShape* module_spacer = new TGeoCompositeShape("module_spacer", Form("((%s:%s)-(%s:%s))", module_deadspace_outer_volume->GetName(), dummyTranslation->GetName(), module_deadspace_inner_volume->GetName(), moduleInnerTranslation->GetName()));

  AObscuration* pixel_deadspace = new AObscuration("deadspace",pixel_spacer);
  pixel_deadspace->SetFillColor(15);

  AObscuration* module_deadspace = new AObscuration("deadspace",module_spacer);
  module_deadspace->SetFillColor(1);


  std::vector< AFocalSurface* > pixelFSA;
  std::vector< std::vector<double> > pixelM;
  std::vector< std::vector< std::vector<double> > > moduleA;
  const int nP=8;
  const int nM=6;
  int mRows = 6;
  int mCols = 6;
  //define the pixel position map: these are half-length multipliers for determining the pixel centrepoint
  double pmX[] = {-7,-5,-3,-1, 1, 3, 5, 7};
  double pmY[] = { 7, 5, 3, 1,-1,-3,-5,-7};
  //define the module position map: these are half-length multipliers for determining the module centre point
  double module_X_posMultiplier[] = {-3,-1,1,3,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-3,-1,1,3};
  double module_Y_posMultiplier[] = {5,5,5,5,3,3,3,3,3,3,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-3,-3,-3,-3,-3,-3,-5,-5,-5,-5};

  std::cout << " module multipliers array sizes: " << std::extent <decltype( module_X_posMultiplier )> ::value << " " << std::extent <decltype( module_Y_posMultiplier )> ::value << std::endl;

  double xP, yP, zP, rP;
  double xM, yM, zM, rM;
  TGeoRotation rb;	  
  TGeoRotation pixelR;
  int pixelCount=0;
  int moduleCount=0;

  AFocalSurface* pixel_focal_surface[32][nP][nP];

  for(int k=0; k!=32; ++k)//loop through modules
    {
      //first get the module centre point
      xM = module_X_posMultiplier[k] * module_deadspace_outer_half_length;
      yM = module_Y_posMultiplier[k] * module_deadspace_outer_half_length;

      for(int i=0; i!=nP; ++i)//loop through pixel rows
	{
	  for(int j=0; j!=nP; ++j)//loop through pixel columns
	    {
	      pixel_focal_surface[k][i][j] = new AFocalSurface(Form("pixel_MN%d_PR%d_PC%d",k,i,j),pixel_volume);
	      //then get the pixel position within the module
	      xP = pmX[j]*pixel_deadspace_outer_half_length;
	      yP = pmY[i]*pixel_deadspace_outer_half_length;
	      zP = radius_of_curvature;
	      pixelR.SetAngles(0,0,0);
	      TGeoTranslation pixelT(xP,yP,zP);
	      TGeoCombiTrans pixelC(pixelT,pixelR);
	      TGeoHMatrix *pixelM = new TGeoHMatrix(pixelC);
	      //Rotate pixel 3.1 degrees about Y-axis and 3.15 degrees about X-axis
	      pixelM->RotateX( module_Y_posMultiplier[k] * 0.5 * (-3.1));
	      pixelM->RotateY( module_X_posMultiplier[k] * 0.5 * 3.15 );
	      //Place rotated pixel at focal plane position in global reference frame
	      double* pixel_coords = pixelM->GetTranslation();
	      pixelM->SetDz(pixel_coords[2] + kRp);//need to move the pixels up to the focal plane position by adding rotation point offset to z-coord
	      fManager->GetTopVolume()->AddNode(pixel_focal_surface[k][i][j], pixelCount, pixelM);
	      fManager->GetTopVolume()->AddNode(pixel_deadspace, pixelCount, pixelM);
	      pixelCount++;
	      //std::cout << "mCols: " << mCols << " i: " << i << " j: " << j << " k: " << k << " moduleCount: " << moduleCount << std::endl;
	      //std::cout << "moduleCount: " << moduleCount << std::endl;
	    }//end pixel columns
	}//end pixel rows
      
      
      zM = radius_of_curvature - module_deadspace_outer_half_height;
      TGeoRotation moduleR;
      moduleR.SetAngles(0,0,0);
      TGeoTranslation moduleT(0,0,zM);//(xM,yM,zM);
      //TGeoTranslation moduleT(0,0,zM);
      TGeoCombiTrans moduleC(moduleT,moduleR);
      TGeoHMatrix *moduleM = new TGeoHMatrix(moduleC);
      //Rotate module 3.1 degrees about Y-axis and 3.15 degrees about X-axis
      moduleM->RotateX( module_Y_posMultiplier[k] * 0.5 * (-3.1)); //module_X_rotation[k]) );
      moduleM->RotateY( module_X_posMultiplier[k] * 0.5 * 3.15 );//module_Y_rotation[k]) );
      //Place rotated module at focal plane position in global reference frame
      double* module_coords = moduleM->GetTranslation();
      moduleM->SetDz(module_coords[2] + kRp);//need to move the modules up to the focal plane position by adding rotation point offset to z-coord
      fManager->GetTopVolume()->AddNode(module_deadspace, moduleCount, moduleM);     
      //std::cout << "mCols: " << mCols << " i: " << i << " j: " << j << " k: " << k << " moduleCount: " << moduleCount << " yM: " << yM << " xM: " << xM << std::endl; 
      moduleCount++;
    }
}


void SST_GATE::AddMAPMTFocalPlane(double deltaZ)
{
  const double kZs = fF/fQ;                 // fF/fQ = separation distance between primary and secondary
  const double kZf = kZs - (1 - fAlpha)*fF; // (1 - fAlpha)*fF = separation distance between secondary & focal plane
  //const double kZf = fZf;

 double innerbox_WX = 49.*mm;
  double innerbox_WY = 49.*mm;
  double innerbox_DZ = 25.8 * mm;//0.5 * 25.8 * mm; 5.*mm;//
  double outerbox_WX = 52. * mm + 6*mm;
  double outerbox_WY = 52. * mm + 6*mm;
  double outerbox_DZ = 0.5 * 25.8 * mm;//0.01*mm;//
  double pix_WX = 49.*mm;
  double pix_WY = 49.*mm;
  double pix_DZ = outerbox_DZ;
 
  std::cout << "box dims: " << innerbox_WX << " " << innerbox_WY << std::endl;
  std::cout << "box dims: " << outerbox_WX << " " << outerbox_WY << std::endl;
  std::cout << "pix dims: " << pix_WX << " " << pix_WY << std::endl;

  double outerbox_halfX = 0.5 * outerbox_WX;
  double outerbox_halfY = 0.5 * outerbox_WY;
  int row = 6;
  int col = 6;

  int nNodesA = 6;
  int nNodesB = 6;
  int boxCount = 0;
  int volCount = 0;
  int nodeCount = 0;
  double rc = 1000*mm; //1m = 1000mm
  double dS = TMath::ASin( outerbox_WX/rc )*TMath::RadToDeg();//3.05778;
  double minTheta = 90 - (2.5 * dS );
  double minPhi =  0 - (2.5 * dS );
  double dTheta = dS ;
  double dPhi = dS ;
  double dx=0., dy=outerbox_DZ, dz=0.;
 

  //void AddAngledBoxSet(TGeoVolume* top, TString name, double width, double height, double depth, double dx=0, double dy=0, double dz=0, int lineCol=1, int fillCol=1)




  //Single:

  TGeoBBox* pixelBox = new TGeoBBox("pixelBox", 0.5*pix_WX, 0.5*pix_WY, pix_DZ);
  TGeoBBox* innerBox = new TGeoBBox("innerBox", 0.5*innerbox_WX, 0.5*innerbox_WY, innerbox_DZ);
  TGeoBBox* outerBox = new TGeoBBox("outerBox", 0.5*outerbox_WX, 0.5*outerbox_WY, outerbox_DZ);
  TGeoBBox* cameraBox = new TGeoBBox("cameraBox", 0.5*450.*mm, 0.5*450.*mm, 48.935*cm*0.5);
  TGeoBBox* baseObscuration = new TGeoBBox("baseObscuration", 0.5*outerbox_WX, 0.5*outerbox_WY, 0.1*cm);


  TGeoTranslation* dummyTranslation = new TGeoTranslation("dummyTranslation", 0., 0., 0.);
  dummyTranslation->RegisterYourself();

  TGeoCompositeShape* spacer = new TGeoCompositeShape("spacer", Form("((%s:%s)-(%s:%s))", outerBox->GetName(), dummyTranslation->GetName(), innerBox->GetName(), dummyTranslation->GetName()));

  AObscuration* deadspace = new AObscuration("deadspace",spacer,0);
  deadspace->SetFillColor(1);
  AObscuration*  base = new AObscuration(" base", baseObscuration,0);
  base->SetFillColor(1);
  AFocalSurface* pixel = new AFocalSurface("pixel", pixelBox);

  std::vector<double> pixel_z_positions;

  for(int j=0; j<nNodesB; j++)
    {
      if(j==0 || j==(nNodesB-1))
	{
	  nNodesA = nNodesB-2; 
	}
      else
	{
	  nNodesA = nNodesB;
	}
      for(int i=0; i<nNodesA; i++)
	{
	  double theta = minTheta + (i*dTheta);//(i * (maxTheta-minTheta) / nNodesA);
	  if(nNodesA ==4)
	    {
	      theta =  minTheta + ((i+1)*dTheta);
	    }
	  double phi = minTheta + (j*dPhi);//(j * (maxPhi-minPhi) / nNodesB);
	  double x = (rc-dy) * sin((theta)*TMath::DegToRad()) * cos((phi)*TMath::DegToRad());
	  double y = ((rc-dy) * sin((theta)*TMath::DegToRad()) * sin((phi)*TMath::DegToRad())) + (kZf - 1*m);
	  double z = ((rc-dy) * cos((theta)*TMath::DegToRad()));
	  TGeoRotation rb;
	  rb.SetAngles(-90+phi,-theta,0);	//rotation defined by Euler angles
	  TGeoTranslation deadspaceT(x,y,z);
	  TGeoTranslation pixelT(x,y,z);
	  TGeoTranslation baseT(x,y-(dy),z);
	  TGeoCombiTrans c1(deadspaceT,rb);
	  TGeoCombiTrans c2(pixelT,rb);
	  TGeoCombiTrans c3(baseT,rb);
	  TGeoHMatrix *deadspaceM = new TGeoHMatrix(c1); 
	  deadspaceM->RotateX(90);
	  TGeoHMatrix *pixelM = new TGeoHMatrix(c2); 
	  pixelM->RotateX(90);
	  TGeoHMatrix *baseM = new TGeoHMatrix(c3); 
	  baseM->RotateX(90);

	  fManager->GetTopVolume()->AddNode(pixel, nodeCount, pixelM);
	  fManager->GetTopVolume()->AddNode(deadspace, nodeCount, deadspaceM);
	  fManager->GetTopVolume()->AddNode(base, nodeCount, baseM);
	  //	  top->AddNode(volumes[volCount],nodeCount, ph);
	  std::cout << i << " " << theta << " " << phi << " " << x << " " << y+dy << " " << z << " " << kZf << " " << kZf - (y+dy) << std::endl;
	  pixel_z_positions.push_back(y);
	  nodeCount ++;
	}
    }

  double lowestZ = *std::min_element(pixel_z_positions.begin(),pixel_z_positions.end());
  double setZ = 0.15*cm + (kZf - lowestZ);

  std::cout << "lowestZ: " << lowestZ << " setZ: " << setZ << std::endl;

  

  /*
  double innerbox_WX = 49.*mm;
  double innerbox_WY = 49.*mm;
  double innerbox_DZ = 5.*mm;//25.8 * mm;//0.5 * 25.8 * mm;
  double outerbox_WX = 52. * mm + 6*mm;
  double outerbox_WY = 52. * mm + 6*mm;
  double outerbox_DZ = 0.01*mm;//0.5 * 25.8 * mm;
  double pix_WX = 49.*mm;
  double pix_WY = 49.*mm;
  double pix_DZ = outerbox_DZ;
 
  std::cout << "box dims: " << innerbox_WX << " " << innerbox_WY << " " << innerbox_DZ << std::endl;
  std::cout << "box dims: " << outerbox_WX << " " << outerbox_WY << " " << outerbox_DZ << std::endl;
  std::cout << "pix dims: " << pix_WX << " " << pix_WY << " " << pix_DZ << std::endl;

  double outerbox_halfX = 0.5 * outerbox_WX;
  double outerbox_halfY = 0.5 * outerbox_WY;
  int row = 6;
  int col = 6;

  //Single:

  TGeoBBox* pixelBox = new TGeoBBox("pixelBox", 0.5*pix_WX, 0.5*pix_WY, pix_DZ);
  TGeoBBox* innerBox = new TGeoBBox("innerBox", 0.5*innerbox_WX, 0.5*innerbox_WY, innerbox_DZ);
  TGeoBBox* outerBox = new TGeoBBox("outerBox", 0.5*outerbox_WX, 0.5*outerbox_WY, outerbox_DZ);

  TGeoTranslation* dummyTranslation = new TGeoTranslation("dummyTranslation", 0., 0., 0.);
  dummyTranslation->RegisterYourself();

  TGeoCompositeShape* spacer = new TGeoCompositeShape("spacer", Form("((%s:%s)-(%s:%s))", outerBox->GetName(), dummyTranslation->GetName(), innerBox->GetName(), dummyTranslation->GetName()));

  AObscuration* deadspace = new AObscuration("deadspace",spacer);
  deadspace->SetFillColor(1);
  AFocalSurface* pixel = new AFocalSurface("pixel", pixelBox);

  //Outer box translations:
  TGeoTranslation* boxTR[row][col];
  
  double box_diagonal = sqrt(outerbox_WX*outerbox_WX + outerbox_WY*outerbox_WY);
  double box_start_radius = (2*box_diagonal) + (0.5*box_diagonal);
  double startX = -1.0 * box_start_radius * TMath::Cos(45*TMath::DegToRad());
  double startY = box_start_radius * TMath::Sin(45*TMath::DegToRad());
  double boxSeparationX = startX;
  double boxSeparationY = startY;
  double box_Xarray[row][col];
  double box_Yarray[row][col];

  std::cout << "startX: " << startX << " startY: " << startY << " outerbox_WX: " << outerbox_WX << " outerbox_WY: " << outerbox_WY << " outerbox_halfX: " << outerbox_halfX << " outerbox_halfY: " << outerbox_halfY << std::endl;

  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  boxTR[i][j] = new TGeoTranslation(Form("boxTR_R%d_C%d", i, j), boxSeparationX, boxSeparationY, kZf);
	  boxTR[i][j]->RegisterYourself();
	  box_Xarray[i][j] =  boxSeparationX;
	  box_Yarray[i][j] = boxSeparationY;
	  std::cout << "boxSeparationX: " << boxSeparationX << " boxSeparationY: " << boxSeparationY << std::endl;
	  boxSeparationX = boxSeparationX + outerbox_WX;	  
	}
      boxSeparationY = boxSeparationY - outerbox_WY;
      boxSeparationX = startX;
    }

  //Inner box / PIXEL translations:
  TGeoTranslation* pixelTR[row][col];
  
  double pixelSeparationX = startX; //pixel starts at the same central position as the outerbox
  double pixelSeparationY = startY; //pixel starts at the same central position as the outerbox
  double pixel_Xarray[row][col];
  double pixel_Yarray[row][col];

  std::cout << "startX: " << startX << " startY: " << startY << " innerbox_WX: " << innerbox_WX << " innerbox_WY: " << innerbox_WY << std::endl;

  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  pixelTR[i][j] = new TGeoTranslation(Form("pixelTR_R%d_C%d", i, j), box_Xarray[i][j], box_Yarray[i][j], kZf);
	  pixelTR[i][j]->RegisterYourself();
	  std::cout << "pixelSeparationX: " << box_Xarray[i][j] << " pixelSeparationY: " << box_Yarray[i][j] << std::endl;
	}
    }

  //implement surface curvature

  AGeoAsphericDisk* idealCamera = new AGeoAsphericDisk("idealCamera", kZf - 1*um, 0, kZf, 0, fRf, 0);
  double sagPar[5] = {-5.0e-3, -1.25e-7, -6.25e-12, -3.90625e-16, -2.734375e-20}; // Yi focal plane surface
  idealCamera->SetPolynomials(5, sagPar, 5, sagPar);
  
  TVector3 optical_axis;
  optical_axis.SetMagThetaPhi(kZf,0.0,TMath::Pi()/2);
  
  TVector3* pixel_positions[row][col];
  TVector3* box_positions[row][col];
  TVector3* pixel_spherical[row][col];
  double radius[row][col];
  double rho[row][col];
  double phi[row][col];
  double theta[row][col];
  double angle_Theta[row][col];
  double angle_Phi[row][col];
  double pixel_Zcorr[row][col];

  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  radius[i][j] = sqrt( (box_Xarray[i][j]*box_Xarray[i][j]) + (box_Yarray[i][j]*box_Yarray[i][j]) );
	  pixel_positions[i][j] = new TVector3(box_Xarray[i][j], box_Yarray[i][j], idealCamera->CalcF1(radius[i][j]));
	  box_positions[i][j] = new TVector3(box_Xarray[i][j], box_Yarray[i][j], idealCamera->CalcF1(radius[i][j]));
	}
    }
  //determine spherical coordindates
  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  rho[i][j] = pixel_positions[i][j]->Mag();
	  theta[i][j] = TMath::ATan(pixel_positions[i][j]->Y()/pixel_positions[i][j]->X());
	  phi[i][j] = (pixel_positions[i][j]->Z()/pixel_positions[i][j]->Mag());
	  angle_Theta[i][j] = pixel_positions[i][j]->Theta()*TMath::RadToDeg();
	  angle_Phi[i][j] = pixel_positions[i][j]->Phi()*TMath::RadToDeg();
	}
    }
  //Determine pixel rotations
  TGeoRotation* pixelRot[row][col];
  
  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  pixel_Zcorr[i][j] = pix_DZ / TMath::Cos(TMath::Abs(angle_Phi[i][j])*TMath::DegToRad());
	  if(i<=2 && j<=2) pixelRot[i][j] = new TGeoRotation("pixelRot",angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
	  if(i<=2 && j>=3) pixelRot[i][j] = new TGeoRotation("pixelRot", angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
	  if(i>=3 && j<=2) pixelRot[i][j] = new TGeoRotation("pixelRot", angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
	  if(i>=3 && j>=3) pixelRot[i][j] = new TGeoRotation("pixelRot", angle_Phi[i][j]-90,-angle_Theta[i][j],-(angle_Phi[i][j]-90));
	  
	}
    }
  //Combine translation and rotation transforms
  TGeoCombiTrans* pixelCTR[row][col];
  double idealFocalPlaneZ=0.0, adjustPixelCentre=0.0, focalPlaneZ=0.0;
  //idealFocalPlaneZ = idealCamera->CalcF1(radius[i][j]);
  //focalPlaneZ = idealFocalPlaneZ + deltaZ; // this will place the top of the central pixel at the focal plane distance. 


  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{
	  pixelCTR[i][j] = new TGeoCombiTrans(Form("pixelCTR_R%d_C%d", i, j), box_Xarray[i][j], box_Yarray[i][j],idealCamera->CalcF1(radius[i][j]), pixelRot[i][j]);
	  //std::cout << "fpZ: " << idealFocalPlaneZ << " deltaZ: " << deltaZ << " finalZ: " << focalPlaneZ << " pixel_Zcorr: " << pixel_Zcorr[i][j] << std::endl;
	}
    }

  int nodeNum=1;

  AObscuration* backplate = new AObscuration("backplate",idealCamera);
  backplate->SetFillColor(4);

  for(int i=0; i<row; i++)
    {
      for(int j=0; j<col; j++)
	{ 
	  if(i==0 && j==0) std::cout << "------- ";// cout << "test 1: " << i << " " << j << endl;
	  else if(i==0 && j==5) std::cout << "-------";// cout << "test 2: " << i << " " << j << endl;
	  else if(i==5 && j==0) std::cout << "------- ";// cout << "test 3: " << i << " " << j << endl;
	  else if(i==5 && j==5) std::cout << "-------";// cout << "test 4: " << i << " " << j << endl;
	  else
	    {
	      //fManager->GetTopVolume()->AddNode(pixel[i][j], 1, pixelCTR[i][j]);
	      //fManager->GetTopVolume()->AddNode(deadspace[i][j], 1, boxTR[i][j]);
	      //fManager->GetTopVolume()->AddNode(pixel, nodeNum, pixelTR[i][j]);
	      //fManager->GetTopVolume()->AddNode(deadspace, nodeNum, boxTR[i][j]);
	      fManager->GetTopVolume()->AddNode(pixel, nodeNum, pixelCTR[i][j]);
	      fManager->GetTopVolume()->AddNode(deadspace, nodeNum, pixelCTR[i][j]);
	      std::cout << box_Xarray[i][j] << "," << box_Yarray[i][j] << " ";
	      nodeNum++;
	    }
	}
      std::cout << std::endl;
    }
  fManager->GetTopVolume()->AddNode(backplate, 1, new TGeoTranslation(0,0,-(idealCamera->CalcF1(radius[5][5]))+(idealCamera->CalcF1(radius[5][5])-4*mm)));
  */
  std::cout << "Successfully added pixels to focal plane." << std::endl;

}

#endif
