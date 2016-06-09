#ifndef SEGMENTED_MIRRORS_H
#define SEGMENTED_MIRRORS_H

#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraph.h"
#include "TPad.h"
#include "TApplication.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TColor.h"
#include "TNtuple.h"

#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoCompositeShape.h"
#include "TGeoPgon.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
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



class SegmentedMirror{
 private:

 protected:
  double fRmin; // The innner radius of a sector
  double fRmax; // The outer radius of a sector
  double fPhimin; // The minimum azimuthal angle
  double fPhimax; // The maximum azimuthal angle 
  double fMargin; // Margin around a mirror
  double fPositionErrorX; // Alignment error along X axis
  double fPositionErrorY; // Alignment error along Y axis
  double fPositionErrorZ; // Alignment error along Z axis
  double fRotationErrorXY; // Rotation error in X-Y plane
  double fRotationErrorSagittal; // Rotation error in sagittal plane
  double fRotationErrorTangential; // Rotation error in tangential plane
  double fRoughness; // Mirror roughness (deg)

 public:
  SegmentedMirror(double rmin, double rmax, double phimin, double phimax);
  virtual ~SegmentedMirror() {}
  virtual AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, bool isPrimary) = 0;
  TGeoCombiTrans* BuildMirrorCombiTrans(AGeoAsphericDisk* disk, bool isPrimary);
  double GetRoughness() const {return fRoughness;}
  void SetMargin(double margin);
  void SetPositionErrors(double x, double y, double z);
  void SetRotationErrors(double xy, double sag, double tan);
  void SetRoughness(double roughness);
};

SegmentedMirror::SegmentedMirror(double rmin, double rmax, double phimin, double phimax)
  : fRmin(rmin), fRmax(rmax), fPhimin(phimin), fPhimax(phimax),
    fMargin(0), fPositionErrorX(0), fPositionErrorY(0), fPositionErrorZ(0),
    fRotationErrorXY(0), fRotationErrorSagittal(0), fRotationErrorTangential(0),
    fRoughness(0)
{
}

TGeoCombiTrans* SegmentedMirror::BuildMirrorCombiTrans(AGeoAsphericDisk* disk, bool isPrimary)
{
  double cphi = (fPhimax + fPhimin)/2.;
  double cr = (fRmax + fRmin)/2.;
  double cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);
  
  TGeoRotation rot1("", -90. + cphi + fRotationErrorXY, fRotationErrorSagittal, 0);
  TGeoRotation rot2("", -cphi, 0., 0.);
  TGeoRotation rot3("", 0., fRotationErrorTangential, 0.);
  TGeoRotation rot4("", cphi, 0., 0.);
  
  TGeoTranslation tr4(cr*TMath::Cos(cphi*TMath::DegToRad()) + fPositionErrorX, cr*TMath::Sin(cphi*TMath::DegToRad()) + fPositionErrorY, cz + fPositionErrorZ);
  
  TGeoCombiTrans* combi = new TGeoCombiTrans(tr4, ((rot4*rot3)*rot2)*rot1);

  return combi;
}

void SegmentedMirror::SetMargin(double margin){
  //set margin around the mirror
  fMargin = margin;
}

void SegmentedMirror::SetPositionErrors(double x, double y, double z){
  //set the position errors along the x,y,z axes
  fPositionErrorX = x;
  fPositionErrorY = y;
  fPositionErrorZ = z;
}

void SegmentedMirror::SetRotationErrors(double xy, double sag, double tan){
  //set the rotation errors in xy,sagittal,tangential planes
  fRotationErrorXY = xy;
  fRotationErrorSagittal = sag;
  fRotationErrorTangential = tan;
}

void SegmentedMirror::SetRoughness(double roughness){
  //set the mirror roughness in units of degrees
  fRoughness = roughness;
}


#endif
