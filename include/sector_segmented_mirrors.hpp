#ifndef SECTOR_SEGMENTED_MIRRORS_H
#define SECTOR_SEGMENTED_MIRRORS_H

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

#include "segmented_mirrors.hpp"

class SectorSegmentedMirrors : virtual public SegmentedMirror
{
private:
  
public:
  SectorSegmentedMirrors(double rmin, double rmax, double phimin, double phimax);
  ~SectorSegmentedMirrors(){}

  AMirror* BuildMirror(const char* name, AGeoAsphericDisk* disk, bool isPrimary);
};

SectorSegmentedMirrors::SectorSegmentedMirrors(double rmin, double rmax, double phimin, double phimax)
  : SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

AMirror* SectorSegmentedMirrors::BuildMirror(const char* name, AGeoAsphericDisk* disk, bool isPrimary)
{
  double zmin = disk->GetOrigin()[2] - disk->GetDZ();
  double zmax = disk->GetOrigin()[2] + disk->GetDZ();
  double dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  double cr = (fRmax + fRmin)/2.;
  double cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // Define the base of the segment shape
  TGeoTubeSeg* seg1 = new TGeoTubeSeg(Form("%s_seg1", name), fRmin + fMargin, fRmax - fMargin, (zmax - zmin)/2., -dphi + 90., dphi + 90.);

  // To be used for the margin cut of the side edges
  TGeoTubeSeg* seg2 = new TGeoTubeSeg(Form("%s_seg2", name), 0, fRmax - fMargin, (zmax - zmin)/2., -dphi + 90., dphi + 90.);

  TGeoTranslation* tr1 = new TGeoTranslation(Form("%s_tr1", name), 0, -cr, (zmax + zmin)/2. - cz);
  tr1->RegisterYourself();

  double d = fMargin/TMath::Sin(dphi*TMath::DegToRad());

  TGeoTranslation* tr2 = new TGeoTranslation(Form("%s_tr2", name), 0, d - cr, (zmax + zmin)/2. - cz);
  tr2->RegisterYourself();

  TGeoTranslation* tr3 = new TGeoTranslation(Form("%s_tr3", name), 0, -cr, -cz);
  tr3->RegisterYourself();
  
  TGeoCompositeShape* comp = new TGeoCompositeShape(Form("%s_comp", name), Form("%s:%s*%s:%s*%s:%s", seg1->GetName(), tr1->GetName(), seg2->GetName(), tr2->GetName(), disk->GetName(), tr3->GetName()));
  
  AMirror* mirror = new AMirror(Form("%s_mirror", name), comp);
  mirror->SetLineColor(36);
  return mirror;
}

#endif
