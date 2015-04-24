// Copyright 2015 Cameron Rulten

//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at

//        http://www.apache.org/licenses/LICENSE-2.0

//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.


// Generic SST-GATE test performance progamme.
// Needs to be changed depending on specific requirements
// To run: ./sst_gate_test 2.0 test_24022014 square 0.0 0.0

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

#include "sst_gate_class.hpp"
#include "sst_gate_globals.hpp"
//#include "sst_gate_test_performance_leds.hpp"
#include "sst_gate_test_performance.hpp"
//#include "sst_gate_test_shadowing.hpp"

using namespace std;

void colourPalette();

int main(int argc, char** argv)
{
  //gStyle->SetPalette(1);

  stringstream get_radius; //get the primary radius
  double radius;
  stringstream get_ofileName; // get the name assigned to the output file
  string ofileName;
  stringstream get_method; // get the ray tracing method
  string method;

  get_radius << argv[1];
  get_radius >> radius;
  get_ofileName << argv[2];
  ofileName = get_ofileName.str();
  get_method << argv[3];
  method = get_method.str();

  double start_pos = strtod(argv[4],NULL); // [mm]
  double depth = strtod(argv[5],NULL);     // [um]

  cout << " running: ./" << argv[0] << " " << radius << " " << ofileName << " " << method << " " << start_pos << " " << depth << endl; 

  //  cout << "New M1 Radius = " << radius << " " << ofileName << endl;

  TApplication *theApp = new TApplication("theApp",&argc,argv);
  
  colourPalette();

  SST_GATE* mySST = new SST_GATE();
  // SST_GATE_LED_PERFORMANCE* test_SST_LEDs = new SST_GATE_LED_PERFORMANCE();
  SST_GATE_PERFORMANCE* test_SST = new SST_GATE_PERFORMANCE();
  //SST_GATE_SHADOWING* test_SSTS = new SST_GATE_SHADOWING();
  
  //const double myRadius=1.0;
  mySST->setNewPrimaryMaxRadius(radius);
  mySST->StartBuildingTelescope();
  mySST->GetManager()->DisableFresnelReflection(1);

  double marginP = 10.*mm;
  double marginS = 0.*mm;
  int set_focal_plane=0; //camera=0, primary=1, secondary=2
  //(if choosing 1 or 2 need to check sst_gate_test_performance.hpp as not fully implemented. i.e. psf should not be calculated)

  AddPrimaryF(*mySST, marginP, false); //set primary as focal plane
  AddSecondaryJ(*mySST, marginS, true, false); //set secondary as monolithic, set secondary as focal plane.

  mySST->AddSecondaryObscuration();
  mySST->AddIdealFocalPlane();
  //mySST->Add_CHEC_MAPMT_FocalPlane(start_pos*mm); //includes pixels
  //mySST->AddMAPMTFocalPlane();
  //mySST->AddMAPMTFocalPlane(start_pos*mm);//, depth*um);

  mySST->AddCameraBody(); //CHEC camera
  mySST->AddCameraLid(start_pos*mm, depth*um); //starting z-position and depth of plane. (must provide unit conversion)
  mySST->AddTelescopeFrame();
  //mySST->AddTelescopeFrame_version2();
  //mySST->AddPrimaryMask();
  //mySST->AddPrimaryDesignFrame();
  //mySST->AddPrimaryTDesignFrame();
  mySST->CloseGeometry();

  cout << "Starting Test Performance..." << endl;
  //test_SST_LEDs->testLEDperformance(*mySST, method, ofileName, false); //set secondary as focal plane
  test_SST->TestPerformance(*mySST, method, ofileName, set_focal_plane); //methods = square or cone
  //test_SSTS->TestShadowing(*mySST, method, ofileName); //methods = square or cone
  // TCanvas *cGeom = new TCanvas("cGeom","cGeom",800,600);
  //mySST->GetManager()->GetTopVolume()->Draw("ogl");

  
  theApp->Run();
  
  cout << "---SST test complete---" << endl;
  
  return 1;
}

void colourPalette()
{
  //example of new colors (greys) and definition of a new palette
  const int NRGBs = 5;
  const int NCont = 999;
  //double stops[NRGBs] = { 0.0, 0.34, 0.61, 0.84, 1.00 };
  //double red[NRGBs] = {0.00, 0.00, 0.87, 1.0, 0.51 };
  //double green[NRGBs] = {0.00, 0.81, 1.0, 0.20, 0.00 };
  //double blue[NRGBs] = {0.51, 1.0, 0.12, 0.00, 0.00 };
  double stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  double red[NRGBs] =   { 0.00, 0.00, 0.87, 1.00, 0.51 };
  double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  double blue[NRGBs] =  { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}
