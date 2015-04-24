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

//
//The following code includes a set of functions for outputing the optical performance results to data files.
//Written by Cameron Rulten 2013-2014 Observatoire de Paris, France.
//

#ifndef SST_GATE_DATA_OUTPUT_H
#define SST_GATE_DATA_OUTPUT_H

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
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
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

#include "sst_gate_class.hpp"
#include "sst_gate_globals.hpp"
#include "sst_gate_plotter.hpp"

class SST_GATE_DATA_OUTPUT{
private:

public:
  SST_GATE_DATA_OUTPUT();
  ~SST_GATE_DATA_OUTPUT();
  
  void write_data(std::vector< std::vector<double> > psf_encircled_in, std::vector< std::vector<double> > psf_ensquared_in, std::vector<double> fov_angles_in, TH2D* spotH_cm[], double effectiveArea[], double angularResSagittal[], double angularResTangential[], double maxRprimary, std::string method, std::string ofile_name);
};

SST_GATE_DATA_OUTPUT::SST_GATE_DATA_OUTPUT()
{
}

SST_GATE_DATA_OUTPUT::~SST_GATE_DATA_OUTPUT()
{
}

void SST_GATE_DATA_OUTPUT::write_data(std::vector< std::vector<double> > psf_encircled_in, std::vector< std::vector<double> > psf_ensquared_in, std::vector<double> fov_angles_in, TH2D* spotH_cm[], double effectiveArea[], double angularResSagittal[], double angularResTangential[], double maxRprimary, std::string method, std::string ofile_name)
{

  //setup output file
  std::ofstream photonsFile;
  std::string output_name = method + "/data_" + ofile_name + ".txt";
  photonsFile.open(output_name.c_str(), std::ios::app);
  
  //print file headers
  photonsFile << "#FOV(deg): Tot.Phot.: PSF.Phot.: PSF(cm): PSF(deg): PSF(arcmin): EffArea(m2): AngResS(arcmin): AngResT(arcmin): M1radius(cm): Method:" << std::endl;
  
  std::vector<double>::iterator psf_encircled_cm_iter, psf_encircled_deg_iter, psf_encircled_arcmin_iter, psf_encircled_photons_iter, psf_encircled_error_iter;
  std::vector<double>::iterator psf_ensquared_cm_iter, psf_ensquared_deg_iter, psf_ensquared_arcmin_iter, psf_ensquared_photons_iter, psf_ensquared_error_iter;
  std::vector<double>::iterator fov_angles_in_iter;

  //set start iterators where array item: 0 = cm, 1 = degrees, 2 = arcmins, 3 = photons, 4 = error on step radius/length
  psf_encircled_cm_iter=psf_encircled_in[0].begin();
  psf_encircled_deg_iter=psf_encircled_in[1].begin();
  psf_encircled_arcmin_iter=psf_encircled_in[2].begin();
  psf_encircled_photons_iter=psf_encircled_in[3].begin();
  psf_encircled_error_iter=psf_encircled_in[4].begin();
  psf_ensquared_cm_iter=psf_ensquared_in[0].begin();
  psf_ensquared_deg_iter=psf_ensquared_in[1].begin();
  psf_ensquared_arcmin_iter=psf_ensquared_in[2].begin();
  psf_ensquared_photons_iter=psf_ensquared_in[3].begin();
  psf_ensquared_error_iter=psf_ensquared_in[4].begin();
  fov_angles_in_iter=fov_angles_in.begin();

  int countA=0;

  for(int i=0; i!=10;++i)
    {
      //output values to file
      photonsFile << (*fov_angles_in_iter) << "\t" << spotH_cm[i]->GetEntries() << "\t" << (*psf_encircled_photons_iter) << "\t" << (*psf_encircled_cm_iter) << "\t" << (*psf_encircled_deg_iter) << "\t" << (*psf_encircled_arcmin_iter) << "\t" << effectiveArea[i] << "\t" << angularResSagittal[i] << "\t" << angularResTangential[i] << "\t" << maxRprimary << "\t" << method << "\t" << (*psf_ensquared_photons_iter) << "\t" << (*psf_ensquared_cm_iter) << "\t" << (*psf_ensquared_deg_iter) << "\t" << (*psf_ensquared_arcmin_iter) << std::endl;
      
      //output values to screen
      std::cout << "\tpsf: \t\tfov: \tphotons:" << std::endl;
      std::cout << "\t" << (*psf_encircled_cm_iter) << "\t" << (*fov_angles_in_iter) << "\t" << (*psf_encircled_photons_iter) << "\t [cm, deg, number]" << std::endl;
      std::cout << "\t" << (*psf_encircled_deg_iter)<< "\t" << (*fov_angles_in_iter) << "\t" << (*psf_encircled_photons_iter) << "\t [deg, deg, number]" <<std::endl;
      std::cout << "\t" << (*psf_encircled_arcmin_iter) << "\t" << (*fov_angles_in_iter) << "\t" << (*psf_encircled_photons_iter) << "\t [arcmin, deg, number]" <<std::endl;
      std::cout << "\t" << (*psf_ensquared_cm_iter) << "\t" << (*fov_angles_in_iter) << "\t" << (*psf_ensquared_photons_iter) << "\t [cm, deg, number]" << std::endl;
      std::cout << "\t" << (*psf_ensquared_deg_iter)<< "\t" << (*fov_angles_in_iter) << "\t" << (*psf_ensquared_photons_iter) << "\t [deg, deg, number]" <<std::endl;
      std::cout << "\t" << (*psf_ensquared_arcmin_iter) << "\t" << (*fov_angles_in_iter) << "\t" << (*psf_ensquared_photons_iter) << "\t [arcmin, deg, number]" <<std::endl;
      std::cout << "\t" << i << std::endl;
      psf_encircled_cm_iter++;      
      psf_encircled_deg_iter++;
      psf_encircled_arcmin_iter++;
      psf_encircled_photons_iter++;
      psf_encircled_error_iter++;
      psf_ensquared_cm_iter++;
      psf_ensquared_deg_iter++;
      psf_ensquared_arcmin_iter++;
      psf_ensquared_photons_iter++;
      psf_ensquared_error_iter++;
      fov_angles_in_iter++;
    }
  photonsFile.close();



}

#endif
