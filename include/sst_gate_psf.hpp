//
//The following code includes a set of functions for calculating and plotting the PSF
//Written by Cameron Rulten 2013-2014 Observatoire de Paris, France.
//
#ifndef SST_GATE_PSF_H
#define SST_GATE_PSF_H

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

class SST_GATE_PSF{
public:
  SST_GATE_PSF();
  ~SST_GATE_PSF();
  
  void calculate_psf(double platescale_cm, std::vector< std::vector<double> > x_coord_array, std::vector< std::vector<double> > y_coord_array,std::vector<double> fov_angles_in, std::string method, std::string ofile_name, std::vector< std::vector<double> > &psf_encircled_out, std::vector< std::vector<double> > &psf_ensquared_out,int centering_method);
  std::pair<double,double> get_encircled_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn, int centering_method);
  std::pair<double,double> get_ensquared_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn, double sqHalfLength, int centering_method);
  double get_mean(std::vector<double> inArray);
  void get_mean_coordinates(std::vector< std::vector<double> > xcoords_all_in, std::vector< std::vector<double> > ycoords_all_in, std::vector<double> &meanX_array_in, std::vector<double> &meanY_array_in);
  void centre_coordinates(std::vector< std::vector<double> > xcoords_all_in, std::vector< std::vector<double> > ycoords_all_in, std::vector< std::vector<double> > &xcoords_centred_in, std::vector< std::vector<double> > &ycoords_centred_in, int centre_method, int tN);
};

SST_GATE_PSF::SST_GATE_PSF()
{
}

SST_GATE_PSF::~SST_GATE_PSF()
{
}

void SST_GATE_PSF::calculate_psf(double platescale_cm, std::vector< std::vector<double> > x_coord_array_centred, std::vector< std::vector<double> > y_coord_array_centred, std::vector<double> fov_angles_in, std::string method, std::string ofile_name, std::vector< std::vector<double> > &psf_encircled_out, std::vector< std::vector<double> > &psf_ensquared_out, int centering_method)
{
  //centering_method must be set to 0 (zero) for the mean-centred method or 1 for the platescale-centred method. Any other value will default to the mean-centering method.
  //basically platescale_cm = gate.deg2dist
  //centre the spot to (0,0) along the x-axis.
  
  //make sure the vectors to be filled are empty
  psf_encircled_out.clear();
  psf_ensquared_out.clear();

  std::vector< std::vector<double> >::iterator xIter;
  std::vector< std::vector<double> >::iterator yIter;
  std::vector<double> psf_encircled_cm;
  std::vector<double> psf_encircled_deg;
  std::vector<double> psf_encircled_arcmin;
  std::vector<double> psf_encircled_error;
  std::vector<double> psf_encircled_photons;
  std::vector<double> psf_ensquared_cm;
  std::vector<double> psf_ensquared_deg;
  std::vector<double> psf_ensquared_arcmin;
  std::vector<double> psf_ensquared_error;
  std::vector<double> psf_ensquared_photons;

  //double psf_encircled_cm[kN], psf_encircled_deg[kN], psf_encircled_arcmin[kN];
  //double psf_encircled_error[kN], psf_encircled_photons[kN];
  double encircled_step=0.0001;
  double ensquared_step=0.0001;
      
  yIter = y_coord_array_centred.begin();
  
  for(xIter=x_coord_array_centred.begin();xIter!=x_coord_array_centred.end();++xIter)
    {
      //cout << "d2size: " << x_coord_array_centred.size() << " " << y_coord_array_centred.size() << " " << (*xIter).size() << " " << (*yIter).size() << endl;
      //calculate encircled psf
      std::pair<double,double> psf_encircled_vals = get_encircled_PSF((*xIter),(*yIter),encircled_step, centering_method); //[cm] [first=psf,second=photons]
      psf_encircled_cm.push_back(psf_encircled_vals.first);                       //[cm]
      psf_encircled_deg.push_back(psf_encircled_vals.first/platescale_cm);        //[deg]
      psf_encircled_arcmin.push_back(psf_encircled_vals.first/platescale_cm*60.0);//[arcmin]
      psf_encircled_photons.push_back(psf_encircled_vals.second);
      psf_encircled_error.push_back(encircled_step/platescale_cm);                //[deg] for [arcmin] * 60.0      
      //calculate ensquared psf
      std::pair<double,double> psf_ensquared_vals = get_ensquared_PSF((*xIter),(*yIter),ensquared_step,0.01, centering_method);//start with a 1mm pixel i.e. 0.1cm
      psf_ensquared_cm.push_back(psf_ensquared_vals.first);                       //[cm]
      psf_ensquared_deg.push_back(psf_ensquared_vals.first/platescale_cm);        //[deg]
      psf_ensquared_arcmin.push_back(psf_ensquared_vals.first/platescale_cm*60.0);//[arcmin]
      psf_ensquared_photons.push_back(psf_ensquared_vals.second);
      psf_ensquared_error.push_back(ensquared_step/platescale_cm);                //[deg] for [arcmin] * 60.0      
      ++yIter;
    }

  //fill output vectors
  //where array item: 0 = cm, 1 = degrees, 2 = armins, 3 = photons, 4 = error on step radius/length
  psf_encircled_out.push_back(psf_encircled_cm);
  psf_encircled_out.push_back(psf_encircled_deg);
  psf_encircled_out.push_back(psf_encircled_arcmin);
  psf_encircled_out.push_back(psf_encircled_photons);
  psf_encircled_out.push_back(psf_encircled_error);
  psf_ensquared_out.push_back(psf_ensquared_cm);
  psf_ensquared_out.push_back(psf_ensquared_deg);
  psf_ensquared_out.push_back(psf_ensquared_arcmin);
  psf_ensquared_out.push_back(psf_ensquared_photons);
  psf_ensquared_out.push_back(psf_ensquared_error);


  std::string output_name;
  int pts = fov_angles_in.size();

  std::cout << std::endl;
  std::cout << "The PSF plotting points total: " << pts << std::endl;
  std::cout << std::endl;

  TGraphErrors* psf_encircled_graph = new TGraphErrors(pts, &fov_angles_in[0], &psf_encircled_deg[0], 0, &psf_encircled_error[0]);//psf_error
  TGraphErrors* psf_ensquared_graph = new TGraphErrors(pts, &fov_angles_in[0], &psf_ensquared_deg[0], 0, &psf_ensquared_error[0]);//psf_error
  TMultiGraph* mg1 = new TMultiGraph();

  TCanvas* plotPSF = new TCanvas("plotPSF", "plotPSF", 900, 700);
  gStyle->SetTitleFontSize(0.04);
  psf_encircled_graph->SetTitle("Evolution of the 80% containment radius with field angle; field angle (degrees);PSF (degrees)");
  psf_encircled_graph->GetXaxis()->SetRangeUser(0.0,4.9);  
  psf_encircled_graph->GetYaxis()->SetRangeUser(0.01,0.1);
  psf_encircled_graph->GetXaxis()->SetTitleSize(0.03);
  psf_encircled_graph->GetYaxis()->SetTitleSize(0.03);
  psf_encircled_graph->GetYaxis()->SetTitleOffset(1.5);
  psf_encircled_graph->GetXaxis()->SetLabelSize(0.028);
  psf_encircled_graph->GetYaxis()->SetLabelSize(0.028);
  psf_encircled_graph->SetMarkerStyle(20);
  psf_encircled_graph->SetMarkerColor(4);
  psf_encircled_graph->Draw("ALP");  
  //mg1->Add(psf_encircled_graph,"ALP");
  psf_ensquared_graph->SetMarkerStyle(21);
  psf_ensquared_graph->SetMarkerColor(3);
  psf_ensquared_graph->Draw("sameLP");
  //mg1->Add(psf_ensquared_graph,"LP");
  //mg1->Draw();
  TGaxis* arcminAxis = new TGaxis(4.9,0.01,4.9,0.1,0.6,6,510,"+L");
  arcminAxis->SetName("arcminAxis");
  arcminAxis->SetTitle("PSF (arcmin)");
  arcminAxis->SetTitleOffset(0.8);
  arcminAxis->SetTitleSize(0.03);
  arcminAxis->SetLabelSize(0.028);
  arcminAxis->SetLabelFont(42);
  arcminAxis->SetTitleFont(42);
  arcminAxis->Draw("same");
  TLegend* leg = new TLegend(0.149554, 0.759287, 0.447545, 0.851412);
  //leg->SetFillStyle(0);
  //leg->SetTextFont(132);
  leg->AddEntry(psf_encircled_graph, "80% encircled (radius)" , "lp");
  leg->AddEntry(psf_ensquared_graph, "80% ensquared (half length)" , "lp");
  leg->Draw("same");
  plotPSF->Update();
  
  output_name = method + "/psf_deg_" + ofile_name + ".png";    //  string output_name;
  plotPSF->SaveAs(output_name.c_str());
  output_name = method + "/psf_deg_" + ofile_name + ".eps";
  plotPSF->SaveAs(output_name.c_str());
  

  
}

//---------  CACULATE PSF SECTION ---------//

std::pair<double,double> SST_GATE_PSF::get_encircled_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn, int centering_method)
{
  //centering_method must be set to 0 (zero) for the mean-centred method or 1 for the platescale-centred method. Any other value will default to the mean-centering method.

  std::vector<double>::iterator x_iter;
  std::vector<double>::iterator y_iter;
  double totalPhotonsX = x_arrayIn.size(); 
  double totalPhotonsY = y_arrayIn.size();

  //cout << "NEW: " << totalPhotonsX << " " << totalPhotonsY << endl;

  double mX = get_mean(x_arrayIn);
  double mY = get_mean(y_arrayIn);
  //std::cout << "encircled mX: " << mX << " mY: " << mY << std::endl;

  //define radius based on mean-centred method or platescale-centred method
  double radius=0.;
  switch(centering_method)
    {
    case 0:
      radius = stepIn;
      break;
    case 1:
      radius = mX + stepIn;
      break;
    default:
      radius = stepIn;
    }

  //  if(centering_method==0) radius = stepIn;
  //else radius = mX + stepIn;

  double ratio = 0.0;
  int photonsInside;

  while(ratio<80)
    {
      double photon_radius=0.0;  
      photonsInside=0;      

      y_iter=y_arrayIn.begin();
      for(x_iter=x_arrayIn.begin();x_iter!=x_arrayIn.end();++x_iter)
	{
	  double x = (*x_iter);
	  double y = (*y_iter);
	  photon_radius = sqrt(x*x + y*y);
	  if(photon_radius<=radius)
	    {
	      //cout << " x: " << x << " y: " << y << " photon_radius: " << photon_radius << endl;
	      photonsInside++;
	    }
	  
	  ++y_iter;
	}
      
      ratio = static_cast<int>(photonsInside/totalPhotonsX*100.0);
      //cout << " photonsInside: " << photonsInside << " photonsInside/totalPhotons :" << ratio << endl;
      radius+=stepIn;
    }

  std::cout << " final radius = " << radius << " photons: " << photonsInside << std::endl;
  std::cout << std::endl;

  return std::make_pair(radius,photonsInside);
}

std::pair<double,double> SST_GATE_PSF::get_ensquared_PSF(std::vector<double> x_arrayIn, std::vector<double> y_arrayIn, double stepIn, double sqHalfLength, int centering_method)
{
  //centering_method must be set to 0 (zero) for the mean-centred method or 1 for the platescale-centred method. Any other value will default to the mean-centering method.
  std::vector<double>::iterator x_iter;
  std::vector<double>::iterator y_iter;
  double totalPhotonsX = x_arrayIn.size(); 
  double totalPhotonsY = y_arrayIn.size();

  //cout << "NEW: " << totalPhotonsX << " " << totalPhotonsY << endl;

  double mX = get_mean(x_arrayIn);
  double mY = get_mean(y_arrayIn);
  //std::cout << "ensquared mX: " << mX << " mY: " << mY << std::endl;

  double squareXminus; // = mX-stepIn;
  double squareXplus;  // = mX+stepIn;
  double squareYminus; // = mY-stepIn;
  double squareYplus;  // = mY+stepIn;

  switch(centering_method)
    {
    case 0:
      squareXminus = 0.0-stepIn;
      squareXplus = 0.0+stepIn;
      squareYminus = 0.0-stepIn;
      squareYplus = 0.0+stepIn;
      break;
    case 1:
      squareXminus = mX-stepIn;
      squareXplus = mX+stepIn;
      squareYminus = mY-stepIn;
      squareYplus = mY+stepIn;
      break;
    default:
      squareXminus = 0.0-stepIn;
      squareXplus = 0.0+stepIn;
      squareYminus = 0.0-stepIn;
      squareYplus = 0.0+stepIn;
    }

  //double radius = mX+stepIn;
  double ratio = 0.0;
  int photonsInside;

  while(ratio<80)
    {
      double photon_radius=0.0;  
      photonsInside=0;      

      y_iter=y_arrayIn.begin();
      for(x_iter=x_arrayIn.begin();x_iter!=x_arrayIn.end();++x_iter)
	{
	  double x = (*x_iter);
	  double y = (*y_iter);
	  switch(centering_method)
	    {
	    case 0:
	      if(x>=(0.0-sqHalfLength) and x<=(0.0+sqHalfLength) and y>=(0.0-sqHalfLength) and y<=(0.0+sqHalfLength))
		{
		  photonsInside++;
		}
	      break;
	    case 1:
	      if(x>=(mX-sqHalfLength) and x<=(mX+sqHalfLength) and y>=(mY-sqHalfLength) and y<=(mY+sqHalfLength))
		{
		  photonsInside++;
		}
	      break;
	    default:
	      if(x>=(0.0-sqHalfLength) and x<=(0.0+sqHalfLength) and y>=(0.0-sqHalfLength) and y<=(0.0+sqHalfLength))
		{
		  photonsInside++;
		}
	    }

	  // if(x>=(mX-sqHalfLength) and x<=(mX+sqHalfLength) and y>=(mY-sqHalfLength) and y<=(mY+sqHalfLength))
	  //   {
	  //     photonsInside++;
	  //   }
	  // // photon_radius = sqrt(x*x + y*y);
	  // // if(photon_radius<=radius)
	  // //   {
	  // //     //cout << " x: " << x << " y: " << y << " photon_radius: " << photon_radius << endl;
	  // //     photonsInside++;
	  // //   }
	  
	  ++y_iter;
	}
      
      ratio = static_cast<int>(photonsInside/totalPhotonsX*100.0);
      //cout << " photonsInside: " << photonsInside << " photonsInside/totalPhotons :" << ratio << endl;
      sqHalfLength+=stepIn;
      //radius+=stepIn;
    }

  std::cout << " final square half length = " << sqHalfLength << " photons: " << photonsInside << std::endl;
  //cout << " final radius = " << radius << " photons: " << photonsInside << endl;
  std::cout << std::endl;

  return std::make_pair(sqHalfLength,photonsInside);
}

double SST_GATE_PSF::get_mean(std::vector<double> inArray)
{
  double sum=0.;
  double mean=0.;

  if(inArray.size()<=0)
    {
      std::cout << "Your array is empty! You cannot perform this calculation!" << std::endl;
      exit(1);
    }
  else
    {
      sum = std::accumulate(inArray.begin(),inArray.end(),0.0);
      mean = sum / inArray.size();
    }
  return mean;
}

void SST_GATE_PSF::get_mean_coordinates(std::vector< std::vector<double> > xcoords_all_in, std::vector< std::vector<double> > ycoords_all_in, std::vector<double> &meanX_array_in, std::vector<double> &meanY_array_in)
{
  std::vector<double>::iterator xcoords_in_iter;
  std::vector<double>::iterator ycoords_in_iter;
  std::vector< std::vector<double> >::iterator xcoords_all_in_iter;
  std::vector< std::vector<double> >::iterator ycoords_all_in_iter;
  double sumX;
  double sumY;
  double meanX;
  double meanY;
   
  ycoords_all_in_iter=ycoords_all_in.begin();
  for(xcoords_all_in_iter=xcoords_all_in.begin(); xcoords_all_in_iter!=xcoords_all_in.end(); ++xcoords_all_in_iter)
    {
      sumX = std::accumulate((*xcoords_all_in_iter).begin(),(*xcoords_all_in_iter).end(),0.0);
      meanX = sumX / (*xcoords_all_in_iter).size();
      sumY = std::accumulate((*ycoords_all_in_iter).begin(),(*ycoords_all_in_iter).end(),0.0);
      meanY = sumY / (*ycoords_all_in_iter).size();
      meanX_array_in.push_back(meanX);
      meanY_array_in.push_back(meanY);
      ycoords_all_in_iter++;
    }
}

void SST_GATE_PSF::centre_coordinates(std::vector< std::vector<double> > xcoords_all_in, std::vector< std::vector<double> > ycoords_all_in, std::vector< std::vector<double> > &xcoords_centred_in, std::vector< std::vector<double> > &ycoords_centred_in, int centre_method, int tN)
{
  //There are two spot centering methods:
  // 0 = mean centred
  // 1 = platescale centred

  double chec_platescale=4.02222222; // cm per degree
  double platescale_corr[] = {0.0, 0.5*chec_platescale,1.0*chec_platescale,1.5*chec_platescale,2.0*chec_platescale,2.5*chec_platescale,3.0*chec_platescale,3.5*chec_platescale,4.0*chec_platescale,4.4*chec_platescale};
  //double platescale_corr[] = {0.5*chec_platescale,1.0*chec_platescale,1.5*chec_platescale,2.0*chec_platescale,2.5*chec_platescale,3.0*chec_platescale,3.5*chec_platescale,4.0*chec_platescale,4.5*chec_platescale,5.0*chec_platescale};
  std::vector<double> newX_array;
  std::vector<double> newY_array;
  double sumX;
  double sumY;
  double meanX;
  double meanY;

  std::vector< std::vector<double> >::iterator xcoords_all_in_iter;
  //std::vector< std::vector<double> >::iterator ycoords_all_in_iter;
  //std::vector<double>::iterator xcoords_in_iter;
  //std::vector<double>::iterator ycoords_in_iter;

  std::vector< std::vector<double> >::iterator ycoords_all_in_iter=ycoords_all_in.begin();
  //for(xcoords_all_in_iter=xcoords_all_in.begin(); xcoords_all_in_iter!=xcoords_all_in.end(); ++xcoords_all_in_iter)
  for(int i=0;i!=tN;++i)
    {
      // sumX = std::accumulate((*xcoords_all_in_iter).begin(),(*xcoords_all_in_iter).end(),0.0);
      // meanX = sumX / (*xcoords_all_in_iter).size();
      // sumY = std::accumulate((*ycoords_all_in_iter).begin(),(*ycoords_all_in_iter).end(),0.0);
      // meanY = sumY / (*ycoords_all_in_iter).size();
      sumX = std::accumulate(xcoords_all_in[i].begin(),xcoords_all_in[i].end(),0.0);
      meanX = sumX / xcoords_all_in[i].size();
      sumY = std::accumulate(ycoords_all_in[i].begin(),ycoords_all_in[i].end(),0.0);
      meanY = sumY / ycoords_all_in[i].size();

      std::vector<double>::iterator xcoords_in_iter;
      //std::vector<double>::iterator ycoords_in_iter = (*ycoords_all_in_iter).begin();
      //for(xcoords_in_iter=(*xcoords_all_in_iter).begin(); xcoords_in_iter!=(*xcoords_all_in_iter).end(); ++xcoords_in_iter)
      std::vector<double>::iterator ycoords_in_iter = ycoords_all_in[i].begin();
      for(xcoords_in_iter=xcoords_all_in[i].begin(); xcoords_in_iter!=xcoords_all_in[i].end(); ++xcoords_in_iter)
       	{
	  double getX = (*xcoords_in_iter);// - meanX;//- meanX_array_in[countA]; //- platescale_corr[countA];//
	  double getY = (*ycoords_in_iter);// - meanY_array_in[countA];
	  double getX_mc = getX-meanX;
	  double getX_psc = getX-platescale_corr[i];
	  if(centre_method==1)
	    {
	      newX_array.push_back(getX_psc);
	      newY_array.push_back(getY);
	    }
	  else
	    {
	      newX_array.push_back(getX_mc);
	      newY_array.push_back(getY);
	    }
	  ycoords_in_iter++;
	}
      xcoords_centred_in.push_back(newX_array);
      ycoords_centred_in.push_back(newY_array);
      newX_array.clear();
      newY_array.clear();
    }

}


#endif
