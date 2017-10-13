/*
The following code includes a set of functions for producing pretty plots
Written by Cameron Rulten 2013 Observatoire de Paris, France.
 */
#ifndef SST_GATE_PLOTTER_H
#define SST_GATE_PLOTTER_H

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
#include "TMultiGraph.h"
#include "TGraph2D.h"
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
#include "TMarker.h"

#include "sst_gate_class.hpp"
#include "sst_gate_globals.hpp"


//using namespace std;

class SST_GATE_PLOTTER{
private:

public:
  SST_GATE_PLOTTER();
  ~SST_GATE_PLOTTER();

  void plot_effArea_angRes(TGraph* graEffArea, TGraph* graSigma[], std::string method, std::string ofile_name, TLegend* leg);
  void create_psf_ellipse(TEllipse ellipseIn[], TH2D* spotH_cm[], int tN, double psf_radius[]);
  void create_pixel_box(TBox* pixelBox[], TH2D* spotH_cm[], int tN, double pixel_half_length);
  void plot_spot(TH2D* spotH_cm[], TBox* pixelBox[], TEllipse* psfEllipse[], int tN, double psf_radius[], std::string method, std::string ofile_name, double xD[], double yD[], double psf_photons[],//
		 double psf[], double effectiveArea[], double angularResSagittal[], double angularResTangential[], double fRpMax, double kDegStep, double fDegStep);
  void plot_spot_ground(TH2D* spotGround_cm[], int tN, std::string method, std::string ofile_name, double kDegStep, double fDegStep);
  void plot_time_xFP(TProfile* hTime_xFP[], int tN);
  void plot_time_yFP(TProfile* hTime_yFP[], int tN);
  void plot_time_xyS(TProfile2D* hTime_xyS[], int tN, std::string method, std::string ofile_name);
  void plot_time_xyP(TProfile2D* hTime_xyP[], int tN, std::string method, std::string ofile_name);
  void plot_time_xyFP_S(TProfile2D* hTime_xyFP_S[], int tN, std::string method, std::string ofile_name);
  void plot_time_xyFP_Total(TProfile2D* hTime_xyFP_Total[], int tN, std::string method, std::string ofile_name);
  //  void plot_time_xyFP_Total(TProfile2D* hTime_xyFP_Total[], int tN);
  void plot_time_xyFP_Delay(TProfile2D* hTime_xyFP_Delay[], int tN, std::string method, std::string ofile_name);
  void plot_mean_time_delay(int tN, double fov_angles[], double delay_times_mean[], double delay_times_spread[], std::string method, std::string ofile_name);
  //void plot_onAxis(TProfile2D* hist_xy[], std::string method, std::string ofile_name);
  //void plot_onAxis(TProfile2D* hist_xy[], std::string method, std::string ofile_name, int axisAngle);
  void plot_onAxis(TProfile2D* hist_xy[], std::string method, std::string ofile_name, int axisAngle, std::string canvasName);
  void plot_rays(const SST_GATE& gate_in, std::vector<TPolyLine3D*> poly_in);
  void plot_rays_multi(const SST_GATE& gate_in, std::vector< std::vector<TPolyLine3D*> > poly_in);
  void plot_focal_plane(std::vector< std::vector<double> > x_arrayIn, std::vector< std::vector<double> > y_arrayIn, std::vector< std::vector<double> > z_arrayIn);
  void plot_shadowing(double totalRaysIn[], double sCameraIn[], double sLidIn[], double sMastsIn[], double sTrussesIn[], double sSecondaryIn[], double fov_angles[], int tN, std::string method, std::string ofile_name);
  void plot_on_v_off(TProfile2D* hist_xy[], std::string method, std::string ofile_name, int on_axis, int off_axis, int mirror_index, double fieldStep, double fDegStep);
};

SST_GATE_PLOTTER::SST_GATE_PLOTTER()
{
}

SST_GATE_PLOTTER::~SST_GATE_PLOTTER()
{
}

void plot_effArea_angRes(TGraph* graEffArea, TGraph* graSigma[], std::string method, std::string ofile_name, TLegend* leg)
{
  TCanvas* canFig5 = new TCanvas("canFig5", "canFig5", 1200, 600);
  canFig5->Divide(2, 1);
  canFig5->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  graEffArea->Draw("apl");
  graEffArea->SetMarkerStyle(25);
  graEffArea->GetXaxis()->SetLimits(0, 5);
  graEffArea->GetYaxis()->SetRangeUser(0, 20);
  

  // PSF is not consistent with the original paper, but the spot diagram at
  // 5 (deg) is consistent with each other by eye comparison. There may be a
  // difference between calculations of RMS in my code and the paper
  canFig5->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  graSigma[0]->Draw("apl");
  graSigma[1]->Draw("pl same");
  graSigma[0]->SetMarkerStyle(25);
  graSigma[1]->SetMarkerStyle(20);
  graSigma[0]->GetXaxis()->SetLimits(0, 5);
  graSigma[0]->GetYaxis()->SetRangeUser(0, 5);
  leg->Draw();

  std::string output_name = method + "/effArea_AngRes_" + ofile_name + ".pdf";

  canFig5->SaveAs(output_name.c_str());
}

void create_psf_ellipse(TEllipse* ellipseIn[], TH2D* spotH_cm[], int tN, double psf_radius[])
{

  double xD[tN], yD[tN];
  
  for(int i=0;i<tN;i++)
    {
      double meanX = spotH_cm[i]->GetMean(1);
      double meanY = spotH_cm[i]->GetMean(2);
      ellipseIn[i] = new TEllipse(meanX,meanY,psf_radius[i],psf_radius[i]);
    }

}

void create_pixel_box(TBox* pixelBox[], TH2D* spotH_cm[], int tN, double pixel_half_length)
{

  for(int i=0;i<tN;i++)
    {
      double meanX = spotH_cm[i]->GetMean(1);
      double meanY = spotH_cm[i]->GetMean(2);
      pixelBox[i] = new TBox(meanX-pixel_half_length*mm, meanY-pixel_half_length*mm, meanX+pixel_half_length*mm, meanY+pixel_half_length*mm);
    }
}

void plot_spot(TH2D* spotH_cm[], TBox* pixelBox[], TEllipse* psfEllipse[], int tN, double psf_radius[], std::string method, std::string ofile_name, double xD[], double yD[], double psf_photons[], double psf[], double effectiveArea[], double angularResSagittal[], double angularResTangential[], double fRpMax, double kDegStep, double fDegStep)
{

  //gStyle->SetPalette(53);
  
  std::string output_name;
  // ofstream photonsFile;
  // std::string output_name = method + "/data_" + ofile_name + ".txt";
  // photonsFile.open(output_name.c_str(),ios::app);
  
  
  //TCanvas* plotSpot_cm = new TCanvas("plotSpot_cm", "plotSpot_cm", 1917, 796);
  TCanvas* plotSpot_cm = new TCanvas("plotSpot_cm", "plotSpot_cm", 1440, 500);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.06);
  //TCanvas* plotSpot_cm = new TCanvas("plotSpot_cm", "plotSpot_cm", 3508, 2480); //optimised for A4 landscape 300dpi
  plotSpot_cm->Divide(5,2);
  TPad *pads[tN];
  std::stringstream get_padNum;
  std::string padNum;
  TPaveStats *mystats[tN];
  TPaveStats *newStats;

  for(int i=0;i!=tN;++i)
    {
      pads[i] = (TPad*)plotSpot_cm->GetPad(i+1);
    }

  for(int i=0;i!=tN;++i)
    {
      pads[i]->SetMargin(0.15,0.15,0.15,0.15);
      //pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }

  plotSpot_cm->Update();

  double angle=0.0; //reset angle
  

  //photonsFile << "# FOVangle(deg): \t photons: \t photons: \t psf(cm): \t effArea(m2): \t angResS(arcmin): \t angResT(arcmin): \t M1radius(cm): \t method:" << std::endl;
  for(int i=0;i!=tN;++i)
    {
      //std::cout << "spot " << i << " title: " << spotH_cm[i]->GetTitle() << std::endl;
      //change pad
      pads[i]->cd();
      pads[i]->SetLogz();
      spotH_cm[i]->GetZaxis()->SetRangeUser(1,spotH_cm[i]->GetMaximum());
      spotH_cm[i]->GetXaxis()->SetLabelSize(0.058);
      spotH_cm[i]->GetYaxis()->SetLabelSize(0.058);
      spotH_cm[i]->GetZaxis()->SetLabelSize(0.058);
      spotH_cm[i]->GetXaxis()->SetLabelOffset(0.01);
      spotH_cm[i]->GetYaxis()->SetLabelOffset(0.01);
      spotH_cm[i]->GetXaxis()->SetTitleSize(0.058);
      spotH_cm[i]->GetYaxis()->SetTitleSize(0.058);
      spotH_cm[i]->GetXaxis()->SetTitleOffset(1.2);
      spotH_cm[i]->GetYaxis()->SetTitleOffset(1.2);
      spotH_cm[i]->SetLineColor(1);
      spotH_cm[i]->SetMarkerStyle(6);
      spotH_cm[i]->SetMarkerColor(4);
      spotH_cm[i]->Draw("colz");
      //To draw zero value bins as lowest colour uncomment this line.
      //spotH[i]->SetMinimum(-1); 
      //std::cout << "mean: " << spotH[i]->GetMean(1) << std::endl;
      xD[i] = spotH_cm[i]->GetMean(1);
      yD[i] = angle;

      //optimise labels and titles
      // spotH_cm[i]->SetLabelSize(0.0425,"X:Y");
      // spotH_cm[i]->SetLabelSize(0.045,"Z");
      // spotH_cm[i]->SetLabelOffset(0.009,"X");
      // spotH_cm[i]->SetLabelOffset(0.008,"Y");
      // spotH_cm[i]->SetLabelOffset(0.008,"Z");
      // spotH_cm[i]->SetTitleSize(0.045,"X:Y");
      // spotH_cm[i]->SetTitleOffset(1.,"X");
      // spotH_cm[i]->SetTitleOffset(1.15,"Y");
     

      //optimise axis range (only when not full camera)
      //spotH_cm[i]->SetAxisRange(spotH_cm[i]->GetMean(1)-0.7,spotH_cm[i]->GetMean(2)+0.7,"X;Y");
      //spotH_cm[i]->SetAxisRange(spotH_cm[i]->GetMean(1)-0.6,spotH_cm[i]->GetMean(1)+0.6,"X");
      //spotH_cm[i]->SetAxisRange(spotH_cm[i]->GetMean(2)-0.6,spotH_cm[i]->GetMean(2)+0.6,"Y");
      
      //draw spot
      //spotH_cm[i]->Draw("colz");

      // //output values to file
      // photonsFile << angle << "\t" << spotH_cm[i]->GetEntries() << "\t" << psf_photons[i] << "\t" << psf[i] << "\t" << effectiveArea[i] << "\t" << angularResSagittal[i] << "\t" << angularResTangential[i] << "\t" << fRpMax << "\t" << method << std::endl;

      //draw pixel box
      pixelBox[i]->SetLineColor(4);
      pixelBox[i]->SetLineWidth(1);
      pixelBox[i]->SetLineStyle(1);
      pixelBox[i]->SetFillStyle(0);
      pixelBox[i]->Draw("sames");
      
      //draw PSF circle
      psfEllipse[i]->SetLineColor(7);
      psfEllipse[i]->SetLineWidth(1);
      psfEllipse[i]->SetFillStyle(0);
      psfEllipse[i]->Draw("sames");

      //step through FOV angles
      if(i==9) angle+=fDegStep;
      else angle+=kDegStep;
    }

  plotSpot_cm->Modified();
  plotSpot_cm->Update();

  // //reposition stats boxes
  // for(int i=0;i<tN;i++)
  //   {
  //     pads[i]->cd();
  //     get_padNum << i+1;
  //     padNum = "stats_" + get_padNum.str();     
  //     newStats = (TPaveStats*)pads[i]->GetPrimitive("stats");
  //     newStats->SetName(padNum.c_str());  
  //     newStats->SetX1NDC(0.213861);
  //     newStats->SetY1NDC(0.689887);
  //     newStats->SetX2NDC(0.734114);
  //     newStats->SetY2NDC(0.947591);
  //     get_padNum.str("");
  //   }

  //close files and save plots
  //photonsFile.close();
  output_name = method + "/spot_" + ofile_name + ".eps";
  plotSpot_cm->SaveAs(output_name.c_str());
  output_name = method + "/spot_" + ofile_name + ".png";
  plotSpot_cm->SaveAs(output_name.c_str());

}

void plot_spot_ground(TH2D* spotGround_cm[], int tN, std::string method, std::string ofile_name, double kDegStep, double fDegStep)
{

  //gStyle->SetPalette(53);
  
  std::string output_name;
  
  //TCanvas* plotSpot_cm = new TCanvas("plotSpot_cm", "plotSpot_cm", 1917, 796);
  TCanvas* plotSpot_ground = new TCanvas("plotSpot_ground", "plotSpot_ground", 1440, 500);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.06);
  //TCanvas* plotSpot_ground = new TCanvas("plotSpot_ground", "plotSpot_ground", 3508, 2480); //optimised for A4 landscape 300dpi
  plotSpot_ground->Divide(5,2);
  TPad *pads[tN];
  std::stringstream get_padNum;
  std::string padNum;
  TPaveStats *mystats[tN];
  TPaveStats *newStats;

  for(int i=0;i!=tN;++i)
    {
      pads[i] = (TPad*)plotSpot_ground->GetPad(i+1);
    }

  for(int i=0;i!=tN;++i)
    {
      pads[i]->SetMargin(0.15,0.15,0.15,0.15);
      //pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }

  plotSpot_ground->Update();

  double angle=0.0; //reset angle
  

  //photonsFile << "# FOVangle(deg): \t photons: \t photons: \t psf(cm): \t effArea(m2): \t angResS(arcmin): \t angResT(arcmin): \t M1radius(cm): \t method:" << std::endl;
  for(int i=0;i!=tN;++i)
    {
      //std::cout << "spot " << i << " title: " << spotGround_cm[i]->GetTitle() << std::endl;
      //change pad
      pads[i]->cd();
      pads[i]->SetLogz();
      spotGround_cm[i]->GetZaxis()->SetRangeUser(1,spotGround_cm[i]->GetMaximum());
      spotGround_cm[i]->GetXaxis()->SetLabelSize(0.058);
      spotGround_cm[i]->GetYaxis()->SetLabelSize(0.058);
      spotGround_cm[i]->GetZaxis()->SetLabelSize(0.058);
      spotGround_cm[i]->GetXaxis()->SetLabelOffset(0.01);
      spotGround_cm[i]->GetYaxis()->SetLabelOffset(0.01);
      spotGround_cm[i]->GetXaxis()->SetTitleSize(0.058);
      spotGround_cm[i]->GetYaxis()->SetTitleSize(0.058);
      spotGround_cm[i]->GetXaxis()->SetTitleOffset(1.2);
      spotGround_cm[i]->GetYaxis()->SetTitleOffset(1.2);
      spotGround_cm[i]->SetLineColor(1);
      spotGround_cm[i]->SetMarkerStyle(6);
      spotGround_cm[i]->SetMarkerColor(4);
      spotGround_cm[i]->Draw("colz");
      //To draw zero value bins as lowest colour uncomment this line.
      //spotH[i]->SetMinimum(-1); 

      //step through FOV angles
      if(i==9) angle+=fDegStep;
      else angle+=kDegStep;
    }

  plotSpot_ground->Modified();
  plotSpot_ground->Update();

  output_name = method + "/ground_" + ofile_name + ".eps";
  plotSpot_ground->SaveAs(output_name.c_str());
  output_name = method + "/ground_" + ofile_name + ".png";
  plotSpot_ground->SaveAs(output_name.c_str());

}





// hTime_xFP[i]->Fill(px3,pt3);                        // photon position in x direction and mean time taken from secondary mirror
// hTime_yFP[i]->Fill(py3,pt3);                        // photon position in y direction and mean time taken from secondary mirror
// hTime_xyFP_S[i]->Fill(px3-x,py3,(pt3-pt2)*1.e9);    // photon position on focal plane and mean time taken from secondary mirror
// hTime_xyS[i]->Fill(px2,py2,(pt2-pt1)*1.e9);         // photon position on secondary mirror and mean time taken from primary mirror
// hTime_xyFP_Total[i]->Fill(px3-x,py3,(pt3-pt1)*1.e9);// photon position on focal plane and mean time taken from primary mirror
// hTime_xyFP_Delay[i]->Fill(px3-x,py3,(pt3-pt1)*1.e9);// photon position on focal plane and mean delay between photon arrival time

void plot_time_xFP(TProfile* hTime_xFP[], int tN)
{
  TCanvas* plotTime_xFP = new TCanvas("plotTime_xFP", "plotTime_xFP", 3508, 2480); //optimised for A4 landscape 300dpi
  plotTime_xFP->Divide(5,2);

  for(int i=0;i<tN;i++)
    {
      //optimise labels and titles
      hTime_xFP[i]->SetMinimum(1.5);
      hTime_xFP[i]->SetMaximum(3.0);
      hTime_xFP[i]->SetLabelSize(0.0425,"X:Y");
      hTime_xFP[i]->SetLabelSize(0.045,"Z");
      hTime_xFP[i]->SetLabelOffset(0.009,"X");
      hTime_xFP[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_xFP[i]->SetTitleSize(0.045,"X:Y");
      hTime_xFP[i]->SetTitleOffset(1.,"X");
      hTime_xFP[i]->SetTitleOffset(1.1,"Y");
      plotTime_xFP->cd(i + 1);
      hTime_xFP[i]->Draw();
    }
}

void plot_time_yFP(TProfile* hTime_yFP[], int tN)
{
  TCanvas* plotTime_yFP = new TCanvas("plotTime_yFP", "plotTime_yFP", 3508, 2480); //optimised for A4 landscape 300dpi
  plotTime_yFP->Divide(5,2);

  for(int i=0;i<tN;i++)
    {
      //optimise labels and titles
      hTime_yFP[i]->SetMinimum(1.5);
      hTime_yFP[i]->SetMaximum(3.0);
      hTime_yFP[i]->SetLabelSize(0.0425,"X:Y");
      hTime_yFP[i]->SetLabelSize(0.045,"Z");
      hTime_yFP[i]->SetLabelOffset(0.009,"X");
      hTime_yFP[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_yFP[i]->SetTitleSize(0.045,"X:Y");
      hTime_yFP[i]->SetTitleOffset(1.,"X");
      hTime_yFP[i]->SetTitleOffset(1.1,"Y");
      plotTime_yFP->cd(i + 1);
      hTime_yFP[i]->Draw();
    }
}

void plot_time_xyFP_S(TProfile2D* hTime_xyFP_S[], int tN, std::string method, std::string ofile_name)
{
  TCanvas* plotTime_xyFP_S = new TCanvas("plotTime_xyFP_S", "plotTime_xyFP_S", 3508, 2480); //optimised for A4 landscape 300dpi
  plotTime_xyFP_S->Divide(5,2);
  TPad *pads[tN];

  for(int i=0;i<tN;i++)
    {
      pads[i] = (TPad*)plotTime_xyFP_S->GetPad(i+1);
    }
  
  for(int i=0;i<tN;i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }
  plotTime_xyFP_S->Update();

  for(int i=0;i<tN;i++)
    {
      pads[i]->cd();
      //optimise labels and titles
      hTime_xyFP_S[i]->SetStats(0);
      hTime_xyFP_S[i]->SetLabelSize(0.04,"X:Y:Z");
      //hTime_xyFP_S[i]->SetLabelSize(0.045,"Z");
      hTime_xyFP_S[i]->SetLabelOffset(0.009,"X");
      hTime_xyFP_S[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_xyFP_S[i]->SetTitleSize(0.045,"X:Y:Z");
      hTime_xyFP_S[i]->SetTitleOffset(1.,"X");
      hTime_xyFP_S[i]->SetTitleOffset(1.1,"Y");
      hTime_xyFP_S[i]->Draw("colz");
    }
  
  std::string output_name = method + "/time_StoFP_" + ofile_name + ".pdf";
  
  plotTime_xyFP_S->SaveAs(output_name.c_str());
}

void plot_time_xyFP_Total(TProfile2D* hTime_xyFP_Total[], int tN, std::string method, std::string ofile_name)
{
  //TCanvas* plotTime_xyFP_Total = new TCanvas("plotTime_xyFP_Total", "plotTime_xyFP_Total", 3508, 2480); //optimised for A4 landscape 300dpi
  //TCanvas* plotTime_xyFP_Total = new TCanvas("plotTime_xyFP_Total", "plotTime_xyFP_Total", 1790, 1240); //optimised for A5 landscape 300dpi
  TCanvas* plotTime_xyFP_Total = new TCanvas("plotTime_xyFP_Total", "plotTime_xyFP_Total", 1917, 796);
  plotTime_xyFP_Total->Divide(5,2);
  TPad *pads[tN];

  for(int i=0;i<tN;i++)
    {
      pads[i] = (TPad*)plotTime_xyFP_Total->GetPad(i+1);
    }
  
  for(int i=0;i<tN;i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }
  plotTime_xyFP_Total->Update();

  for(int i=0;i<tN;i++)
    {
      pads[i]->cd();
      //optimise labels and titles
      hTime_xyFP_Total[i]->SetStats(0);
      hTime_xyFP_Total[i]->SetLabelSize(0.04,"X:Y:Z");
      //hTime_xyFP_Total[i]->SetLabelSize(0.045,"Z");
      hTime_xyFP_Total[i]->SetLabelOffset(0.009,"X");
      hTime_xyFP_Total[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_xyFP_Total[i]->SetTitleSize(0.045,"X:Y:Z");
      hTime_xyFP_Total[i]->SetTitleOffset(1.,"X");
      hTime_xyFP_Total[i]->SetTitleOffset(1.1,"Y");
      hTime_xyFP_Total[i]->Draw("colz");
    }

  std::string output_name = method + "/time_xyFP_total_" + ofile_name + ".png";
  
  plotTime_xyFP_Total->SaveAs(output_name.c_str());

}


void plot_time_xyFP_Delay(TProfile2D* hTime_xyFP_Delay[], int tN, std::string method, std::string ofile_name)
{
  TCanvas* plotTime_xyFP_Delay = new TCanvas("plotTime_xyFP_Delay", "plotTime_xyFP_Delay", 3508, 2480); //optimised for A4 landscape 300dpi

  plotTime_xyFP_Delay->Divide(5,2);
  TPad *pads[tN];

  for(int i=0;i<tN;i++)
    {
      pads[i] = (TPad*)plotTime_xyFP_Delay->GetPad(i+1);
    }
  
  for(int i=0;i<tN;i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }
  plotTime_xyFP_Delay->Update();

  for(int i=0;i<tN;i++)
    {
      pads[i]->cd();
      //optimise labels and titles
      hTime_xyFP_Delay[i]->SetStats(0);
      hTime_xyFP_Delay[i]->SetLabelSize(0.04,"X:Y:Z");
      //hTime_xyFP_Delay[i]->SetLabelSize(0.045,"Z");
      hTime_xyFP_Delay[i]->SetLabelOffset(0.009,"X");
      hTime_xyFP_Delay[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_xyFP_Delay[i]->SetTitleSize(0.045,"X:Y:Z");
      hTime_xyFP_Delay[i]->SetTitleOffset(1.,"X");
      hTime_xyFP_Delay[i]->SetTitleOffset(1.1,"Y");
      hTime_xyFP_Delay[i]->Draw("colz");
    }

  std::string output_name = method + "/time_delay_spot_" + ofile_name + ".pdf";
  
  plotTime_xyFP_Delay->SaveAs(output_name.c_str());
}

void plot_mean_time_delay(int tN, double fov_angles[], double delay_times_mean[], double delay_times_spread[],std::string method, std::string ofile_name)
{
  TGraphErrors* gr_delay = new TGraphErrors(tN, fov_angles, delay_times_mean, 0, delay_times_spread);
  TCanvas* c_gr_delay = new TCanvas("c_gr_delay","c_gr_delay",800,600);
  c_gr_delay->SetMargin(0.15, 0.85, 0.15, 0.85);
  gr_delay->SetTitle("; field angle (degrees); time (ns)");//Mean photon delay between first and last photon arrivals
  gr_delay->GetXaxis()->SetLabelSize(0.04);
  gr_delay->GetYaxis()->SetLabelSize(0.04);
  //gr_delay->GetXaxis()->SetLabelOffset(0.008);
  //gr_delay->GetYaxis()->SetLabelOffset(0.008);
  gr_delay->GetXaxis()->SetTitleSize(0.04);
  gr_delay->GetYaxis()->SetTitleSize(0.04);
  gr_delay->GetXaxis()->SetTitleOffset(1.1);
  gr_delay->GetYaxis()->SetTitleOffset(1.1);
  gr_delay->SetMarkerColor(4);
  gr_delay->SetMarkerStyle(20);
  gr_delay->GetXaxis()->SetLimits(0.0,5.0);
  gr_delay->GetYaxis()->SetRangeUser(0.0,1.0);
  gr_delay->Draw("ALPE");
  c_gr_delay->Update();

  std::string output_name = method + "/time_delay_fov_" + ofile_name + "_STDEV.eps";
  
  c_gr_delay->SaveAs(output_name.c_str());
  output_name = method + "/time_delay_fov_" + ofile_name + "_STDEV.pdf";
  c_gr_delay->SaveAs(output_name.c_str());
  
}

void plot_time_xyS(TProfile2D* hTime_xyS[], int tN, std::string method, std::string ofile_name)
{
  //TCanvas* plotTime_xyS = new TCanvas("plotTime_xyS", "plotTime_xyS", 3508, 2480); //optimised for A4 landscape 300dpi
  //TCanvas* plotTime_xyS = new TCanvas("plotTime_xyS", "plotTime_xyS", 1790, 1240); //optimised for A5 landscape 300dpi
  TCanvas* plotTime_xyS = new TCanvas("plotTime_xyS", "plotTime_xyS", 1440, 500);
  plotTime_xyS->Divide(5,2);
  TPad *pads[tN];

  for(int i=0;i<tN;i++)
    {
      pads[i] = (TPad*)plotTime_xyS->GetPad(i+1);
    }
  
  for(int i=0;i<tN;i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }
  plotTime_xyS->Update();

  for(int i=0;i<tN;i++)
    {
      pads[i]->cd();
      pads[i]->SetLogz();
      //optimise labels and titles
      hTime_xyS[i]->SetStats(0);
      hTime_xyS[i]->SetLabelSize(0.04,"X:Y:Z");
      //hTime_xyS[i]->SetLabelSize(0.045,"Z");
      hTime_xyS[i]->SetLabelOffset(0.009,"X");
      hTime_xyS[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_xyS[i]->SetTitleSize(0.045,"X:Y:Z");
      hTime_xyS[i]->SetTitleOffset(1.,"X");
      hTime_xyS[i]->SetTitleOffset(1.1,"Y");
      hTime_xyS[i]->Draw("colz");
    }

  std::string output_name = method + "/time_secondary_mirror_" + ofile_name + ".pdf";
  
  plotTime_xyS->SaveAs(output_name.c_str());
}

void plot_time_xyP(TProfile2D* hTime_xyP[], int tN, std::string method, std::string ofile_name)
{
  //TCanvas* plotTime_xyP = new TCanvas("plotTime_xyP", "plotTime_xyP", 3508, 2480); //optimised for A4 landscape 300dpi
  //TCanvas* plotTime_xyP = new TCanvas("plotTime_xyP", "plotTime_xyP", 1790, 1240); //optimised for A4 landscape 300dpi
  TCanvas* plotTime_xyP = new TCanvas("plotTime_xyP", "plotTime_xyP", 1440, 500);
  plotTime_xyP->Divide(5,2);
  TPad *pads[tN];

  for(int i=0;i<tN;i++)
    {
      pads[i] = (TPad*)plotTime_xyP->GetPad(i+1);
    }
  
  for(int i=0;i<tN;i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }
  plotTime_xyP->Update();

  for(int i=0;i<tN;i++)
    {
      pads[i]->cd();
      pads[i]->SetLogz();
      //optimise labels and titles
      hTime_xyP[i]->SetStats(0);
      hTime_xyP[i]->SetLabelSize(0.04,"X:Y:Z");
      //hTime_xyP[i]->SetLabelSize(0.045,"Z");
      hTime_xyP[i]->SetLabelOffset(0.009,"X");
      hTime_xyP[i]->SetLabelOffset(0.008,"Y:Z");
      hTime_xyP[i]->SetTitleSize(0.045,"X:Y:Z");
      hTime_xyP[i]->SetTitleOffset(1.,"X");
      hTime_xyP[i]->SetTitleOffset(1.1,"Y");
      hTime_xyP[i]->Draw("colz");
    }

  std::string output_name = method + "/time_primary_mirror_" + ofile_name + ".pdf";
  
  plotTime_xyP->SaveAs(output_name.c_str());
}

void plot_onAxis(TProfile2D* hist_xy[], std::string method, std::string ofile_name, int axisAngle, std::string canvasName)
{
  //on axis is array element 5
  TCanvas* plot_onAxis = new TCanvas(canvasName.c_str(),canvasName.c_str(), 400, 400);

  TEllipse *e1 = new TEllipse(0.0,0.0,100.0,100.0);//secondary circumference
  //TEllipse *e1 = new TEllipse(0.0,0.0,200.0,200.0);//primary circumference
  TEllipse *e2 = new TEllipse(0.0,0.0,60.0*mm,60*mm);// realistic hole size  //new TEllipse(0.0,0.0,18.0,18.0);//focal plane circumference
  TEllipse *e3 = new TEllipse(0.0,0.0,150.0*mm,150.0*mm);//max allowable hole size  //new TEllipse(0.0,0.0,11.8496,11.8496);//plot projected radius of led onto secondary
  TMarker *centreP = new TMarker(0.0,0.0,20);//plot centre point
  e1->SetLineColor(15);
  e1->SetLineWidth(1);
  e1->SetFillStyle(0);
  e2->SetLineColor(2);
  e2->SetLineWidth(1);
  e2->SetFillStyle(0);
  e3->SetLineWidth(1);
  e3->SetFillStyle(0);
  centreP->SetMarkerColor(1);
  centreP->SetMarkerSize(1);

  //gPad->SetRightMargin(0.15);
  //gPad->SetLogz();
  hist_xy[axisAngle]->SetTitle("");
  //optimise labels and titles
  hist_xy[axisAngle]->SetStats(0);
  hist_xy[axisAngle]->SetLabelSize(0.04,"X:Y:Z");
  //hist_xy[axisAngle]->SetLabelSize(0.045,"Z");
  hist_xy[axisAngle]->SetLabelOffset(0.009,"X");
  hist_xy[axisAngle]->SetLabelOffset(0.008,"Y:Z");
  hist_xy[axisAngle]->SetTitleSize(0.045,"X:Y:Z");
  hist_xy[axisAngle]->SetTitleOffset(1.,"X");
  hist_xy[axisAngle]->SetTitleOffset(1.1,"Y");
  hist_xy[axisAngle]->SetMarkerColor(4);
  hist_xy[axisAngle]->Draw();
  // int colN=2;
  // for(int i=9;i<10;i++)
  //   {
  //     hist_xy[i]->SetMarkerColor(colN);
  //     hist_xy[i]->Draw("same");
  //     colN++;
  //   }
  e1->Draw("same");
  e2->Draw("same");
  e3->Draw("same");
  centreP->Draw("same");
  plot_onAxis->Update();

  int countAng=-4;
  for(int i=0;i!=9;++i)
    {
      std::cout << "TH2 Total events: " << hist_xy[i]->Integral() << " angle: " << countAng << std::endl;
      countAng++;
    }

  std::string output_name = method + "/sst-gate_secondary_4.5_illumination_" + ofile_name + ".eps";
  
  //plot_onAxis->SaveAs(output_name.c_str());

  output_name = method + "/sst-gate_secondary_4.5_illumination_" + ofile_name + ".png";
  
  //plot_onAxis->SaveAs(output_name.c_str());
}

void plot_rays(const SST_GATE& gate_in, std::vector<TPolyLine3D*> poly_in)
{
  std::vector<TPolyLine3D*>::iterator poly_in_iter;

  TCanvas *cRays = new TCanvas("cRays","cRays",800,600);
  //TCanvas *cPoly = new TCanvas("cPoly","cPoly",1790, 1240);
  //cRays->SetFillColor(1);
  gate_in.GetManager()->GetTopVolume()->Draw("ogl");

  for(poly_in_iter=poly_in.begin();poly_in_iter!=poly_in.end();++poly_in_iter)
    {
      //      (*poly_in_iter)->SetLineColor(2);
      //(*poly_in_iter)->SetLineStyle(2);
      (*poly_in_iter)->Draw("same");
    }
  cRays->Update();
}

void plot_rays_multi(const SST_GATE& gate_in, std::vector< std::vector<TPolyLine3D*> > poly_in)
{
  std::vector<TPolyLine3D*>::iterator poly_iter;
  std::vector< std::vector<TPolyLine3D*> >::iterator poly_in_iter;

  TCanvas *cRaysMulti = new TCanvas("cRaysMulti","cRaysMulti",800,600);
  cRaysMulti->Divide(5,2);

  int canvasCount=1;
  
  for(poly_in_iter=poly_in.begin(); poly_in_iter!= poly_in.end();++poly_in_iter)
    {
      //gate_in.GetManager()->GetTopVolume()->Draw();
      cRaysMulti->cd(canvasCount);
      
      for(poly_iter=(*poly_in_iter).begin();poly_iter!=(*poly_in_iter).end();++poly_iter)
	{
	  (*poly_iter)->SetLineColor(2);
	  (*poly_iter)->SetLineStyle(2);
	  (*poly_iter)->Draw("same");
	}
      canvasCount++;
    }
}

void plot_shadowing(double totalRaysIn[], double sCameraIn[], double sLidIn[], double sMastsIn[], double sTrussesIn[], double sSecondaryIn[], double fov_angles[], int tN, std::string method,  std::string ofile_name)
{
  //TGraphErrors* gr_delay = new TGraphErrors(tN, fov_angles, delay_times_mean, 0, delay_times_spread);

  TGraph* gr_totalRays = new TGraph(tN, fov_angles, totalRaysIn);
  TGraph* grS_camera = new TGraph(tN, fov_angles, sCameraIn);
  TGraph* grS_lid = new TGraph(tN, fov_angles, sLidIn);
  TGraph* grS_masts = new TGraph(tN, fov_angles, sMastsIn);
  TGraph* grS_trusses = new TGraph(tN, fov_angles, sTrussesIn);
  TGraph* grS_secondary = new TGraph(tN, fov_angles, sSecondaryIn);

  TMultiGraph* mg1 = new TMultiGraph();
  TCanvas* cShadowingMulti = new TCanvas("cShadowingMulti","cShadowingMulti", 800, 600);
  gr_totalRays->SetTitle("Fractional loss of photons due to shadowing; fov (degrees); fractional loss");
  gr_totalRays->GetXaxis()->SetLabelSize(0.03);
  gr_totalRays->GetYaxis()->SetLabelSize(0.03);
  gr_totalRays->GetXaxis()->SetLabelOffset(0.008);
  gr_totalRays->GetYaxis()->SetLabelOffset(0.008);
  gr_totalRays->GetXaxis()->SetTitleSize(0.03);
  gr_totalRays->GetYaxis()->SetTitleSize(0.03);
  gr_totalRays->GetXaxis()->SetTitleOffset(1.1);
  gr_totalRays->GetYaxis()->SetTitleOffset(1.2);
  gr_totalRays->SetMaximum(1.0);
  gr_totalRays->SetMinimum(0.0);
  gr_totalRays->SetMarkerColor(4);
  gr_totalRays->SetMarkerStyle(20);
  //gr_totalRays->Draw("ALPE");
  //gr_totalRays->Draw("ALP");
  mg1->Add(gr_totalRays,"ALP");
  grS_camera->SetMarkerColor(1);
  grS_camera->SetMarkerStyle(20);
  mg1->Add(grS_camera,"LP");
  grS_lid->SetMarkerColor(2);
  grS_lid->SetMarkerStyle(20);
  mg1->Add(grS_lid,"LP");
  grS_masts->SetMarkerColor(3);
  grS_masts->SetMarkerStyle(20);
  mg1->Add(grS_masts,"LP");
  grS_trusses->SetMarkerColor(4);
  grS_trusses->SetMarkerStyle(20);
  mg1->Add(grS_trusses,"LP");
  grS_secondary->SetMarkerColor(5);
  grS_secondary->SetMarkerStyle(20);
  mg1->Add(grS_secondary,"LP");
  mg1->Draw();

  TLegend* leg = new TLegend(0.15, 0.65, 0.5, 0.85);
  leg->SetFillStyle(0);
  leg->SetTextFont(132);
  leg->AddEntry(gr_totalRays, "totals rays fired" , "lp");
  leg->AddEntry(grS_camera, "rays stopped by camera envelope" , "lp");
  leg->AddEntry(grS_lid, "rays stopped by camera lid envelope" , "lp");
  leg->AddEntry(grS_masts, "rays stopped by masts" , "lp");
  leg->AddEntry(grS_trusses, "rays stopped by trusses" , "lp");
  leg->AddEntry(grS_secondary, "rays stopped by secondary mirror" , "lp");
  leg->Draw();

  std::string output_name = method + "/shadowing_" + ofile_name + ".png";
  
  cShadowingMulti->SaveAs(output_name.c_str());

}

void plot_focal_plane(std::vector< std::vector<double> > x_arrayIn, std::vector< std::vector<double> > y_arrayIn, std::vector< std::vector<double> > z_arrayIn)
{ 
  std::vector< std::vector<double> >::iterator x_iter;
  std::vector< std::vector<double> >::iterator y_iter;
  std::vector< std::vector<double> >::iterator z_iter;
  int NN=0;

  std::cout << "x_arrayIn.size(): " << x_arrayIn.size() << std::endl;

  TGraph2D* fp_graph2D[x_arrayIn.size()];
  //TGraph* fp_graph[x_arrayIn.size()];
  y_iter=y_arrayIn.begin();
  z_iter=z_arrayIn.begin();

  for(x_iter=x_arrayIn.begin(); x_iter!=x_arrayIn.end(); ++x_iter)
    {
      fp_graph2D[NN] = new TGraph2D(Form("fp_gragh_%d", NN), Form("fp_gragh_%d", NN), (*x_iter).size(), &(*x_iter)[0], &(*y_iter)[0], &(*z_iter)[0]);
      //fp_graph[NN] = new TGraph( (*x_iter).size(), &(*x_iter)[0], &(*y_iter)[0]);
      std::cout << (*x_iter).size() << " " << (*y_iter).size() << " " << (*z_iter).size() << std::endl;
      ++y_iter;
      ++z_iter;
      NN++;
    }

  TCanvas* cFocalPlane = new TCanvas("cFocalPlane","cFocalPlane",1917, 796);
  cFocalPlane->Divide(5,2);

  TPad *pads[x_arrayIn.size()];
  
  for(int i=0;i<x_arrayIn.size();i++)
    {
      pads[i] = (TPad*)cFocalPlane->GetPad(i+1);
    }
  
  for(int i=0;i<x_arrayIn.size();i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }
  cFocalPlane->Update();


  int ccount=1;
  NN=0;

  for(int i=0; i<x_arrayIn.size(); i++)
    {
      cFocalPlane->cd(ccount);
      pads[i]->SetLogz();
      //optimise labels and titles
      // fp_graph[i]->GetXaxis()->SetLabelSize(0.04);
      // fp_graph[i]->GetYaxis()->SetLabelSize(0.04);
      // fp_graph[i]->GetZaxis()->SetLabelSize(0.04);
      // fp_graph[i]->GetXaxis()->SetLabelOffset(0.009);
      // fp_graph[i]->GetYaxis()->SetLabelOffset(0.008);
      // fp_graph[i]->GetZaxis()->SetLabelOffset(0.008);
      // fp_graph[i]->GetXaxis()->SetTitleSize(0.045);
      // fp_graph[i]->GetYaxis()->SetTitleSize(0.045);
      // fp_graph[i]->GetZaxis()->SetTitleSize(0.045);
      // fp_graph[i]->GetXaxis()->SetTitleOffset(1.);
      // fp_graph[i]->GetYaxis()->SetTitleOffset(1.1);
      // fp_graph[i]->SetMarkerStyle(6);
      // fp_graph[i]->SetMarkerSize(1);
      // fp_graph[i]->GetXaxis()->SetLimits(-1.0,1.0);
      // fp_graph[i]->SetMaximum(1.0);
      // fp_graph[i]->SetMinimum(-1.0);
      // //fp_graph[i]->SetMarkerColor();      
      // fp_graph[i]->Draw("AP");
      //TGraph 2D
      fp_graph2D[i]->GetXaxis()->SetLabelSize(0.04);
      fp_graph2D[i]->GetYaxis()->SetLabelSize(0.04);
      fp_graph2D[i]->GetZaxis()->SetLabelSize(0.04);
      fp_graph2D[i]->GetXaxis()->SetLabelOffset(0.009);
      fp_graph2D[i]->GetYaxis()->SetLabelOffset(0.008);
      fp_graph2D[i]->GetZaxis()->SetLabelOffset(0.008);
      fp_graph2D[i]->GetXaxis()->SetTitleSize(0.045);
      fp_graph2D[i]->GetYaxis()->SetTitleSize(0.045);
      fp_graph2D[i]->GetZaxis()->SetTitleSize(0.045);
      fp_graph2D[i]->GetXaxis()->SetTitleOffset(1.);
      fp_graph2D[i]->GetYaxis()->SetTitleOffset(1.1);
      fp_graph2D[i]->SetMarkerStyle(6);
      fp_graph2D[i]->SetMarkerSize(1);
      fp_graph2D[i]->GetXaxis()->SetLimits(-1.0,1.0);
      fp_graph2D[i]->SetMaximum(1.0);
      fp_graph2D[i]->SetMinimum(-1.0);
      //fp_graph2D[i]->SetMarkerColor();      
      fp_graph2D[i]->Draw("APcolz");
      ccount++;
    }


 

}

void plot_on_v_off(TProfile2D* hist_xy[], std::string method, std::string ofile_name, int on_axis, int off_axis, int mirror_index, double fieldStep, double fDegStep)
{
  std::cout << "inputs: " << fieldStep << " " << fDegStep << " " << on_axis << " " << off_axis << std::endl;
  double off_field_angle=0.;
  if(off_axis==9)
    {
      off_field_angle = ((off_axis-1)*fieldStep) + fDegStep;
      std::cout << "off field angle: " << off_field_angle << " " << fieldStep << " " << fDegStep << std::endl;
    }
  else off_field_angle = off_axis*fieldStep;

  std::stringstream get_off_field_angle;
  std::string off_field_angle_name;
  get_off_field_angle <<  off_field_angle;
  off_field_angle_name=get_off_field_angle.str();

  std::cout << "off field angle name: " << off_field_angle_name << " " << get_off_field_angle.str() << " " << off_field_angle << " " << off_axis << std::endl;

  //on axis is array element 5
  TCanvas* plot_mirror_on_v_off = new TCanvas("plot_mirror_on_v_off","plot_mirror_on_v_off", 1200, 600);
  gStyle->SetTitleFontSize(0.04);
  std::string mirror_name, plot_title; 
  TEllipse *e1;
  TEllipse *e2;
  TEllipse *e3;
  TEllipse *e4;  

  switch(mirror_index)
    {
    case 0: //primary
      e1 = new TEllipse(0.0,0.0,200.0,200.0);//set circumference
      e4 = new TEllipse(0.0,0.0,65.0,65.0);//set inner radius
      e1->SetLineColor(15);
      e1->SetLineWidth(2);
      e1->SetFillStyle(0);
      e4->SetLineColor(15);
      e4->SetLineWidth(2);
      e4->SetFillStyle(0);
      mirror_name = "primary";
      break;
    case 1: //secondary
      e1 = new TEllipse(0.0,0.0,100.0,100.0);       // set circumference
      e2 = new TEllipse(0.0,0.0,60.0*mm,60*mm);     // realistic hole siz
      e3 = new TEllipse(0.0,0.0,150.0*mm,150.0*mm); // max allowable hole size
      e1->SetLineColor(15);
      e1->SetLineWidth(1);
      e1->SetFillStyle(0);
      e2->SetLineColor(2);
      e2->SetLineWidth(1);
      e2->SetFillStyle(0);
      e3->SetLineWidth(1);
      e3->SetFillStyle(0);
      mirror_name = "secondary";
      break;
    }
  //TEllipse *e2 = new TEllipse(0.0,0.0,18.0,18.0);//focal plane circumference
  //TEllipse *e3 = new TEllipse(0.0,0.0,11.8496,11.8496);//plot projected radius of led onto secondary
  TMarker *centreP = new TMarker(0.0,0.0,20);//plot centre point
  
  centreP->SetMarkerColor(1);
  centreP->SetMarkerSize(1);

  //gPad->SetRightMargin(0.15);
  //gPad->SetLogz();

  plot_mirror_on_v_off->Divide(2,1);
  plot_mirror_on_v_off->cd(1);
  plot_title = "#theta_{field angle} = 0 (degrees)";
  //plot_title = "on-axis (#theta = 0 degrees) observation for " + mirror_name + " mirror.";
  hist_xy[on_axis]->SetTitle(plot_title.c_str());
  //optimise labels and titles
  hist_xy[on_axis]->SetStats(0);
  hist_xy[on_axis]->SetLabelSize(0.035,"X:Y:Z");//0.028 //0.035 for publication
  //hist_xy[on_axis]->SetLabelSize(0.045,"Z");
  hist_xy[on_axis]->SetLabelOffset(0.009,"X");
  hist_xy[on_axis]->SetLabelOffset(0.008,"Y:Z");
  hist_xy[on_axis]->SetTitleSize(0.04,"X:Y:Z");//0.03 //0.04 for publication
  hist_xy[on_axis]->SetTitleOffset(1.,"X");
  hist_xy[on_axis]->SetTitleOffset(1.3,"Y");
  hist_xy[on_axis]->SetMarkerColor(4);
  hist_xy[on_axis]->Draw();
  // int colN=2;
  // for(int i=9;i<10;i++)
  //   {
  //     hist_xy[i]->SetMarkerColor(colN);
  //     hist_xy[i]->Draw("same");
  //     colN++;
  //   }
  switch(mirror_index)
    {
    case 0:
      e1->Draw("same");
      e4->Draw("same");
      centreP->Draw("same");
      break;
    case 1:
      e1->Draw("same");
      //e2->Draw("same");
      //e3->Draw("same");
      centreP->Draw("same");
      break;
    }
  

  plot_mirror_on_v_off->cd(2);
  plot_title = "#theta_{field angle} = " + off_field_angle_name +  " (degrees)";
  //plot_title = "off-axis (#theta = 4 degrees) observation for " + mirror_name + " mirror.";
  hist_xy[off_axis]->SetTitle(plot_title.c_str());
  //optimise labels and titles
  hist_xy[off_axis]->SetStats(0);
  hist_xy[off_axis]->SetLabelSize(0.035,"X:Y:Z");//0.028 //0.035 for publication
  //hist_xy[off_axis]->SetLabelSize(0.045,"Z");
  hist_xy[off_axis]->SetLabelOffset(0.009,"X");
  hist_xy[off_axis]->SetLabelOffset(0.008,"Y:Z");
  hist_xy[off_axis]->SetTitleSize(0.04,"X:Y:Z");//0.03 //0.04 for publication
  hist_xy[off_axis]->SetTitleOffset(1.,"X");
  hist_xy[off_axis]->SetTitleOffset(1.3,"Y");
  hist_xy[off_axis]->SetMarkerColor(4);
  hist_xy[off_axis]->Draw();
  // int colN=2;
  // for(int i=9;i<10;i++)
  //   {
  //     hist_xy[i]->SetMarkerColor(colN);
  //     hist_xy[i]->Draw("same");
  //     colN++;
  //   }
    
  switch(mirror_index)
    {
    case 0:
      e1->Draw("same");
      e4->Draw("same");
      centreP->Draw("same");
      break;
    case 1:
      e1->Draw("same");
      //e2->Draw("same");
      //e3->Draw("same");
      centreP->Draw("same");
      break;
    }
  
  plot_mirror_on_v_off->Update();

  int countAng=-4;
  for(int i=0;i!=9;++i)
    {
      std::cout << "TH2 Total events: " << hist_xy[i]->Integral() << " angle: " << countAng << std::endl;
      countAng++;
    }

  std::string output_name = method + "/sst-gate_" + mirror_name + "_on_v_off_" + off_field_angle_name + "_degrees_" +  ofile_name + ".eps";
  
  plot_mirror_on_v_off->SaveAs(output_name.c_str());

  output_name = method + "/sst-gate_" + mirror_name + "_on_v_off_" + off_field_angle_name + "_degrees_"  + ofile_name + ".png";
  
  plot_mirror_on_v_off->SaveAs(output_name.c_str());
}



#endif


/*
void plot_time(TProfile* hTime_xFP[], TProfile* hTime_yFP[], TProfile2D* hTime_xyFP_S[], TProfile2D* hTime_xyS[], TProfile2D* hTime_xyFP_Total[], TProfile2D* hTime_xyFP_Delay[], int tN, std::string method, std::string ofile_name, double xD[], double yD[], double fRpMax, double kDegStep)
{

  TCanvas* plotSpot_cm = new TCanvas("plotSpot_cm", "plotSpot_cm", 3508, 2480); //optimised for A4 landscape 300dpi
  plotSpot_cm->Divide(5, 2);
  TPad *pads[tN];
  std::stringstream get_padNum;
  std::string padNum;
  TPaveStats *mystats[tN];
  TPaveStats *newStats;

  for(int i=0;i<tN;i++)
    {
      pads[i] = (TPad*)plotSpot_cm->GetPad(i+1);
    }

  for(int i=0;i<tN;i++)
    {
      pads[i]->SetRightMargin(0.15);
      pads[i]->Draw();
    }

  plotSpot_cm->Update();

  double angle=0.0; //reset angle

  for(int i=0;i<tN;i++)
    {
      //change pad
      pads[i]->cd();
           
      //To draw zero value bins as lowest colour uncomment this line.
      //spotH[i]->SetMinimum(-1); 
      //std::cout << "mean: " << spotH[i]->GetMean(1) << std::endl;
      xD[i] = spotH_cm[i]->GetMean(1);
      yD[i] = angle;

      //optimise labels and titles
      spotH_cm[i]->SetLabelSize(0.0425,"X:Y");
      spotH_cm[i]->SetLabelSize(0.045,"Z");
      spotH_cm[i]->SetLabelOffset(0.009,"X");
      spotH_cm[i]->SetLabelOffset(0.008,"Y:Z");
      spotH_cm[i]->SetTitleSize(0.045,"X:Y");
      spotH_cm[i]->SetTitleOffset(1.,"X");
      spotH_cm[i]->SetTitleOffset(1.1,"Y");

      //optimise axis range
      spotH_cm[i]->SetAxisRange(spotH_cm[i]->GetMean(1)-0.6,spotH_cm[i]->GetMean(1)+0.6,"X;Y");
      //spotH_cm[i]->SetAxisRange(spotH_cm[i]->GetMean(1)-0.6,spotH_cm[i]->GetMean(1)+0.6,"X");
      //spotH_cm[i]->SetAxisRange(spotH_cm[i]->GetMean(2)-0.6,spotH_cm[i]->GetMean(2)+0.6,"Y");
      
      //draw spot
      spotH_cm[i]->Draw("colz");

      //step through FOV angles
      angle+=kDegStep;
    }

  plotSpot_cm->Modified();
  plotSpot_cm->Update();

  //reposition stats boxes
  for(int i=0;i<tN;i++)
    {
      pads[i]->cd();
      get_padNum << i+1;
      padNum = "stats_" + get_padNum.str();     
      newStats = (TPaveStats*)pads[i]->GetPrimitive("stats");
      newStats->SetName(padNum.c_str());  
      newStats->SetX1NDC(0.213861);
      newStats->SetY1NDC(0.689887);
      newStats->SetX2NDC(0.734114);
      newStats->SetY2NDC(0.947591);
      get_padNum.str("");
    }

  //close files and save plots
  output_name = method + "/spot_" + ofile_name + ".pdf";
  plotSpot_cm->SaveAs(output_name.c_str());
  output_name = method + "/spot_" + ofile_name + ".jpg";
  plotSpot_cm->SaveAs(output_name.c_str());

}
*/
