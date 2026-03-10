#ifndef UTILITIES_H
#define UTILITIES_H

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <fstream>
#include "../Armory/AnalysisLibrary.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

const int N = 6;

const std::vector<double> peakPositions = {0.0, 1.982, 3.552, 3.630, 3.920, 5.255};

void readFitParameters(const TString &fileName, std::vector<std::vector<double>> &params) {
  std::ifstream inFile(fileName);
  if (!inFile.is_open()) {
    std::cerr << "Failed to open fit parameters file!" << std::endl;
    return;
  }

  std::string line;
  while (std::getline(inFile, line)) {
    std::istringstream iss(line);
    std::string temp;
    int detectorIndex;

    iss >> temp >> detectorIndex >> temp;

    std::vector<double> fitParams;
    double param;
    while (iss >> param) {
      fitParams.push_back(param);
    }

    if (detectorIndex >= params.size()) {
      params.resize(detectorIndex + 1);
    }
    params[detectorIndex] = fitParams;
  }
  inFile.close();
}

double Peak(double *dim, double *par){

  double  x       = dim[0];

  double  area    = par[0];
  double  cent    = par[1];
  double  sigma   = par[2];

  return area/(sigma*TMath::Sqrt(2*TMath::Pi())) * TMath::Gaus(x,cent,sigma);

}

double FitNPeaks(double *dim, double *par){

  double  x       = dim[0];

  double  val     = 0;

  double  sigma   = par[0];
  double  p0      = par[1];
  double  p1      = par[2];

  double  *PeakPar    = new double[3];

  for(int i=0;i<N;i++){
    PeakPar[0]  = par[3+i*2]; // area of the peak
    if(i==0)
      PeakPar[1] = par[4]; // Center (mean) of the peak
    else
      PeakPar[1] = par[4+i*2] + par[4];
    PeakPar[2]  = sigma;
    val += Peak(dim,PeakPar);
  }

  val += p0 + p1*x;
    
  return val;
}

double calculateArea(double amplitude, double width) {
  return amplitude * TMath::Sqrt(2 * TMath::Pi()) * width;
}

void fitSpectra(TH1F* hist, int detectorId) {

  TF1* fitFunc = new TF1("fitFunc", FitNPeaks, -1, 7, 3 + N * 2); // 6 peaks and 3 background parameters
    
  // Set the initial parameters for the background (linear)
  fitFunc->SetParameter(0, 0.1);
  fitFunc->SetParLimits(0, 0.05, 0.15);
  fitFunc->SetParameter(1, 0);
  fitFunc->SetParameter(2, 0);
    
  // Set parameter limits for the peaks and background
  for (int i = 0; i < N; i++) {
    fitFunc->SetParLimits(3 + i * 2, 0, 100000);
    fitFunc->SetParLimits(4 + i * 2, peakPositions[i] - 0.2, peakPositions[i] + 0.2);
    fitFunc->SetParName(3+2*i,Form("Area%i",i+1));        
    fitFunc->SetParName(4+2*i,Form("Cent%i",i+1));
  }

  if (detectorId == 0) {
    fitFunc->FixParameter(4, peakPositions[0]);
  } else if (detectorId == 1) {
    fitFunc->FixParameter(4, peakPositions[0]);
    fitFunc->FixParameter(6, peakPositions[1]);
    fitFunc->FixParameter(8, peakPositions[2]);
  } else if (detectorId == 2) {
    fitFunc->FixParameter(4, peakPositions[0]);
    fitFunc->FixParameter(6, peakPositions[1]);
    fitFunc->FixParameter(8, peakPositions[2]);
    fitFunc->FixParameter(10, peakPositions[3]);
    fitFunc->FixParameter(12, peakPositions[4]);
  } else if (detectorId == 3) {
    fitFunc->FixParameter(4, peakPositions[0]);
    fitFunc->FixParameter(6, peakPositions[1]);
    fitFunc->FixParameter(8, peakPositions[2]);
    fitFunc->FixParameter(10, peakPositions[3]);
    fitFunc->FixParameter(12, peakPositions[4]);
    fitFunc->FixParameter(14, peakPositions[5]);
  } else {
    for (int i = 0; i < N; i++) {
      fitFunc->FixParameter(4 + i * 2, peakPositions[i]);
    }
  }

  hist->Fit(fitFunc, "R");

  std::cout << "Fitting Results for Detector " << detectorId << ":\n";
  for (int i = 0; i < N; i++) {
    double peakPosition = fitFunc->GetParameter(i * 2);
    double amplitude = fitFunc->GetParameter(i * 2 + 1);
    double width = fitFunc->GetParameter(i * 2 + 2);
            
    double area = calculateArea(amplitude, width);
            
    std::cout << "Peak " << i << " - Position: " << peakPosition << " MeV, Area: " << area << std::endl;
  }
    
  TCanvas* c1 = new TCanvas("c1", "Fitting Results", 800, 600);
  hist->Draw();
  fitFunc->Draw("same");
    
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hist, "Data", "l");
  legend->AddEntry(fitFunc, "Fit", "l");
  legend->Draw();
    
  c1->Update();
  c1->SaveAs(Form("17O_analysis_ExCorr/Fit_det%d.png", detectorId));
  c1->SaveAs(Form("17O_analysis_ExCorr/Fit_det%d.root", detectorId));
}


#endif
