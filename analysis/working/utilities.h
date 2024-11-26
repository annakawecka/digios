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

#define N 5

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
    PeakPar[0]  = par[3+i*2];
    if(i==0)
      PeakPar[1] = par[4];
    else
      PeakPar[1] = par[4+i*2] + par[4];
    PeakPar[2]  = sigma;
    val += Peak(dim,PeakPar);
  }

  val += p0 + p1*x;
    
  return val;
}


#endif
