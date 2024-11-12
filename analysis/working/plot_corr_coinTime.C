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

void plot_corr_coinTime() {
  std::vector<std::vector<double>> fitParams;
  //readFitParameters("coinTimeCorr_17O/coinTime_x_fit_parameters.txt", fitParams);
  readFitParameters("coinTimeCorr_17F/coinTime_x_fit_parameters.txt", fitParams);

  if (fitParams.empty()) {
    std::cerr << "No parameters loaded from file!" << std::endl;
    return;
  }

  TChain * chain = new TChain("tree");
  //chain->Add("trace_run055-066.root"); // 17O
  chain->Add("trace_run022-052.root"); // 17F

  const int nDet = fitParams.size();
  std::vector<TH1F*> correctedHistVec;

  for (int i = 0; i < nDet; ++i) {
    TString histName;
    histName.Form("corrected_coinTime_det%d", i);
    correctedHistVec.push_back(new TH1F(histName, histName, 400, -100, 200));
  }

  for (int i = 0; i < nDet; ++i) {
    if (fitParams[i].empty()) {
      std::cerr << "No fit parameters for detector " << i << ". Skipping." << std::endl;
      continue;
    }

    TString expression;
    if (fitParams[i].size() == 3) {
      expression.Form("coinTime - (%f + %f * x[%d] + %f * x[%d]**2) >> corrected_coinTime_det%d",
		      fitParams[i][0],
		      fitParams[i][1],
		      i,
		      fitParams[i][2],
		      i,
		      i);
    } else if (fitParams[i].size() == 4) {
      expression.Form("coinTime - (%f + %f * x[%d] + %f * x[%d]**2 + %f * x[%d]**3) >> corrected_coinTime_det%d",
		      fitParams[i][0],
		      fitParams[i][1],
		      i,
		      fitParams[i][2],
		      i,
		      fitParams[i][3],
		      i,
		      i);
    } else if (fitParams[i].size() == 5) {
      expression.Form("coinTime - (%f + %f * x[%d] + %f * x[%d]**2 + %f * x[%d]**3 + %f * x[%d]**4) >> corrected_coinTime_det%d",
		      fitParams[i][0],
		      fitParams[i][1],
		      i,
		      fitParams[i][2],
		      i,
		      fitParams[i][3],
		      i,
		      fitParams[i][4],
		      i,
		      i);
    } else {
      std::cerr << "Unexpected number of fit parameters for detector " << i << ". Skipping." << std::endl;
      continue;
    }

    chain->Draw(expression, Form("x[%d] > -0.95 && x[%d] < 0.95", i, i), "goff");
  }

  TCanvas *cAllDetectors = new TCanvas("cAllDetectors", "Corrected CoinTime for All Detectors", 1200, 800);
  cAllDetectors->Divide(6, 4);

  for (int i = 0; i < nDet; ++i) {
    cAllDetectors->cd(i + 1);
    correctedHistVec[i]->Draw();
  }

  cAllDetectors->SaveAs("Corrected_CoinTime_AllDetectors.png");

  TH1F *combinedHist = (TH1F*)correctedHistVec[0]->Clone("combinedHist");
  combinedHist->Reset();
  for (int i = 0; i < nDet; ++i) {
    combinedHist->Add(correctedHistVec[i]);
  }

  TCanvas *cCombined = new TCanvas("cCombined", "Combined Corrected CoinTime", 800, 600);
  combinedHist->Draw();
  cCombined->SaveAs("Combined_Corrected_CoinTime.png");

  for (auto hist : correctedHistVec) {
    delete hist;
  }
}
