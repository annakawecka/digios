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

TString rdtCutFile = "rdtCuts_17O.root";

TObjArray * cutList;
Bool_t isCutFileOpen;
int numCut;

TCutG* cutG;

bool rdtgate = false;
bool xgate = false;
bool cointimegate = false;

const int nDet = 24;

int n_coin = 0;
int n_coin_rdt = 0;

int recoil_n = 0;

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

void analysis_17O(){
  //================================= coinTime fit parameters

  std::vector<std::vector<double>> fitParams;
  readFitParameters("coinTimeCorr_17O/coinTime_x_fit_parameters.txt", fitParams);

  if (fitParams.empty()) {
    std::cerr << "No parameters loaded from file!" << std::endl;
    return;
  }

  //================================= rdt cut

  TFile * fCut = new TFile(rdtCutFile);
  isCutFileOpen = fCut->IsOpen();
  if(!isCutFileOpen) {
    printf( "Failed to open rdt-cutfile : %s\n" , rdtCutFile.Data());
    rdtCutFile = "";
  }
  numCut = 0 ;

  if( isCutFileOpen ){
    cutList = (TObjArray *) fCut->FindObjectAny("cutList");
    numCut = cutList->GetEntries();
    printf("=========== found %d cutG in %s \n", numCut, fCut->GetName());

    for(int i = 0; i < numCut ; i++){
      printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
	     cutList->At(i)->GetName(),
	     ((TCutG*)cutList->At(i))->GetVarX(),
	     ((TCutG*)cutList->At(i))->GetVarY(),
	     ((TCutG*)cutList->At(i))->GetN());
    }
  }

  //================================= reading the tree

  TChain *chain = new TChain("tree");
  chain->Add("trace_run055-066.root");

  std::vector<TH1F*> correctedCoinTime;
  std::vector<TH1F*> correctedCoinTimeXgate;
  std::vector<TH1F*> correctedCoinTimeRDTCoin;
  std::vector<TH1F*> correctedCoinTimeXgateRDTCoin;
  std::vector<TH2F*> recoilDEE;

  TH1F* x_rdt_coinTime_gatedEx = new TH1F("x_rdt_coinTime_gatedEx", "Ex gated on x, recoils and coinTime", 200, -2, 12);
  TH1F* rdt_coinTime_gatedEx = new TH1F("rdt_coinTime_gatedEx", "Ex gated on recoils and coinTime", 200, -2, 12);
  TH1F* coinTime_gatedEx = new TH1F("coinTime_gatedEx", "Ex gated on coinTime", 200, -2, 12);
  TH1F* rdt_gatedEx = new TH1F("rdt_gatedEx", "Ex gated on recoils", 200, -2, 12);

  for (int i = 0; i < nDet; ++i) {
    TString histName;
    histName.Form("corrected_coinTime_det%d", i);
    correctedCoinTime.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("corrected_coinTime_det%d_x_gate", i);
    correctedCoinTimeXgate.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("corrected_coinTime_det%d_rdt_coincidence", i);
    correctedCoinTimeRDTCoin.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("corrected_coinTime_det%d_x_gate_rdt_coincidence", i);
    correctedCoinTimeXgateRDTCoin.push_back(new TH1F(histName, histName, 400, -100, 200));
  }

  for (int i = 0; i < 4; ++i) {
    TString histName;
    histName.Form("recoilDEE_%d_%d", i, i+1);
    recoilDEE.push_back(new TH2F(histName, histName, 2096, 0, 10000, 2096, 0, 10000));
  }

  // Set the branch addresses
  Float_t rdt[8], Ex, coinTime, x[24];
  Int_t detID;
  Float_t coinTimeCorr;
  chain->SetBranchAddress("rdt", rdt);
  chain->SetBranchAddress("Ex", &Ex);
  chain->SetBranchAddress("coinTime", &coinTime);
  chain->SetBranchAddress("x", x);
  chain->SetBranchAddress("detID", &detID);

  // Loop over the events and apply the cut
  Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    rdtgate = false;
    xgate = false;
    cointimegate = false;

    chain->GetEntry(i);

    if( isCutFileOpen ){
      for(int i = 0 ; i < numCut ; i++ ){
	cutG = (TCutG *)cutList->At(i) ;
	if(cutG->IsInside(rdt[2*i],rdt[2*i+1])) {
	  rdtgate = true;
	  recoil_n = 2*i;
	  break; /// only one is enough
	}
      }
    }

    if ( x[detID] < 0.95 && x[detID] > -0.95)
      xgate = true;

    if ( fitParams[detID].size() == 3 ) {
      coinTimeCorr = coinTime - ( fitParams[detID][0] + fitParams[detID][1] * x[detID] + fitParams[detID][2] * TMath::Power(x[detID], 2) );
    } else if ( fitParams[detID].size() == 4 ) {
      coinTimeCorr = coinTime - ( fitParams[detID][0] + fitParams[detID][1] * x[detID] + fitParams[detID][2] * TMath::Power(x[detID], 2) + fitParams[detID][3] * TMath::Power(x[detID], 3));
    } else if ( fitParams[detID].size() == 5 ) {
      coinTimeCorr = coinTime - ( fitParams[detID][0] + fitParams[detID][1] * x[detID] + fitParams[detID][2] * TMath::Power(x[detID], 2) + fitParams[detID][3] * TMath::Power(x[detID], 3) + fitParams[detID][4] * TMath::Power(x[detID], 4));
    }
    else {
      std::cerr << "Unexpected number of fit parameters for detector " << detID << " (" << fitParams[detID].size() << "). Skipping." << std::endl;
      continue;
    }

    if (coinTimeCorr > -20 && coinTimeCorr < 15) {
      cointimegate = true;

      coinTime_gatedEx->Fill(Ex);
    }

    correctedCoinTime[detID]->Fill(coinTimeCorr);
    n_coin++;

    if (xgate)
      correctedCoinTimeXgate[detID]->Fill(coinTimeCorr);

    if (rdtgate) {
      correctedCoinTimeRDTCoin[detID]->Fill(coinTimeCorr);
      recoilDEE[recoil_n/2]->Fill(rdt[recoil_n], rdt[recoil_n+1]);
      n_coin_rdt++;

      rdt_gatedEx->Fill(Ex);
    }

    if (rdtgate && cointimegate)
      rdt_coinTime_gatedEx->Fill(Ex);

    if (xgate && rdtgate) {
      correctedCoinTimeXgateRDTCoin[detID]->Fill(coinTimeCorr);

      if (cointimegate)
	x_rdt_coinTime_gatedEx->Fill(Ex);
    }

  }

  TCanvas *cAllDetectors = new TCanvas("cAllDetectors", "Corrected CoinTime for All Detectors", 1200, 800);
  cAllDetectors->Divide(6, 4);

  for (int i = 0; i < nDet; ++i) {
    cAllDetectors->cd(i + 1);
    correctedCoinTime[i]->Draw();
  }

  cAllDetectors->SaveAs("17O_analysis/Corrected_CoinTime_AllDetectors.png");

  TCanvas *cAllDetectorsXgate = new TCanvas("cAllDetectorsXgate", "Corrected CoinTime for All Detectors with x gate", 1200, 800);
  cAllDetectorsXgate->Divide(6, 4);

  for (int i = 0; i < nDet; ++i) {
    cAllDetectorsXgate->cd(i + 1);
    correctedCoinTimeXgate[i]->Draw();
  }

  cAllDetectorsXgate->SaveAs("17O_analysis/Corrected_CoinTime_AllDetectors_Xgate.png");

  TCanvas *cAllDetectorsRDT = new TCanvas("cAllDetectorsRDT", "Corrected CoinTime for All Detectors (RDT coin)", 1200, 800);
  cAllDetectorsRDT->Divide(6, 4);

  for (int i = 0; i < nDet; ++i) {
    cAllDetectorsRDT->cd(i + 1);
    correctedCoinTimeRDTCoin[i]->Draw();
  }

  cAllDetectorsRDT->SaveAs("17O_analysis/Corrected_CoinTime_AllDetectors_RDTCoin.png");

  TCanvas *cAllDetectorsXgateRDT = new TCanvas("cAllDetectorsXgateRDT", "Corrected CoinTime for All Detectors (RDT coin && x gate)", 1200, 800);
  cAllDetectorsXgateRDT->Divide(6, 4);

  for (int i = 0; i < nDet; ++i) {
    cAllDetectorsXgateRDT->cd(i + 1);
    correctedCoinTimeXgateRDTCoin[i]->Draw();
  }

  cAllDetectorsRDT->SaveAs("17O_analysis/Corrected_CoinTime_AllDetectors_XgateRDTCoin.png");

  std::cout << "coin: " << n_coin << " coin_rdt: " << n_coin_rdt << std::endl;

  TCanvas *cRDT = new TCanvas("cRDT", "RDT dE-E gated", 1200, 300);
  cRDT->Divide(4);

  for (int i = 0; i < 4; ++i) {
    cRDT->cd(i + 1);
    recoilDEE[i]->Draw();
  }

  cRDT->SaveAs("17O_analysis/RDT_DEE_gated.png");

  TH1F *combinedHist = (TH1F*)correctedCoinTimeXgateRDTCoin[0]->Clone("combinedHist");
  combinedHist->Reset();
  for (int i = 0; i < nDet; ++i) {
    combinedHist->Add(correctedCoinTimeXgateRDTCoin[i]);
  }

  TCanvas *cCombined = new TCanvas("cCombined", "Combined Corrected CoinTime with recoil && x gate", 800, 600);
  combinedHist->Draw();
  cCombined->SaveAs("17O_analysis/Combined_Corrected_CoinTime_XgateRDTCoin.png");

  TCanvas *cGatedEx = new TCanvas("cGatedEx", "Ex gated on x, recoils and coinTime", 800, 600);
  x_rdt_coinTime_gatedEx->Draw();
  cGatedEx->SaveAs("17O_analysis/Ex_x_recoil_coinTime_gated.png");

  TCanvas *cRDTCoinTimeGatedEx = new TCanvas("cRDTCoinTimeGatedEx", "Ex gated on recoils and coinTime", 800, 600);
  rdt_coinTime_gatedEx->Draw();
  cRDTCoinTimeGatedEx->SaveAs("17O_analysis/Ex_recoil_coinTime_gated.png");

  TCanvas *cCoinTimeGatedEx = new TCanvas("cCoinTimeGatedEx", "Ex gated on coinTime", 800, 600);
  coinTime_gatedEx->Draw();
  cCoinTimeGatedEx->SaveAs("17O_analysis/Ex_coinTime_gated.png");

  TCanvas *cRDTGatedEx = new TCanvas("cRDTGatedEx", "Ex gated on rdt", 800, 600);
  rdt_gatedEx->Draw();
  cRDTGatedEx->SaveAs("17O_analysis/Ex_rdt_gated.png");
}
