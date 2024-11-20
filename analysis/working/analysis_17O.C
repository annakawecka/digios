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

bool plothist = false;

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
  TH1F* Ex_nogates = new TH1F("Ex_nogates", "Ex, no gates", 200, -2, 12);

  TH1F* Ex_d0 = new TH1F("Ex_d0", "Ex, det 0, 6, 12, 18", 200, -2, 12);
  TH1F* Ex_d1 = new TH1F("Ex_d1", "Ex, det 1, 7, 13, 19", 200, -2, 12);
  TH1F* Ex_d2 = new TH1F("Ex_d2", "Ex, det 2, 8, 14, 20", 200, -2, 12);
  TH1F* Ex_d3 = new TH1F("Ex_d3", "Ex, det 3, 9, 15, 21", 200, -2, 12);
  TH1F* Ex_d4 = new TH1F("Ex_d4", "Ex, det 4, 10, 16, 22", 200, -2, 12);
  TH1F* Ex_d5 = new TH1F("Ex_d5", "Ex, det 5, 11, 17, 23", 200, -2, 12);

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

      if (cointimegate) {
	x_rdt_coinTime_gatedEx->Fill(Ex);

	if (detID % 6 == 0)
	  Ex_d0->Fill(Ex);
	if (detID % 6 == 1)
	  Ex_d1->Fill(Ex);
	if (detID % 6 == 2)
	  Ex_d2->Fill(Ex);
	if (detID % 6 == 3)
	  Ex_d3->Fill(Ex);
	if (detID % 6 == 4)
	  Ex_d4->Fill(Ex);
	if (detID % 6 == 5)
	  Ex_d5->Fill(Ex);
      }
    }

    Ex_nogates->Fill(Ex);

  }

  //================================= fitting 17O

  TCanvas *cExdet = new TCanvas("cExdet", "Ex for different detector rings", 1000, 800);
  cExdet->Divide(3, 2);

  cExdet->cd(1);
  Ex_d0->Draw();

  cExdet->cd(2);
  Ex_d1->Draw();

  cExdet->cd(3);
  Ex_d2->Draw();

  cExdet->cd(4);
  Ex_d3->Draw();

  cExdet->cd(5);
  Ex_d4->Draw();

  cExdet->cd(6);
  Ex_d5->Draw();

  TF1 *f1 = new TF1("Fit",FitNPeaks,-1,7,2*N+3);

  f1->SetParameter(0,0.1);
  f1->SetParLimits(0,0.05,0.15);
  f1->SetParameter(1, 0);
  f1->SetParameter(2, 0); 

  f1->SetParName(0,"sigma");
  f1->SetParName(1,"p0");
  f1->SetParName(2,"p1");
  for(int i=0;i<N;i++){
    f1->SetParName(3+2*i,Form("Area%i",i+1));        
    f1->SetParName(4+2*i,Form("Cent%i",i+1));
  }

  // Peak 1
  f1->SetParameter(3,10);
  f1->SetParLimits(3,0,1000);
  f1->SetParameter(4,0);
  //f1->SetParLimits(4,-0.5,0.5);
  // Peak 2
  f1->SetParameter(5,10);
  f1->SetParLimits(5,0,1000);
  f1->SetParameter(6,2.0);
  //f1->SetParLimits(6,1.8,2.5);
  //f1->FixParameter(6,0.9376);
  // Peak 3
  f1->SetParameter(7,10);
  f1->SetParLimits(7,0,1000);
  f1->SetParameter(8,3.7);
  //f1->SetParLimits(8,2.9,3.9);
  //f1->FixParameter(8,1.0416);
  // Peak 4
  //f1->SetParameter(9,10);
  //f1->FixParameter(9,0);
  //f1->SetParLimits(9,0,1000);
  //f1->SetParameter(10,3.8);
  //f1->SetParLimits(10,3.1,3.9);
  //f1->FixParameter(10,3.9);
  // Peak 5
  f1->SetParameter(9,10);
  f1->SetParLimits(9,0,1000);
  f1->SetParameter(10,4.1);
  //f1->SetParLimits(10,3.9,4.3);
  // Peak 6
  f1->SetParameter(11,10);
  f1->SetParLimits(11,0,1000);
  f1->SetParameter(12,5.5);
  //f1->SetParLimits(12,5.3,6.0);
  //f1->FixParameter(12,1.1214);

  
  TCanvas *cExd5 = new TCanvas("cExd5", "Ex, det 5...", 800, 600);
  f1->SetNpx(1000);

  Ex_d5->Fit(f1, "RB0");

  Ex_d5->SetBinErrorOption(TH1::kPoisson);
  
  Ex_d5->Draw("e1");
  Ex_d5->SetMarkerStyle(20);
  Ex_d5->SetMarkerColor(kRed);
  Ex_d5->SetMarkerSize(1);
  f1->SetLineWidth(4);
  f1->Draw("same");

  printf("Detectors 5, 11, 17, 23 \n");

  TF1 *indpeak[N];
  for(int i=0;i<N;i++){
    indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
    indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
    if(i==0)
      indpeak[i]->SetParameter(1,f1->GetParameter(4));
    else
      indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
    indpeak[i]->SetParameter(2,f1->GetParameter(0));
    indpeak[i]->SetLineStyle(7);
    indpeak[i]->Draw("same");
    if (i==0)
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
    else
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
  }

  TCanvas *cExd4 = new TCanvas("cExd4", "Ex, det 4...", 800, 600);
  f1->SetNpx(1000);

  Ex_d4->Fit(f1, "RB0");

  Ex_d4->SetBinErrorOption(TH1::kPoisson);

  Ex_d4->Draw("e1");
  Ex_d4->SetMarkerStyle(20);
  Ex_d4->SetMarkerColor(kRed);
  Ex_d4->SetMarkerSize(1);
  f1->SetLineWidth(4);
  f1->Draw("same");

  printf("Detectors 4, 10, 16, 22 \n");

  for(int i=0;i<N;i++){
    indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
    indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
    if(i==0)
      indpeak[i]->SetParameter(1,f1->GetParameter(4));
    else
      indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
    indpeak[i]->SetParameter(2,f1->GetParameter(0));
    indpeak[i]->SetLineStyle(7);
    indpeak[i]->Draw("same");
    if (i==0)
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
    else
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
  }

  TCanvas *cExd3 = new TCanvas("cExd3", "Ex, det 3...", 800, 600);
  f1->SetNpx(1000);

  Ex_d3->Fit(f1, "RB0");

  Ex_d3->SetBinErrorOption(TH1::kPoisson);

  Ex_d3->Draw("e1");
  Ex_d3->SetMarkerStyle(20);
  Ex_d3->SetMarkerColor(kRed);
  Ex_d3->SetMarkerSize(1);
  f1->SetLineWidth(4);
  f1->Draw("same");

  printf("Detectors 3, 9, 15, 21 \n");

  for(int i=0;i<N;i++){
    indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
    indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
    if(i==0)
      indpeak[i]->SetParameter(1,f1->GetParameter(4));
    else
      indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
    indpeak[i]->SetParameter(2,f1->GetParameter(0));
    indpeak[i]->SetLineStyle(7);
    indpeak[i]->Draw("same");
    if (i==0)
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
    else
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
  }

  TCanvas *cExd2 = new TCanvas("cExd2", "Ex, det 2...", 800, 600);
  f1->SetNpx(1000);

  Ex_d2->Fit(f1, "RB0");

  Ex_d2->SetBinErrorOption(TH1::kPoisson);

  Ex_d2->Draw("e1");
  Ex_d2->SetMarkerStyle(20);
  Ex_d2->SetMarkerColor(kRed);
  Ex_d2->SetMarkerSize(1);
  f1->SetLineWidth(4);
  f1->Draw("same");

  printf("Detectors 2, 8, 14, 20 \n");

  for(int i=0;i<N;i++){
    indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
    indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
    if(i==0)
      indpeak[i]->SetParameter(1,f1->GetParameter(4));
    else
      indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
    indpeak[i]->SetParameter(2,f1->GetParameter(0));
    indpeak[i]->SetLineStyle(7);
    indpeak[i]->Draw("same");
    if (i==0)
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
    else
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
  }

  TCanvas *cExd1 = new TCanvas("cExd1", "Ex, det 1...", 800, 600);
  f1->SetNpx(1000);

  Ex_d1->Fit(f1, "RB0");

  Ex_d1->SetBinErrorOption(TH1::kPoisson);

  Ex_d1->Draw("e1");
  Ex_d1->SetMarkerStyle(20);
  Ex_d1->SetMarkerColor(kRed);
  Ex_d1->SetMarkerSize(1);
  f1->SetLineWidth(4);
  f1->Draw("same");

  printf("Detectors 1, 7, 13, 19 \n");

  for(int i=0;i<N;i++){
    indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
    indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
    if(i==0)
      indpeak[i]->SetParameter(1,f1->GetParameter(4));
    else
      indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
    indpeak[i]->SetParameter(2,f1->GetParameter(0));
    indpeak[i]->SetLineStyle(7);
    indpeak[i]->Draw("same");
    if (i==0)
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
    else
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
  }

  TCanvas *cExd0 = new TCanvas("cExd0", "Ex, det 0...", 800, 600);
  f1->SetNpx(1000);

  Ex_d0->Fit(f1, "RB0");

  Ex_d0->SetBinErrorOption(TH1::kPoisson);

  Ex_d0->Draw("e1");
  Ex_d0->SetMarkerStyle(20);
  Ex_d0->SetMarkerColor(kRed);
  Ex_d0->SetMarkerSize(1);
  f1->SetLineWidth(4);
  //f1->Draw("same");

  printf("Detectors 0, 6, 12, 18 \n");

  for(int i=0;i<N;i++){
    indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
    indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
    if(i==0)
      indpeak[i]->SetParameter(1,f1->GetParameter(4));
    else
      indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
    indpeak[i]->SetParameter(2,f1->GetParameter(0));
    indpeak[i]->SetLineStyle(7);
    indpeak[i]->Draw("same");
    if (i==0)
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
    else
      printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
  }
  

  //================================= plotting and saving histograms

  if (plothist)  {

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

    TCanvas *cExnogates = new TCanvas("cExnogates", "Ex, no gates", 800, 600);
    Ex_nogates->Draw();
    cExnogates->SaveAs("17O_analysis/Ex_nogates.png");

  }


}
