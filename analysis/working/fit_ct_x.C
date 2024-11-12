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

void saveFitParameters(const TString &fileName, int detIndex, TF1 *fitFunc) {
  std::ofstream outFile;
  outFile.open(fileName, std::ios::app);

  if (outFile.is_open()) {
    // Print detector index and fit parameters (p0, p1, p2, etc.)
    outFile << "Detector " << detIndex << ": ";
    for (int i = 0; i < fitFunc->GetNpar(); i++) {
      outFile << fitFunc->GetParameter(i) << " ";
    }
    outFile << std::endl;
  } else {
    std::cerr << "Failed to open file for saving fit parameters!" << std::endl;
  }

  outFile.close();
}

void fit_ct_x(){

  //========================================= read data files

  TChain * chain = new TChain("tree");
  chain->Add("trace_run055-066.root"); // 17O
  //chain->Add("trace_run022-052.root"); // 17F

  //========================================= detector geometry
  printf("======================= loading parameters files .... \n");
  string detGeoFileName = "detectorGeo.txt";
  printf("loading detector geometery : %s.", detGeoFileName.c_str());

  DetGeo detGeo;

  TMacro * haha = new TMacro();
  if( haha->ReadFile(detGeoFileName.c_str()) > 0 ) {

    detGeo = LoadDetectorGeo(haha);

    PrintDetGeo(detGeo);

    printf("... done.\n");
  }else{
    printf("... fail\n");
    return;
  }

  double length = detGeo.detLength;
  vector<double> pos = detGeo.detPos;

  int colDet = detGeo.nDet;
  int rowDet = detGeo.mDet;

  int nDet = colDet * rowDet;

  delete haha;

  int coinTimeRange[3] = {400, -100, 200}; // bin, min, max
  double xRange[3] = {200, -1.5, 1.5};

  TH2F ** w = new TH2F*[nDet];
  TH2F ** correctedHist = new TH2F*[nDet];

  TFile * cutFile = new TFile("coinTime_x_cuts.root", "recreate");
  TCutG * cut = NULL;
  TObjArray * cutList = new TObjArray();

  //========================================= canvas

  Int_t Div[2] = {colDet,rowDet};  //x,y
  Int_t size[2] = {800, 800}; //x,y
  TCanvas * cCoinTimeX = new TCanvas("cCoinTimeX", "cCoinTimeX", 0, 0, size[0]*Div[0], size[1]*Div[1]);
  cCoinTimeX->Divide(Div[0],Div[1]);

  for( int i = 1; i <= Div[0]*Div[1] ; i++){
    cCoinTimeX->cd(i)->SetGrid();
  }

  gStyle->SetOptStat(1111);
  gStyle->SetStatY(1.0);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);

  if(cCoinTimeX->GetShowEditor()  )cCoinTimeX->ToggleEditor();
  if(cCoinTimeX->GetShowToolBar() )cCoinTimeX->ToggleToolBar();

  std::ofstream fitParamFile("coinTime_x_fit_parameters.txt");

  if (!fitParamFile.is_open()) {
    std::cerr << "Error: Could not open coinTime_x_fit_parameters.txt for writing!" << std::endl;
    return;
  }

  TH2F ** q = new TH2F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("ctx%d", i);
    q[i] = new TH2F(name, name, xRange[0], xRange[1], xRange[2], coinTimeRange[0], coinTimeRange[1], coinTimeRange[2]);
    q[i]->SetXTitle(name);

    TString expression;
    expression.Form("coinTime:x[%d] >> ctx%d" , i, i);

    cCoinTimeX->cd(i+1);
    chain->Draw(expression, "" , "scat");
    cCoinTimeX->Update();
    gSystem->ProcessEvents();

  }

  cCoinTimeX->SaveAs(Form("coinTime_x_17O.pdf"));
  cCoinTimeX->SaveAs(Form("coinTime_x_17O.png"));

  //cCoinTimeX->SaveAs(Form("coinTime_x_17F.pdf"));
  //cCoinTimeX->SaveAs(Form("coinTime_x_17F.png"));

  for (int i = 0; i < nDet; i++) {
    TString canvasName;
    canvasName.Form("cSingleDet%d", i);
    TCanvas *cSingleDet = new TCanvas(canvasName, canvasName, 800, 800);
    if (!cSingleDet->GetShowToolBar()) cSingleDet->ToggleToolBar();

    q[i]->Draw("scat");

    cSingleDet->Modified();
    cSingleDet->Update();

    gPad->WaitPrimitive();

    cut = (TCutG*) gROOT->FindObject("CUTG");

    if (cut) {
      TString name;
      name.Form("cut%d", i);
      cut->SetName(name);
      cut->SetVarX(Form("x[%d]", i));
      cut->SetVarY("coinTime");
      cut->SetLineColor(i + 1);
      cutList->Add(cut);

      printf(" cut-%d \n", i);

      TString fit_type;
      std::cout << "Enter the fitting function type (e.g., \"pol2\", \"pol4\"): ";
      std::cin >> fit_type;

      TF1 *fitFunc = new TF1("fitFunc", fit_type, -0.95, 0.95);
      fitFunc->SetLineColor(kRed);

      TString histName;
      histName.Form("w%d", i);
      w[i] = new TH2F(histName, histName, xRange[0], xRange[1], xRange[2], coinTimeRange[0], coinTimeRange[1], coinTimeRange[2]);
      w[i]->SetXTitle(Form("x[%d]", i));
      w[i]->SetYTitle("coinTime");

      TString expression;
      expression.Form("coinTime:x[%d] >> %s", i, histName.Data());
      chain->Draw(expression, name, "scat");
      cSingleDet->cd();

      w[i]->Fit(fitFunc, "R");
      saveFitParameters("coinTime_x_fit_parameters.txt", i, fitFunc);

      w[i]->Draw("scat");
      fitFunc->Draw("same");

      cSingleDet->Modified();
      cSingleDet->Update();

      cSingleDet->SaveAs(Form("coinTime_x_17O_Det%d_withFit.png", i));

      canvasName.Form("cSingleDetCorr%d", i);
      TCanvas *cSingleDetCorr = new TCanvas(canvasName, canvasName, 800, 800);

      TString correctedHistName;
      correctedHistName.Form("corrected_q%d", i);
      correctedHist[i] = new TH2F(correctedHistName, correctedHistName, xRange[0], xRange[1], xRange[2], coinTimeRange[0], coinTimeRange[1], coinTimeRange[2]);
      correctedHist[i]->SetXTitle(Form("x[%d]", i));
      correctedHist[i]->SetYTitle("coinTime");

      TString correctionExpression = Form("coinTime - ( %f + %f * x[%d] + %f * x[%d]**2 )",
                                          fitFunc->GetParameter(0),
                                          fitFunc->GetParameter(1),
					  i,
                                          fitFunc->GetParameter(2),
					  i);

      if (fitFunc->GetNpar() == 4) {
        correctionExpression = Form("coinTime - ( %f + %f * x[%d] + %f * x[%d]**2 + %f * x[%d]**3 )",
				    fitFunc->GetParameter(0),
				    fitFunc->GetParameter(1),
				    i,
				    fitFunc->GetParameter(2),
				    i,
				    fitFunc->GetParameter(3),
				    i);
      }
      if (fitFunc->GetNpar() == 5) {
        correctionExpression = Form("coinTime - ( %f + %f * x[%d] + %f * x[%d]**2 + %f * x[%d]**3 + %f * x[%d]**4 )",
				    fitFunc->GetParameter(0),
				    fitFunc->GetParameter(1),
				    i,
				    fitFunc->GetParameter(2),
				    i,
				    fitFunc->GetParameter(3),
				    i,
				    fitFunc->GetParameter(4),
				    i);
      }


      expression = Form("(%s) : x[%d] >> %s", correctionExpression.Data(), i, correctedHistName.Data());
      chain->Draw(expression, "", "scat");

      correctedHist[i]->Draw("scat");

      TString canvasName;
      canvasName.Form("rescaled_spectrum_coinTime_x_Det%d.png", i);
      cSingleDetCorr->SaveAs(canvasName);


    } else {
      printf("No cut created for %d-th plot. Skipping.\n", i);
    }
  }

  TCanvas *cCombinedProj = new TCanvas("cCombinedProj", "Combined Projection onto Y-axis", 800, 600);
  cCombinedProj->SetGrid();

  TH1F *combinedProj = nullptr;

  for (int i = 0; i < nDet; i++) {
    if (correctedHist[i]) {  // Make sure the histogram exists
      if (!combinedProj) {
	combinedProj = (TH1F*) correctedHist[i]->ProjectionY("combinedProj", 1, -1);
      } else {
	combinedProj->Add(correctedHist[i]->ProjectionY("", 1, -1));
      }
    }
  }

  if (combinedProj) {
    combinedProj->SetLineColor(kBlack);
    combinedProj->SetLineWidth(2);
    combinedProj->Draw();
  }

  cCombinedProj->SaveAs("combined_projection_Y.png");

  // cSingleDet->SaveAs(Form("coinTime_x_17O_Det%d.png", i));


  delete[] q;
}
