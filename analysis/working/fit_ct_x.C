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

  TH2F ** q = new TH2F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("ctx%d", i);
    q[i] = new TH2F(name, name, xRange[0], xRange[1], xRange[2], coinTimeRange[0], coinTimeRange[1], coinTimeRange[2]);
    q[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("coinTime:x[%d] >> ctx%d" , i, i);
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("e[%d] > 0", i);
      
    cCoinTimeX->cd(i+1);
    chain->Draw(expression, "" , "scat");
    cCoinTimeX->Update();
    gSystem->ProcessEvents();

  }

  cCoinTimeX->SaveAs(Form("coinTime_x_17O.pdf"));
  cCoinTimeX->SaveAs(Form("coinTime_x_17O.png"));
}
