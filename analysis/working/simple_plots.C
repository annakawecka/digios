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

void simple_plots(){

  //========================================= read data files
  
  TChain * chain = new TChain("tree");
  chain->Add("gen_run022-052.root");

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

  int energyRange[3] = {200, -2, 40}; // bin, min, max
  double threshold = 0.2;

  //========================================= canvas

  Int_t Div[2] = {colDet,rowDet};  //x,y
  Int_t size[2] = {800, 800}; //x,y
  TCanvas * cAlpha = new TCanvas("cAlpha", "cAlpha", 0, 0, size[0]*Div[0], size[1]*Div[1]);
  cAlpha->Divide(Div[0],Div[1]);
   
  for( int i = 1; i <= Div[0]*Div[1] ; i++){
    cAlpha->cd(i)->SetGrid();
  }
   
  gStyle->SetOptStat(1111);
  gStyle->SetStatY(1.0);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.1);
   
  if(cAlpha->GetShowEditor()  )cAlpha->ToggleEditor();
  if(cAlpha->GetShowToolBar() )cAlpha->ToggleToolBar();

  TString * gate = new TString[nDet];
  
  //========================================= energy
  TH1F ** q = new TH1F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("q%d", i);
    q[i] = new TH1F(name, name, energyRange[0], energyRange[1], energyRange[2]);
    q[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("e[%d] >> q%d" ,i, i);
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("e[%d] > 0", i);
      
    cAlpha->cd(i+1);
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
  }
   
  //----------- 1, pause for save Canvas
  int dummy = 0;
  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("e.pdf"));

  //========================================= xn
  energyRange[0] = 600;
  energyRange[1] = -9000;
  energyRange[2] = 13000;
  
  TH1F ** r = new TH1F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("r%d", i);
    r[i] = new TH1F(name, name, energyRange[0], energyRange[1], energyRange[2]);
    r[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("xn[%d] >> r%d" ,i, i);
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("e[%d] > 0", i);
      
    cAlpha->cd(i+1);
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
  }

  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("xn.pdf"));

  //========================================= xf
  
  TH1F ** s = new TH1F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("s%d", i);
    s[i] = new TH1F(name, name, energyRange[0], energyRange[1], energyRange[2]);
    s[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("xf[%d] >> s%d" ,i, i);
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("e[%d] > 0", i);
      
    cAlpha->cd(i+1);
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
  }

  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("xf.pdf"));

  //========================================= xn vs xf
  
  TH2F ** t = new TH2F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("t%d", i);
    t[i] = new TH2F(name, name, energyRange[0], energyRange[1], energyRange[2], energyRange[0], energyRange[1], energyRange[2]);
    t[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("xn[%d]:xf[%d] >> t%d", i, i, i);
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("e[%d] > 0", i);
      
    cAlpha->cd(i+1);
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
  }

  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("xn_vs_xf.pdf"));

  //========================================= Ex
  
  TH1F ** u = new TH1F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("u%d", i);
    u[i] = new TH1F(name, name, 200, 0., 10.);
    u[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("Ex >> u%d" , i);
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    gate[i].Form("detID == %d", i);
      
    cAlpha->cd(i + 1)->Clear();
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
  }

  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("Ex.pdf"));
  
  //========================================= e vs z

  TH2F ** v = new TH2F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("v%d", i);
    v[i] = new TH2F(name, name, 400, -600, -100, 100, 0., 35.);
    v[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("e[%d]:z[%d] >> v%d", i, i, i);
    gate[i].Form("");
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("z[%d] > -600 && z[%d] < -100", i, i);
      
    cAlpha->cd(i + 1)->Clear();
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
  }

  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("e_vs_z.pdf"));

  //========================================= e vs z close-up
  int jj = 0;
  TH2F ** w = new TH2F*[nDet];
  for( int i = 0; i < nDet; i ++){
    if (jj == 6)
      jj = 0;
    TString name;
    name.Form("w%d", i);
    w[i] = new TH2F(name, name, 400, pos[jj]-100., pos[jj]+50., 100, 0., 35.);
    w[i]->SetXTitle(name);
      
    TString expression;
    expression.Form("e[%d]:z[%d] >> w%d", i, i, i);
    gate[i].Form("");
    //gate[i].Form("ring[%d]==0 && !TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i, i);
    //gate[i].Form("!TMath::IsNaN(xf[%d]) && !TMath::IsNaN(xn[%d])", i, i);
    //gate[i].Form("z[%d] > -600 && z[%d] < -100", i, i);
      
    cAlpha->cd(i + 1)->Clear();
    chain->Draw(expression, gate[i] , "");
    cAlpha->Update();
    gSystem->ProcessEvents();
    jj++;
  }

  cAlpha->Update();
  gSystem->ProcessEvents();

  cAlpha->SaveAs(Form("e_vs_z_closeup.pdf"));

  for (int i = 0; i < nDet; i++) {
    delete q[i];
    delete r[i];
    delete s[i];
    delete t[i];
    delete u[i];
    delete v[i];
    delete w[i];
  }
  
  delete[] q;
  delete[] r;
  delete[] s;
  delete[] t;
  delete[] u;
  delete[] v;
  delete[] w;

}
