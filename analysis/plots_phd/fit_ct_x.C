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

// code for correcting coincidence time vs x
// fits coinTime vs x and plots the corrected spectrum

void SetThesisStyle() {
  gStyle->SetOptStat(0);                 // Hide stat box
  gStyle->SetOptFit(0);                  // Hide fit box unless wanted

  gStyle->SetPalette(kRainBow);          // Nicer color palette (choose one)

  // Canvas
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);

  // Pads
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Frame
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameBorderMode(0);

  // Fonts
  gStyle->SetTextFont(42);               // Helvetica (thesis-friendly)
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");

  // Margins
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.05);

  // Marker & line styles
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.25);

  gStyle->SetImageScaling(3.0);
}

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

  int colors[] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+2,
                kOrange+7, kAzure+4, kTeal+4, kPink+7};

  bool oxygen = true;
  
  std::string isotope = oxygen ? "17O" : "17F";
  std::string folderName = "plots_" + isotope + "/";

  //========================================= read data files

  TChain * chain = new TChain("tree");
  if (oxygen)
    chain->Add("trace_run055-066_corrected.root"); // 17O
  else
    chain->Add("trace_run022-052_corrected.root"); // 17F

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

  std::string filePath = folderName + "coinTime_x_cuts.root";

  TFile * cutFile = new TFile(filePath.c_str(), "recreate");
  TCutG * cut = NULL;
  TObjArray * cutList = new TObjArray();

  //========================================= canvas

  SetThesisStyle();

  Int_t Div[2] = {colDet,rowDet};  //x,y
  Int_t size[2] = {1200, 1200}; //x,y
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

  filePath = folderName + "coinTime_x_fit_parameters_test.txt";

  std::ofstream fitParamFile(filePath.c_str());

  if (!fitParamFile.is_open()) {
    std::cerr << "Error: Could not open coinTime_x_fit_parameters.txt for writing!" << std::endl;
    return;
  }

  TH2F ** q = new TH2F*[nDet];
  for( int i = 0; i < nDet; i ++){
    TString name;
    name.Form("ctx%d", i);
    q[i] = new TH2F(name, name, xRange[0], xRange[1], xRange[2], coinTimeRange[0], coinTimeRange[1], coinTimeRange[2]);

    //q[i]->SetTitle(Form("Detector %d; Position x_{%d} (relative units); Coincidence Time (ns)", i, i));

    q[i]->SetXTitle(Form("Position x_{%d} (relative units)", i));
    q[i]->SetYTitle("Coincidence Time (ns)");
    q[i]->SetTitle("");
    
    //q[i]->GetXaxis()->SetLineWidth(1); // x-axis line
    //q[i]->GetYaxis()->SetLineWidth(1); // y-axis line
    q[i]->GetXaxis()->SetTickLength(0.02);
    q[i]->GetYaxis()->SetTickLength(0.02);
    q[i]->GetXaxis()->SetLabelSize(0.045);
    q[i]->GetYaxis()->SetLabelSize(0.045);
    q[i]->GetXaxis()->SetLabelOffset(0.02); // default ~0.005–0.01
    q[i]->GetYaxis()->SetLabelOffset(0.02);
    q[i]->GetXaxis()->SetTitleOffset(1.3);
    q[i]->SetTitleSize(0.1, "");

    TString expression;
    expression.Form("coinTime:x[%d] >> ctx%d" , i, i);

    cCoinTimeX->cd(i+1);
    gStyle->SetOptStat(000);
    chain->Draw(expression, "" , "scat");
    cCoinTimeX->Update();
    gSystem->ProcessEvents();

    TLatex mainTitle;
    mainTitle.SetTextFont(42);   // Helvetica
    mainTitle.SetTextSize(0.08); // Bigger main title
    mainTitle.SetTextAlign(22);  // Centered
    gPad->SetTopMargin(0.1);    // Increase top margin to make space for title
    mainTitle.DrawLatexNDC(0.5, 0.95, Form("Detector %d", i));

  }
  

  cCoinTimeX->SaveAs((folderName + "coinTime_x_" + isotope + ".pdf").c_str());
  cCoinTimeX->SaveAs((folderName + "coinTime_x_" + isotope + ".png").c_str());

  for (int i = 0; i < nDet; i++) {
    TString canvasName;
    canvasName.Form("cSingleDet%d", i);
    TCanvas *cSingleDet = new TCanvas(canvasName, canvasName, 800, 600);
    if (!cSingleDet->GetShowToolBar()) cSingleDet->ToggleToolBar();

    q[i]->Draw("scat");

    

    TLatex mainTitle;
    mainTitle.SetTextFont(42);   // Helvetica
    mainTitle.SetTextSize(0.08); // Bigger main title
    mainTitle.SetTextAlign(22);  // Centered
    gPad->SetTopMargin(0.15);    // Increase top margin to make space for title
    mainTitle.DrawLatexNDC(0.5, 0.95, Form("Detector %d", i));

    cSingleDet->Modified();
    cSingleDet->Update();

    //cSingleDet->SaveAs(Form((folderName + "coinTime_x_" + isotope + "_Det%d_noFit.png").c_str(), i), "png300");
    cSingleDet->SaveAs(Form("%scoinTime_x_%s_Det%d_noFit.png", folderName.c_str(), isotope.c_str(), i));
    cSingleDet->SaveAs(Form("%scoinTime_x_%s_Det%d_noFit.pdf", folderName.c_str(), isotope.c_str(), i));

    //gPad->WaitPrimitive();
    cSingleDet->WaitPrimitive("CUTG");
    

    cut = (TCutG*) gROOT->FindObject("CUTG");

    if (cut) {
      TString name;
      name.Form("cut%d", i);
      cut->SetName(name);
      cut->SetVarX(Form("x[%d]", i));
      cut->SetVarY("coinTime");
      cut->SetLineColor(colors[i % 9]);
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
      //w[i]->SetTitle(Form("Detector %d; Position x_{%d} (relative units); Coincidence Time (ns)", i, i));

      w[i]->SetXTitle(Form("Position x_{%d} (relative units)", i));
      w[i]->SetYTitle("Coincidence Time (ns)");
      //w[i]->GetXaxis()->SetLineWidth(1); // x-axis line
      //w[i]->GetYaxis()->SetLineWidth(1); // y-axis line
      w[i]->GetXaxis()->SetTickLength(0.02);
      w[i]->GetYaxis()->SetTickLength(0.02);
      w[i]->GetXaxis()->SetLabelSize(0.045);
      w[i]->GetYaxis()->SetLabelSize(0.045);
      w[i]->GetXaxis()->SetLabelOffset(0.02); // default ~0.005–0.01
      w[i]->GetYaxis()->SetLabelOffset(0.02);
      w[i]->GetXaxis()->SetTitleOffset(1.3);
      w[i]->SetTitle("");

      TString expression;
      expression.Form("coinTime:x[%d] >> %s", i, histName.Data());
      chain->Draw(expression, name, "scat");
      cSingleDet->cd();

      w[i]->Fit(fitFunc, "R");
      saveFitParameters((folderName + "coinTime_x_fit_parameters.txt").c_str(), i, fitFunc);

      w[i]->Draw("scat");
      fitFunc->Draw("same");

      TLatex mainTitle;
      mainTitle.SetTextFont(42);   // Helvetica
      mainTitle.SetTextSize(0.08); // Bigger main title
      mainTitle.SetTextAlign(22);  // Centered
      gPad->SetTopMargin(0.15);    // Increase top margin to make space for title
      mainTitle.DrawLatexNDC(0.5, 0.95, Form("Detector %d", i));

      cSingleDet->Modified();
      cSingleDet->Update();

      //cSingleDet->SaveAs(Form("coinTime_x_17O_Det%d_withFit.png", i));
      //cSingleDet->SaveAs(Form((folderName + "coinTime_x_" + isotope + "_Det%d_withFit.png").c_str(), i), "png300");
      cSingleDet->SaveAs(Form("%scoinTime_x_%s_Det%d_withFit.png", folderName.c_str(), isotope.c_str(), i));
      cSingleDet->SaveAs(Form("%scoinTime_x_%s_Det%d_withFit.pdf", folderName.c_str(), isotope.c_str(), i));

      canvasName.Form("cSingleDetCorr%d", i);
      TCanvas *cSingleDetCorr = new TCanvas(canvasName, canvasName, 800, 600);

      TString correctedHistName;
      correctedHistName.Form("corrected_q%d", i);
      correctedHist[i] = new TH2F(correctedHistName, correctedHistName, xRange[0], xRange[1], xRange[2], coinTimeRange[0], coinTimeRange[1], coinTimeRange[2]);
      correctedHist[i]->SetXTitle(Form("x[%d]", i));
      correctedHist[i]->SetYTitle("coinTime");
      //correctedHist[i]->SetTitle(Form("Corrected Coincidence Time - Detector %d; x_%d; Corrected CoinTime (ns)", i, i));

      correctedHist[i]->SetXTitle(Form("Position x_{%d} (relative units)", i));
      correctedHist[i]->SetYTitle("Corrected CoinTime (ns)");

      
      
      //correctedHist[i]->GetXaxis()->SetLineWidth(1); // x-axis line
      //correctedHist[i]->GetYaxis()->SetLineWidth(1); // y-axis line
      correctedHist[i]->GetXaxis()->SetTickLength(0.02);
      correctedHist[i]->GetYaxis()->SetTickLength(0.02);
      correctedHist[i]->GetXaxis()->SetLabelSize(0.045);
      correctedHist[i]->GetYaxis()->SetLabelSize(0.045);
      correctedHist[i]->GetXaxis()->SetLabelOffset(0.02); // default ~0.005–0.01
      correctedHist[i]->GetYaxis()->SetLabelOffset(0.02);
      correctedHist[i]->GetXaxis()->SetTitleOffset(1.3);
      correctedHist[i]->SetTitle("");


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

      mainTitle.SetTextFont(42);   // Helvetica
      mainTitle.SetTextSize(0.08); // Bigger main title
      mainTitle.SetTextAlign(22);  // Centered
      gPad->SetTopMargin(0.15);    // Increase top margin to make space for title
      mainTitle.DrawLatexNDC(0.5, 0.95, Form("Corrected Coincidence Time - Detector %d", i));

      

      TString canvasName;
      canvasName.Form((folderName + "rescaled_spectrum_coinTime_x_Det%d.png").c_str(), i);
      cSingleDetCorr->SaveAs(canvasName);
      canvasName.Form((folderName + "rescaled_spectrum_coinTime_x_Det%d.pdf").c_str(), i);
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
    //combinedProj->SetTitle("Combined Corrected Coincidence Time;Time (ns);Counts");
    combinedProj->SetXTitle(Form("Time (ns)"));
    combinedProj->SetYTitle("Counts");
    combinedProj->SetTitle("");
    
    combinedProj->SetLineWidth(3);
    combinedProj->SetFillColorAlpha(kBlue, 0.2);
    combinedProj->GetXaxis()->SetTitleOffset(1.3);
    combinedProj->Draw();

    TLatex mainTitle;
    mainTitle.SetTextFont(42);   // Helvetica
    mainTitle.SetTextSize(0.08); // Bigger main title
    mainTitle.SetTextAlign(22);  // Centered
    gPad->SetTopMargin(0.15);    // Increase top margin to make space for title
    mainTitle.DrawLatexNDC(0.5, 0.95, Form("Combined Corrected Coincidence Time"));
    
  }

  cCombinedProj->SaveAs((folderName + "combined_projection_Y.png").c_str());
  cCombinedProj->SaveAs((folderName + "combined_projection_Y.pdf").c_str());

  // cSingleDet->SaveAs(Form("coinTime_x_17O_Det%d.png", i));


  delete[] q;
}
