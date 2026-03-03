#include <TH2F.h>
#include <TFile.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TString.h>
#include <TObjArray.h>
#include <TSystem.h>

// make graphic cuts for rdt for 17F or 17O

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
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.15);

  // Marker & line styles
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.25);

  gStyle->SetImageScaling(3.0);
}

void CutCreator(){

  SetThesisStyle();
	
  printf("================ Graphic Cut Creator for RDT ============== \n");

  bool separate = true;
   
  TChain * chain = new TChain("tree");
  //chain->Add("trace_run022-052_corrected.root");   // 17F
  chain->Add("trace_run055-066_corrected.root"); // 17O
  //chain->Add("data/gen_run49.root");
  //chain->Add("data/gen_run50.root");
   
  chain->GetListOfFiles()->Print();
  //gPad->SetEditable(kTRUE);
   
  TString varX, varY, tag;
		
  gStyle->SetOptStat(11111);
	
  TCanvas * cCutCreator = new TCanvas("cCutCreator", "RDT Cut Creator", 100, 100, 800, 800);
  if( !cCutCreator->GetShowToolBar() ) cCutCreator->ToggleToolBar();
	
   
  TFile * cutFile = new TFile("rdtCuts_17O_up_line.root", "recreate");
  cCutCreator->Update();
	
  TCutG * cut = NULL;
  TObjArray * cutList = new TObjArray();
	
   
  TString expression[10];

  int kkk = 0;
  
  if (separate) {
    for (Int_t i = 0; i < 8; i++) {

      if( i % 2 == 0  ){

	printf("======== make a graphic cut on the plot, %d-th cut: ", i );

	varX.Form("rdt[%d]",i);
	varY.Form("rdt[%d]",i+1);

	expression[i].Form("%s:%s>>h(2096, 0, 10000, 2096, 0, 10000)", 
			   varY.Data(),
			   varX.Data());

	chain->Draw(expression[i], "", "colz");

	TH2F *h = (TH2F*)gPad->GetPrimitive("h");
	if (h) {
	  h->GetXaxis()->SetNdivisions(505);   // reduce x labels
	  h->GetXaxis()->SetLabelSize(0.04);  // smaller label size
	  h->GetXaxis()->SetLabelOffset(0.02); // push labels away
	  h->GetYaxis()->SetLabelSize(0.04);
	  h->GetYaxis()->SetLabelOffset(0.02);

	  h->GetXaxis()->SetTitle("E (arbitrary units)");
	  h->GetYaxis()->SetTitle("#DeltaE (arbitrary units)");

	  h->GetXaxis()->SetTitleSize(0.05);
	  h->GetYaxis()->SetTitleSize(0.05);
	  h->GetXaxis()->SetLabelSize(0.04);
	  h->GetYaxis()->SetLabelSize(0.04);

	  h->GetXaxis()->SetTitleOffset(1.3);
	  h->GetYaxis()->SetTitleOffset(1.8);

	  h->SetTitle("");

	  // Draw main title with TLatex
	  TLatex mainTitle;
	  mainTitle.SetTextFont(42);    // Helvetica
	  mainTitle.SetTextSize(0.07);  // Bigger title
	  mainTitle.SetTextAlign(22);   // Centered
	  gPad->SetTopMargin(0.15);     // make space for title

	  TString titleStr;
	  titleStr.Form("#Delta E vs E plot, segment %d", kkk);
	  mainTitle.DrawLatexNDC(0.5, 0.95, titleStr);
	  kkk++;
	}

	gStyle->SetOptStat(000);
     
	cCutCreator->Modified();
	cCutCreator->Update();

	gPad->WaitPrimitive();

	cut = (TCutG*) gROOT->FindObject("CUTG");

	if (cut) {
	  TString name;
	  name.Form("cut%d", i);
	  cut->SetName(name);
	  cut->SetVarX(varX.Data());
	  cut->SetVarY(varY.Data());
	  cut->SetTitle(tag);
	  cut->SetLineColor(kBlack);
	  cutList->Add(cut);

	  printf(" cut-%d \n", i);
	} else {
	  printf("No cut created for %d-th plot. Skipping.\n", i);
	}

	cCutCreator->SaveAs(Form("plots_17O/rdt_%d_%d_17O_up_line.png", i, i+1));
	cCutCreator->SaveAs(Form("plots_17O/rdt_%d_%d_17O_up_line.pdf", i, i+1));
      }

    }
  }

  /*
  TCanvas *cCombined = new TCanvas("cCombined", "Combined (dE vs E) Plot", 100, 100, 800, 800);
  if (!cCombined->GetShowToolBar()) cCombined->ToggleToolBar();

  // Draw all (dE, E) pairs on the same canvas using "colz" for a density plot
  for (int i = 0; i < 4; i++) {
    TString varX, varY;
    varX.Form("rdt[%d]", 2 * i);  // E variable (E1, E2, ...)
    varY.Form("rdt[%d]", 2 * i+1);      // dE variable (dE1, dE2, ...)

    // Draw with "colz" for a density plot
    TString expression;
    expression.Form("%s:%s", varY.Data(), varX.Data());
    chain->Draw(expression, "", i == 0 ? "colz" : "colz same");
  }

  // Create one cut on the combined plot
  cCombined->Update();
  printf("Draw a single cut on the combined plot and double-click to close it.\n");
  gPad->WaitPrimitive();

  TCutG *combinedCut = (TCutG*)gROOT->FindObject("CUTG");

  gStyle->SetOptStat(11111);
  
  if (combinedCut) {
    combinedCut->SetName("combinedCut");
    cutList->Add(combinedCut);
    cutFile->cd();
    combinedCut->Write();
  } else {
    printf("No combined cut created.\n");
  }

  cCutCreator->SaveAs(Form("plots/rdt_combined_tight.png"));
  */
	
  cutList->Write("cutList", TObject::kSingleKey);
	
  printf("====> saved cuts into rdtCuts.root\n");
	
}
