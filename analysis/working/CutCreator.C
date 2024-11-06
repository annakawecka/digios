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

void CutCreator(){
	
  printf("================ Graphic Cut Creator for RDT ============== \n");

  bool separate = true;
   
  TChain * chain = new TChain("tree");
  chain->Add("gen_run022-052.root");
  //chain->Add("data/gen_run49.root");
  //chain->Add("data/gen_run50.root");
   
  chain->GetListOfFiles()->Print();
  //gPad->SetEditable(kTRUE);
   
  TString varX, varY, tag;
		
  gStyle->SetOptStat(11111);
	
  TCanvas * cCutCreator = new TCanvas("cCutCreator", "RDT Cut Creator", 100, 100, 800, 800);
  if( !cCutCreator->GetShowToolBar() ) cCutCreator->ToggleToolBar();
	
   
  TFile * cutFile = new TFile("rdtCuts_tight.root", "recreate");
  cCutCreator->Update();
	
  TCutG * cut = NULL;
  TObjArray * cutList = new TObjArray();
	
   
  TString expression[10];
  
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
	  cut->SetLineColor(i+1);
	  cutList->Add(cut);

	  printf(" cut-%d \n", i);
	} else {
	  printf("No cut created for %d-th plot. Skipping.\n", i);
	}

	cCutCreator->SaveAs(Form("plots/rdt_%d_%d_tight.png", i, i+1));
      }

    }
  }
  
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
	
  cutList->Write("cutList", TObject::kSingleKey);
	
  printf("====> saved cuts into rdtCuts.root\n");
	
}
