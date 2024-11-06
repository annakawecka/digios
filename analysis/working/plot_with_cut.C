#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCutG.h>

TString rdtCutFile1 = "rdtCuts.root";

TObjArray * cutList1;
Bool_t isCutFileOpen1;
int numCut1;

TCutG* cutG1;

TString rdtCutFile2 = "rdtCuts_tight.root";

TObjArray * cutList2;
Bool_t isCutFileOpen2;
int numCut2;

TCutG* cutG2;

void plot_with_cut() {

  //================================= Cut 1
  
  TFile * fCut1 = new TFile(rdtCutFile1);
  isCutFileOpen1 = fCut1->IsOpen();
  if(!isCutFileOpen1) {
    printf( "Failed to open rdt-cutfile 1 : %s\n" , rdtCutFile1.Data());
    rdtCutFile1 = "";
  }
  numCut1 = 0 ;

  if( isCutFileOpen1 ){
    cutList1 = (TObjArray *) fCut1->FindObjectAny("cutList");
    numCut1 = cutList1->GetEntries();
    printf("=========== found %d cutG in %s \n", numCut1, fCut1->GetName());

    for(int i = 0; i < numCut1 ; i++){
      printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
	     cutList1->At(i)->GetName(),
	     ((TCutG*)cutList1->At(i))->GetVarX(),
	     ((TCutG*)cutList1->At(i))->GetVarY(),
	     ((TCutG*)cutList1->At(i))->GetN());
    }
  }

  cutG1 = (TCutG *)cutList1->At(4);

  //================================= Cut 2 (tight)

  TFile * fCut2 = new TFile(rdtCutFile2);
  isCutFileOpen2 = fCut2->IsOpen();
  if(!isCutFileOpen2) {
    printf( "Failed to open rdt-cutfile 2 : %s\n" , rdtCutFile2.Data());
    rdtCutFile2 = "";
  }
  numCut2 = 0 ;

  if( isCutFileOpen2 ){
    cutList2 = (TObjArray *) fCut2->FindObjectAny("cutList");
    numCut2 = cutList2->GetEntries();
    printf("=========== found %d cutG in %s \n", numCut2, fCut2->GetName());

    for(int i = 0; i < numCut2 ; i++){
      printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
	     cutList2->At(i)->GetName(),
	     ((TCutG*)cutList2->At(i))->GetVarX(),
	     ((TCutG*)cutList2->At(i))->GetVarY(),
	     ((TCutG*)cutList2->At(i))->GetN());
    }
  }

  cutG2 = (TCutG *)cutList2->At(4) ;

  //================================= reading the tree

  TChain *chain = new TChain("tree");
  chain->Add("gen_run022-052.root");

  TCanvas *c = new TCanvas("c", "Plot Ex with Cut", 800, 800);

  TH1F *histEx1 = new TH1F("histEx1", "Ex with Cuts1", 200, -2, 12);
  TH1F *histEx2 = new TH1F("histEx2", "Ex with Cuts2", 200, -2, 12);

  // Set the branch addresses
  Float_t rdt[8], Ex;
  chain->SetBranchAddress("rdt", rdt);
  chain->SetBranchAddress("Ex", &Ex);

  // Loop over the events and apply the cut
  Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    chain->GetEntry(i);

    if (cutG1->IsInside(rdt[0], rdt[1]) ||
	cutG1->IsInside(rdt[2], rdt[3]) ||
	cutG1->IsInside(rdt[4], rdt[5]) ||
	cutG1->IsInside(rdt[6], rdt[7])) {

      histEx1->Fill(Ex);
    }

    if (cutG2->IsInside(rdt[0], rdt[1]) ||
	cutG2->IsInside(rdt[2], rdt[3]) ||
	cutG2->IsInside(rdt[4], rdt[5]) ||
	cutG2->IsInside(rdt[6], rdt[7])) {

      histEx2->Fill(Ex);
    }
  }

  histEx1->Draw();

  // Optionally, save the plot
  c->SaveAs("plots/Ex_with_cut1.png");

  histEx2->Draw();

  c->SaveAs("plots/Ex_with_cut2.png");

}
