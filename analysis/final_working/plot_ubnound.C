#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <TLatex.h>

std::vector<TGraph*> graphsDWBA;
std::vector<TGraph*> graphsDWBA2;
TString ellLabels[6] = {"1.041 MeV", "3.062 MeV", "4.116 MeV", "4.753 MeV", "4.963 MeV", "6.136 MeV"};

void plot_ubnound(){

  TString filename = Form("../corrected_working/DWBA_17F_AK.root");
  TFile *file;
  file = TFile::Open(filename);

  if (!file || file->IsZombie()) {
    std::cerr << "Error opening ROOT file!" << std::endl;
    return;
  }

  TObjArray* objArray = (TObjArray*)file->Get("qList");
  if (!objArray) {
    std::cerr << "Error: Object array not found in file " << std::endl;
    file->Close();
    return;
  }

  for (int i = 0; i < objArray->GetEntries(); ++i) {
    TObject* obj = objArray->At(i);
    if (obj && obj->InheritsFrom(TGraph::Class())) {
      TGraph* graph = (TGraph*)obj;
      graphsDWBA.push_back(graph);

      std::cout << "Found TGraph: " << graph->GetName() << std::endl;
    }
  }

  std::vector<int> colors = {629, 8, 596, 418, 801, 905, 8, 9, 1};

  graphsDWBA[1]->SetLineColor(colors[0]); // 1.041
  graphsDWBA[5]->SetLineColor(colors[1]); // 3.062
  graphsDWBA[20]->SetLineColor(colors[8]); // 4.116
  graphsDWBA[14]->SetLineColor(colors[3]); // 4.753
  graphsDWBA[2]->SetLineColor(colors[4]); // 4.963
  graphsDWBA[1]->SetLineWidth(2);
  graphsDWBA[5]->SetLineWidth(2);
  graphsDWBA[20]->SetLineWidth(2);
  graphsDWBA[14]->SetLineWidth(2);
  graphsDWBA[2]->SetLineWidth(2);

  TString filename2 = Form("../corrected_working/DWBA_17F_AK_he.root");
  TFile *file2;
  file2 = TFile::Open(filename2);

  if (!file2 || file2->IsZombie()) {
    std::cerr << "Error opening ROOT file!" << std::endl;
    return;
  }

  TObjArray* objArray2 = (TObjArray*)file2->Get("qList");
  if (!objArray2) {
    std::cerr << "Error: Object array not found in file " << std::endl;
    file2->Close();
    return;
  }

  for (int i = 0; i < objArray2->GetEntries(); ++i) {
    TObject* obj = objArray2->At(i);
    if (obj && obj->InheritsFrom(TGraph::Class())) {
      TGraph* graph = (TGraph*)obj;
      graphsDWBA2.push_back(graph);

      std::cout << "Found TGraph: " << graph->GetName() << std::endl;
    }
  }

  graphsDWBA2[1]->SetLineColor(colors[7]); // 6.136
  graphsDWBA2[1]->SetLineWidth(2);

  TCanvas* canvas = new TCanvas("cCanvas",
                                "DWBA calculations",
                                1400, 1000);

  TLegend* leg = new TLegend(0.65, 0.2, 0.88, 0.5);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);

  leg->AddEntry(graphsDWBA[1], ellLabels[0], "l");
  leg->AddEntry(graphsDWBA[5], ellLabels[1], "l");
  leg->AddEntry(graphsDWBA[20], ellLabels[2], "l");
  leg->AddEntry(graphsDWBA[14], ellLabels[3], "l");
  leg->AddEntry(graphsDWBA[2], ellLabels[4], "l");
  leg->AddEntry(graphsDWBA2[1], ellLabels[5], "l");

  graphsDWBA[1]->Draw("AL");
  graphsDWBA[1]->GetHistogram()->GetXaxis()->SetRangeUser(0, 40);
  graphsDWBA[1]->SetTitle(Form(";#theta_{CM} (deg);d#sigma/d#Omega (mb/sr)"));

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.15); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.15);

  graphsDWBA[1]->GetXaxis()->SetTitleSize(0.044);
  graphsDWBA[1]->GetYaxis()->SetTitleSize(0.048);

  graphsDWBA[1]->GetXaxis()->SetTitleOffset(1.2);
  graphsDWBA[1]->GetYaxis()->SetTitleOffset(1.4);

  graphsDWBA[1]->GetXaxis()->SetLabelSize(0.045);
  graphsDWBA[1]->GetYaxis()->SetLabelSize(0.045);

  graphsDWBA[1]->GetXaxis()->SetLabelOffset(0.01);
  graphsDWBA[1]->GetYaxis()->SetLabelOffset(0.0098);

  graphsDWBA[5]->Draw("L same");
  graphsDWBA[20]->Draw("L same");
  graphsDWBA[14]->Draw("L same");
  graphsDWBA[2]->Draw("L same");
  graphsDWBA2[1]->Draw("L same");

  leg->Draw();
  gPad->SetLogy();

  canvas->SaveAs("DWBA_unbound.pdf");
  canvas->SaveAs("DWBA_unbound.png");
}
