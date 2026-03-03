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
TString ellLabels[5] = {"#it{l} = 0", "#it{l} = 1", "#it{l} = 2", "#it{l} = 3", "#it{l} = 4"};

void compare_DWBA(){

  TString filename = Form("DWBA_17O_AK.root");
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

  for (int i = 0; i <= 4 && i < graphsDWBA.size(); ++i) {
    graphsDWBA[i]->SetLineColor(colors[i]);
    graphsDWBA[i]->SetLineWidth(2);
  }


  TCanvas* canvas = new TCanvas("cCanvas",
                                "DWBA calculations",
                                1400, 2000);

  canvas->Divide(1, 2, 0.0, 0.0);

  canvas->cd(1);

  TLegend* leg = new TLegend(0.65, 0.55, 0.88, 0.8);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetTextFont(42);

  leg->AddEntry(graphsDWBA[0], ellLabels[0], "l");
  leg->AddEntry(graphsDWBA[1], ellLabels[1], "l");
  leg->AddEntry(graphsDWBA[2], ellLabels[2], "l");
  leg->AddEntry(graphsDWBA[4], ellLabels[4], "l");

  graphsDWBA[0]->Draw("AL");
  graphsDWBA[0]->GetHistogram()->GetXaxis()->SetRangeUser(0, 180);
  graphsDWBA[0]->SetTitle(Form(";;d#sigma/d#Omega (mb/sr)"));

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.0); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.15);

  graphsDWBA[0]->GetXaxis()->SetTitleSize(0.044);
  graphsDWBA[0]->GetYaxis()->SetTitleSize(0.048);

  graphsDWBA[0]->GetXaxis()->SetTitleOffset(1.2);
  graphsDWBA[0]->GetYaxis()->SetTitleOffset(1.4);

  graphsDWBA[0]->GetXaxis()->SetLabelSize(0.0);
  graphsDWBA[0]->GetYaxis()->SetLabelSize(0.045);

  graphsDWBA[0]->GetXaxis()->SetLabelOffset(0.01);
  graphsDWBA[0]->GetYaxis()->SetLabelOffset(0.0098);

  for (int i = 1; i <= 4 && i < graphsDWBA.size(); ++i) {
    if (i != 3)
      graphsDWBA[i]->Draw("L same");
  }

  leg->Draw();

  canvas->cd(2);
  
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);

  TH1F* frame = new TH1F("frame_log", "", 100, 0, 180);
  frame->SetMinimum(1e-3);   // positive minimum for log scale
  frame->SetMaximum(40);    // optional upper limit
  frame->GetXaxis()->SetTitle("#theta_{CM} (deg)");
  frame->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
  frame->Draw();

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.2); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.0);

  frame->GetXaxis()->SetTitleSize(0.048);
  frame->GetYaxis()->SetTitleSize(0.044);

  frame->GetXaxis()->SetTitleOffset(1.28);
  frame->GetYaxis()->SetTitleOffset(1.5);

  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);

  frame->GetXaxis()->SetLabelOffset(0.01);
  frame->GetYaxis()->SetLabelOffset(0.01);

  TLegend* leg2 = new TLegend(0.65, 0.7, 0.88, 0.90);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetTextFont(42);

  leg2->AddEntry(graphsDWBA[0], ellLabels[0], "l");
  leg2->AddEntry(graphsDWBA[1], ellLabels[1], "l");
  leg2->AddEntry(graphsDWBA[2], ellLabels[2], "l");
  leg2->AddEntry(graphsDWBA[4], ellLabels[4], "l");

  leg2->Draw();

  graphsDWBA[0]->Draw("L same");
  graphsDWBA[0]->GetHistogram()->GetXaxis()->SetRangeUser(0, 180);

  for (int i = 1; i < 5 && i < graphsDWBA.size(); ++i) {
    if (i != 3)
      graphsDWBA[i]->Draw("L same");
  }

  gPad->SetLogy();
  gStyle->SetOptStat(0);

  canvas->SaveAs("DWBA_plots2.pdf");
  canvas->SaveAs("DWBA_plots2.png");

  if (graphsDWBA.size() <= 44) {
    std::cerr << "Not enough graphs for c2!" << std::endl;
    return;
  }

  TCanvas* c2 = new TCanvas("c2",
                                "DWBA comp",
                                2200, 1000);
  
  c2->Divide(2, 1);

  c2->cd(1);


  graphsDWBA[26]->SetLineColor(colors[1]);
  graphsDWBA[26]->SetLineWidth(2);
  graphsDWBA[27]->SetLineColor(colors[0]);
  graphsDWBA[27]->SetLineWidth(2);
  graphsDWBA[28]->SetLineColor(colors[2]);
  graphsDWBA[28]->SetLineWidth(2);
  graphsDWBA[44]->SetLineColor(colors[3]);
  graphsDWBA[44]->SetLineWidth(2);

  gPad->SetLogy();
  gStyle->SetOptStat(0);

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.15); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.1);

  TH1F* frame2 = new TH1F("frame2", "", 100, 0, 180);
  frame2->SetMinimum(1e-3);   // positive minimum for log scale
  frame2->SetMaximum(100);
  frame2->GetXaxis()->SetTitle("#theta_{CM} (deg)");
  frame2->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
  frame2->Draw();

  frame2->GetXaxis()->SetTitleSize(0.048);
  frame2->GetYaxis()->SetTitleSize(0.044);
  frame2->GetXaxis()->SetTitleOffset(1.28);
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetLabelOffset(0.01);
  frame2->GetYaxis()->SetLabelOffset(0.01);

  graphsDWBA[27]->Draw("L same");
  graphsDWBA[44]->Draw("L same");

  TLegend* leg3 = new TLegend(0.65, 0.7, 0.88, 0.88);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetTextFont(42);
  leg3->AddEntry(graphsDWBA[27], "#it{l} = 0, 2^{+}", "l");
  leg3->AddEntry(graphsDWBA[44], "#it{l} = 2, 2^{+}", "l");
  leg3->Draw();

  c2->cd(2);
  gPad->SetLogy();
  gStyle->SetOptStat(0);

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.15); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.1);

  TH1F* frame3 = new TH1F("frame3", "", 100, 0, 180);
  frame3->SetMinimum(1e-3);   // positive minimum for log scale
  frame3->SetMaximum(100);
  frame3->GetXaxis()->SetTitle("#theta_{CM} (deg)");
  frame3->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
  frame3->Draw();

  frame3->GetXaxis()->SetTitleSize(0.048);
  frame3->GetYaxis()->SetTitleSize(0.044);
  frame3->GetXaxis()->SetTitleOffset(1.28);
  frame3->GetYaxis()->SetTitleOffset(1.5);
  frame3->GetXaxis()->SetLabelSize(0.04);
  frame3->GetYaxis()->SetLabelSize(0.04);
  frame3->GetXaxis()->SetLabelOffset(0.01);
  frame3->GetYaxis()->SetLabelOffset(0.01);
  
  graphsDWBA[27]->Draw("L same");
  graphsDWBA[28]->Draw("L same");

  TLegend* leg4 = new TLegend(0.65, 0.7, 0.88, 0.88);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);
  leg4->SetTextFont(42);
  leg4->AddEntry(graphsDWBA[27], "#it{l} = 0, 2^{+}", "l");
  leg4->AddEntry(graphsDWBA[28], "#it{l} = 0, 3^{+}", "l");
  leg4->Draw();

  c2->SaveAs("DWBA_comp.pdf");
  c2->SaveAs("DWBA_comp.png");
}
