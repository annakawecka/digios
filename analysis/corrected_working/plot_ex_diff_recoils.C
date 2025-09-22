#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>

void plot_ex_diff_recoils() {
  TFile* file1 = TFile::Open("plots_17F/Ex_x_recoil_coinTime_gated.root");
  TFile* file2 = TFile::Open("plots_17O/Ex_x_recoil_coinTime_gated.root");
  TFile* file3 = TFile::Open("plots_17F/Ex_x_diffrecoil_coinTime_gated.root");

  if (!file1 || !file2) {
    std::cerr << "Error opening files" << std::endl;
    return;
  }

  TH1* hist1_orig = dynamic_cast<TH1*>(file1->Get("x_rdt_coinTime_gatedEx"));
  TH1* hist2_orig = dynamic_cast<TH1*>(file2->Get("x_rdt_coinTime_gatedEx"));
  TH1* hist3_orig = dynamic_cast<TH1*>(file3->Get("x_diffrdt_coinTime_gatedEx"));

  if (!hist1_orig || !hist2_orig || !hist3_orig) {
    std::cerr << "Error retrieving histograms" << std::endl;
    return;
  }

  TH1* hist1 = (TH1*)hist1_orig->Clone("hist1");
  hist1->SetDirectory(0);
  TH1* hist2 = (TH1*)hist2_orig->Clone("hist2");
  hist2->SetDirectory(0);
  TH1* hist3 = (TH1*)hist3_orig->Clone("hist2");
  hist3->SetDirectory(0);


  Int_t kMyBlue  = kViolet + 8; //TColor::GetColor(31, 119, 180);  // #1f77b4
  Int_t kMyRed   = TColor::GetColor(255, 127, 14);  // #ff7f0e
  Int_t kMyGreen = TColor::GetColor(44, 160, 44);   // #2ca02c

  Int_t kBrightBlue  = TColor::GetColor(0, 102, 204);   // mocny niebieski
  Int_t kBrightRed   = TColor::GetColor(255, 51, 0);    // intensywny czerwony
  Int_t kBrightGreen = TColor::GetColor(0, 153, 51);

  Int_t kVividBlue = TColor::GetColor(0, 128, 255);
  Int_t kDarkBlue = TColor::GetColor(0, 102, 204);

  file1->Close();
  file2->Close();
  file3->Close();

  TH1* hist2_shifted = (TH1*)hist2->Clone("hist2_shifted");
  hist2_shifted->Reset();

  double shift = 1.04155; // shifted to match the 0+ T=1 state 
  for (int i = 1; i <= hist2->GetNbinsX(); ++i) {
    double x = hist2->GetBinCenter(i);
    double y = hist2->GetBinContent(i);
    double ex = hist2->GetBinError(i);
    int new_bin = hist2_shifted->FindBin(x + shift);

    if (new_bin >= 1 && new_bin <= hist2_shifted->GetNbinsX()) {
      hist2_shifted->SetBinContent(new_bin, y);
      hist2_shifted->SetBinError(new_bin, ex);
    }
  }

  TLegend *leg = new TLegend(0.55, 0.7, 0.92, 0.88);
  leg->AddEntry(hist1, "17F(d,p) with 18F recoils", "l");
  leg->AddEntry(hist2_shifted, "17O(d,p) with 18O recoils", "l");
  leg->AddEntry(hist3, "17F(d,p) with 17O recoils", "l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.038);

  TCanvas* canvas = new TCanvas("canvas", "Histograms", 2400, 1600);
  canvas->SetCanvasSize(2400, 1600);
  canvas->cd();
  
  hist1->SetLineColor(kDarkBlue); 		 // 17F
  hist2_shifted->SetLineColor(kBrightGreen); // 17O
  hist2->SetLineColor(kBlue);
  hist3->SetLineColor(kBrightRed); // 17F with 17O recoils

  hist1->SetLineStyle(1);//(9); 		 // 17F
  hist2_shifted->SetLineStyle(1);//(7); // 17O
  hist3->SetLineStyle(1); // 17F with 17O recoils

  hist1->SetLineWidth(4); 		 // 17F
  hist2_shifted->SetLineWidth(4); // 17O
  hist3->SetLineWidth(4); // 17F with 17O recoils

  std::cout << hist1->Integral() << std::endl << hist2->Integral()  << std::endl << hist2_shifted->Integral() << std::endl;

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.15); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.1);

  hist2_shifted->GetXaxis()->SetTitleSize(0.044);
  hist2_shifted->GetYaxis()->SetTitleSize(0.044);

  hist2_shifted->GetXaxis()->SetTitleOffset(1.2);
  hist2_shifted->GetYaxis()->SetTitleOffset(1.5);

  hist2_shifted->GetXaxis()->SetLabelSize(0.04);
  hist2_shifted->GetYaxis()->SetLabelSize(0.04);

  hist2_shifted->GetXaxis()->SetLabelOffset(0.01);
  hist2_shifted->GetYaxis()->SetLabelOffset(0.01);
    
  hist2_shifted->Draw("HIST");
  //hist2->Draw("HIST");
  hist1->Draw("HIST SAME");
  hist3->Draw("HIST SAME");

  hist2_shifted->SetTitle("Excitation energy spectrum;E_{x} (MeV);Counts / 70 keV");

  gStyle->SetOptStat(000);

  std::cout << "number of bins: " << hist1->GetNbinsX() << " bins/range = " << 14000. / hist1->GetNbinsX() <<  std::endl;

  double x_line = 5.6071; // 5.6 MeV
  double y_min = gPad->GetUymin();
  double y_max = gPad->GetUymax();

  // Tworzenie linii
  TLine *vline = new TLine(x_line, y_min, x_line, y_max);
  vline->SetLineColor(kBlack);
  vline->SetLineStyle(2); // kreskowana
  vline->SetLineWidth(4);
  vline->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextFont(42);          // 42 = Helvetica, 62 = bold Helvetica
  latex.SetTextColor(kRed);
  latex.SetTextColorAlpha(kGray+2, 0.4);
  latex.DrawLatex(0.7, 0.4, "PRELIMINARY");

  leg->AddEntry(vline, "S_{p} of 18F = 5.61 MeV", "l"); // "l" = linia

  leg->Draw();

  canvas->SaveAs("ex_18O_18F_17O_4.png");
}
