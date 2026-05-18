#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>

void plot_F_O() {
  TFile* file1 = TFile::Open("../corrected_working/plots_17F/Ex_x_recoil_coinTime_gated.root");
  TFile* file2 = TFile::Open("../corrected_working/plots_17O/Ex_x_recoil_coinTime_gated.root");
  TFile* file3 = TFile::Open("../corrected_working/plots_17F/Ex_x_diffrecoil_coinTime_gated.root");
  TFile* file4 = TFile::Open("../corrected_working/plots_17F/Ex_x_Nrecoil_coinTime_gated.root");

  if (!file1 || !file2 || !file3 || !file4) {
    std::cerr << "Error opening files" << std::endl;
    return;
  }

  TH1* hist1_orig = dynamic_cast<TH1*>(file1->Get("x_rdt_coinTime_gatedEx"));
  TH1* hist2_orig = dynamic_cast<TH1*>(file2->Get("x_rdt_coinTime_gatedEx"));
  TH1* hist3_orig = dynamic_cast<TH1*>(file3->Get("x_diffrdt_coinTime_gatedEx"));
  TH1* hist4_orig = dynamic_cast<TH1*>(file4->Get("x_Nrdt_coinTime_gatedEx"));

  if (!hist1_orig || !hist2_orig || !hist3_orig || !hist4_orig) {
    std::cerr << "Error retrieving histograms" << std::endl;
    return;
  }

  TH1* hist1 = (TH1*)hist1_orig->Clone("hist1"); // 17F beam, 18F recoils
  hist1->SetDirectory(0);
  TH1* hist2 = (TH1*)hist2_orig->Clone("hist2"); // 17O beam, 18O recoils
  hist2->SetDirectory(0);
  TH1* hist3 = (TH1*)hist3_orig->Clone("hist2");  // 17F beam, 17O recoils
  hist3->SetDirectory(0);
  TH1* hist4 = (TH1*)hist4_orig->Clone("hist2");  // 17F beam, 17O recoils
  hist4->SetDirectory(0);
  TH1* hist5 = (TH1*)hist1_orig->Clone("hist5"); // 17F beam, 18F recoils
  hist5->SetDirectory(0);


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
  file4->Close();

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

  TLegend *leg = new TLegend(0.7, 0.67, 0.92, 0.898);
  leg->AddEntry(hist1, "^{17}F(d,p)", "l");
  leg->AddEntry(hist2_shifted, "^{17}O(d,p)", "l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.038);

  TCanvas* canvas = new TCanvas("canvas", "Histograms", 2400, 1600);
  canvas->SetCanvasSize(2400, 1200*2);
  canvas->Divide(1, 2, 0, 0);
  canvas->cd(1);

  double norm_scale = 3.423160411;

  hist1->Scale(norm_scale);
  hist3->Scale(norm_scale);
  hist3->Scale(1.42); // efficiecy of detection of recoil after p emission
  //hist2_shifted->Scale(1./norm_scale);
  hist4->Scale(norm_scale);
  hist5->Scale(norm_scale);
  hist4->Scale(1.42);

  hist1->Add(hist3);
  hist1->Add(hist4);
  
  hist1->SetLineColor(kBrightGreen);//(kDarkBlue); 		 // 17F
  hist2_shifted->SetLineColor(kBrightRed);//(kBrightGreen); // 17O
  hist2->SetLineColor(kBlue);
  hist3->SetLineColor(kDarkBlue); // 17F with 17O recoils
  hist4->SetLineColor(kBlack);
  hist4->SetLineColorAlpha(kBlack, 0.2);
  hist5->SetLineColor(kGreen);

  hist1->SetLineStyle(1);//(9); 		 // 17F
  hist2_shifted->SetLineStyle(1);//(7); // 17O
  hist3->SetLineStyle(1); // 17F with 17O recoils
  hist4->SetLineStyle(1);
  hist5->SetLineStyle(1);

  hist2_shifted->GetYaxis()->SetRangeUser(0.01, 1870);

  hist1->SetLineWidth(1); 		 // 17F
  hist2_shifted->SetLineWidth(1); // 17O
  hist3->SetLineWidth(1); // 17F with 17O recoils
  hist4->SetLineWidth(1);

  std::cout << hist1->Integral() << std::endl << hist2->Integral()  << std::endl << hist2_shifted->Integral() << std::endl;

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.0); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.01);  // wąski prawy margines
  gPad->SetTopMargin(0.1);

  hist2_shifted->GetXaxis()->SetTitleSize(0.044);
  hist2_shifted->GetYaxis()->SetTitleSize(0.044);

  hist2_shifted->GetXaxis()->SetTitleOffset(1.2);
  hist2_shifted->GetYaxis()->SetTitleOffset(1.5);

  hist2_shifted->GetXaxis()->SetLabelSize(0.04);
  hist2_shifted->GetYaxis()->SetLabelSize(0.04);

  hist2_shifted->GetXaxis()->SetLabelOffset(0.01);
  hist2_shifted->GetYaxis()->SetLabelOffset(0.01);

  hist5->GetXaxis()->SetTitleSize(0.044);
  hist5->GetYaxis()->SetTitleSize(0.044);
  hist5->GetXaxis()->SetTitleOffset(1.2);
  hist5->GetYaxis()->SetTitleOffset(1.5);
  hist5->GetXaxis()->SetLabelSize(0.04);
  hist5->GetYaxis()->SetLabelSize(0.04);
  hist5->GetXaxis()->SetLabelOffset(0.01);
  hist5->GetYaxis()->SetLabelOffset(0.01);
    
  hist2_shifted->Draw("HIST");
  //hist2->Draw("HIST");
  hist1->Draw("HIST SAME");
  //hist3->Draw("HIST SAME");
  hist2_shifted->Draw("HIST SAME");

  //hist4->Draw("HIST SAME");

  hist2_shifted->SetTitle("Excitation energy spectrum;E_{x} (MeV);Counts / 70 keV");

  gStyle->SetOptStat(000);

  std::cout << "number of bins: " << hist1->GetNbinsX() << " bins/range = " << 14000. / hist1->GetNbinsX() <<  std::endl;

  double x_line = 5.6071; // 5.6 MeV
  double y_min = gPad->GetUymin();
  double y_max = gPad->GetUymax();

  double x_line_alpha = 4.415; // 5.6 MeV

  // Tworzenie linii
  TLine *vline = new TLine(x_line, y_min, x_line, y_max);
  vline->SetLineColor(kBlack);
  vline->SetLineStyle(2); // kreskowana
  vline->SetLineWidth(1);
  vline->Draw();

  TLine *vline2 = new TLine(x_line_alpha, y_min, x_line_alpha, y_max);
  vline2->SetLineColorAlpha(kBlack, 0.3);
  vline2->SetLineStyle(2); // kreskowana
  vline2->SetLineWidth(1);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextFont(42);          // 42 = Helvetica, 62 = bold Helvetica
  latex.SetTextColor(kRed);
  latex.SetTextColorAlpha(kGray+2, 0.4);
  //latex.DrawLatex(0.7, 0.5, "PRELIMINARY");

  leg->AddEntry(vline, "S_{p} of ^{18}F = 5.61 MeV", "l"); // "l" = linia
  //leg->AddEntry(vline2, "S_{#alpha} of ^{18}F = 4.415 MeV", "l"); // "l" = linia

  leg->Draw();

  canvas->cd(2);

  hist4->Scale(7);

  hist5->Draw("HIST");
  hist3->Draw("HIST SAME");
  hist4->Draw("HIST SAME");

  vline->Draw();
  vline2->Draw();

  hist5->SetTitle(";E_{x} (MeV);Counts / 70 keV");
  hist5->GetYaxis()->SetRangeUser(0, 1870);

  TLegend *leg2 = new TLegend(0.55, 0.7, 0.92, 0.898);
  leg2->AddEntry(hist5, "^{17}F(d,p) with ^{18}F recoils", "l");
  leg2->AddEntry(hist3, "^{17}F(d,p) with ^{17}O recoils", "l");
  leg2->AddEntry(hist4, "^{17}F(d,p) with ^{14}N recoils x 7", "l");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextSize(0.038);

  leg2->AddEntry(vline, "S_{p} of ^{18}F = 5.61 MeV", "l"); // "l" = linia
  leg2->AddEntry(vline2, "-Q_{#alpha} of ^{18}F = 4.415 MeV", "l"); // "l" = linia

  leg2->Draw();

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.18); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.01);  // wąski prawy margines
  gPad->SetTopMargin(0.0);
  

  canvas->SaveAs("ex_F_O.png");
  canvas->SaveAs("ex_F_O.pdf");
}
