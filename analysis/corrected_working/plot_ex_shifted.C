#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>

void plot_ex_shifted() {
    TFile* file1 = TFile::Open("plots_17F/Ex_x_recoil_coinTime_gated.root");
    TFile* file2 = TFile::Open("plots_17O/Ex_x_recoil_coinTime_gated.root");

    if (!file1 || !file2) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }

    TH1* hist1_orig = dynamic_cast<TH1*>(file1->Get("x_rdt_coinTime_gatedEx"));
    TH1* hist2_orig = dynamic_cast<TH1*>(file2->Get("x_rdt_coinTime_gatedEx"));

    if (!hist1_orig || !hist2_orig) {
        std::cerr << "Error retrieving histograms" << std::endl;
        return;
    }

    TH1* hist1 = (TH1*)hist1_orig->Clone("hist1");
    hist1->SetDirectory(0);
    TH1* hist2 = (TH1*)hist2_orig->Clone("hist2");
    hist2->SetDirectory(0);

    file1->Close();
    file2->Close();

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

    TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);
    hist1->SetLineColor(kRed); 		 // 17F
    hist2_shifted->SetLineColor(kBlack); // 17O
    hist2->SetLineColor(kBlue);

    std::cout << hist1->Integral() << std::endl << hist2->Integral()  << std::endl << hist2_shifted->Integral() << std::endl;
    
    hist2_shifted->Draw("HIST");
    //hist2->Draw("HIST");
    hist1->Draw("HIST SAME");

    canvas->BuildLegend();
}
