#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

void fitting_170() {
    TFile *file = TFile::Open("rings_17O.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::vector<std::string> histNames = {"Ex_d0", "Ex_d1", "Ex_d2", "Ex_d3", "Ex_d4", "Ex_d5"};

    std::vector<std::vector<double>> peakPositions = {
        {0.15},
        {0.15, 1.982 + 0.15, 3.44},
        {1.982 + 0.15, 3.55 + 0.15, 3.63 + 0.15, 3.92 + 0.15},
        {1.982 + 0.15, 3.55 + 0.15, 3.63 + 0.15, 3.92 + 0.15, 5.255 + 0.15},
        {1.982 + 0.15, 3.55 + 0.15, 3.63 + 0.15, 3.92 + 0.15, 5.255 + 0.15, 6.4, 7.20},
        {1.982 + 0.15, 3.55 + 0.15, 3.63 + 0.15, 3.92 + 0.15, 5.255 + 0.15}
    };

    const double sigma = 0.09;

    for (size_t i = 0; i < histNames.size(); ++i) {
        TH1 *hist = dynamic_cast<TH1*>(file->Get(histNames[i].c_str()));
        if (!hist) {
            std::cerr << "Histogram " << histNames[i] << " not found!" << std::endl;
            continue;
        }

        int nPeaks = peakPositions[i].size();
        std::string funcExpr = "pol1"; // First-degree polynomial
        for (int j = 0; j < nPeaks; ++j) {
            funcExpr += "+gaus(" + std::to_string(3 * j + 2) + ")"; // Gaussians
        }

        TF1 *fitFunc = new TF1(("fitFunc_" + histNames[i]).c_str(), funcExpr.c_str(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

        fitFunc->SetParameter(0, 0.0); // Intercept
        fitFunc->SetParameter(1, 0.0); // Slope

	fitFunc->SetParName(0, "Background_Intercept");
        fitFunc->SetParName(1, "Background_Slope");

        for (int j = 0; j < nPeaks; ++j) {
            fitFunc->SetParameter(3 * j + 2, hist->GetMaximum() * 0.5); // Amplitude
            fitFunc->SetParameter(3 * j + 3, peakPositions[i][j]);     // Mean
            fitFunc->SetParameter(3 * j + 4, sigma);                   // Sigma

            fitFunc->SetParLimits(3 * j + 3, peakPositions[i][j] - 0.1, peakPositions[i][j] + 0.1);
            fitFunc->SetParLimits(3 * j + 4, 0.04, sigma);

	    fitFunc->SetParName(3 * j + 2, ("Peak" + std::to_string(j + 1) + "_Amplitude").c_str());
            fitFunc->SetParName(3 * j + 3, ("Peak" + std::to_string(j + 1) + "_Position").c_str());
            fitFunc->SetParName(3 * j + 4, ("Peak" + std::to_string(j + 1) + "_Sigma").c_str());
        }

	fitFunc->SetNpx(1000);

        hist->Fit(fitFunc, "R");

	fitFunc->SetLineColor(kRed);
	fitFunc->SetLineStyle(2);

        TCanvas *canvas = new TCanvas(("canvas_" + histNames[i]).c_str(), histNames[i].c_str(), 1600, 800);
	canvas->Divide(2, 1);

	canvas->cd(1);
	hist->SetStats(kFALSE);
        hist->Draw();
	fitFunc->Draw("SAME");
	canvas->Update();

        canvas->cd(2);
        TPaveText *stats = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
        stats->SetTextSize(0.04);

	stats->AddText(Form("Background Intercept: %.3f", fitFunc->GetParameter(0)));
        stats->AddText(Form("Background Slope: %.3f", fitFunc->GetParameter(1)));

        for (int j = 0; j < nPeaks; ++j) {
            double amplitude = fitFunc->GetParameter(3 * j + 2);
            double peakPosition = fitFunc->GetParameter(3 * j + 3);
            double sigma = fitFunc->GetParameter(3 * j + 4);

            double integral = amplitude * TMath::Sqrt(2 * TMath::Pi()) * sigma;

            stats->AddText(Form("Peak %d Amplitude: %.3f", j + 1, amplitude));
            stats->AddText(Form("Peak %d Position: %.3f", j + 1, peakPosition));
            stats->AddText(Form("Peak %d Sigma: %.3f", j + 1, sigma));
            stats->AddText(Form("Peak %d Integral: %.3f", j + 1, integral));
        }

        stats->Draw();
        canvas->Update();

        canvas->SaveAs((histNames[i] + "_fitting.png").c_str());

	std::ofstream outFile("17O_analysis/fit_parameters_17O.txt", std::ios::app);  // Open file in append mode
	if (outFile.is_open()) {
	  outFile << "Histogram: " << histNames[i] << std::endl;
	  outFile << "Background Intercept: " << fitFunc->GetParameter(0) << std::endl;
	  outFile << "Background Slope: " << fitFunc->GetParameter(1) << std::endl;

	  for (int j = 0; j < nPeaks; ++j) {
	    double amplitude = fitFunc->GetParameter(3 * j + 2);  // Peak amplitude
	    double peakPosition = fitFunc->GetParameter(3 * j + 3); // Peak position
	    double sigma = fitFunc->GetParameter(3 * j + 4);      // Peak sigma
	    double integral = amplitude * TMath::Sqrt(2 * TMath::Pi()) * sigma;  // Integral

	    outFile << "Peak " << (j + 1) << " Amplitude: " << amplitude << std::endl;
	    outFile << "Peak " << (j + 1) << " Position: " << peakPosition << std::endl;
	    outFile << "Peak " << (j + 1) << " Sigma: " << sigma << std::endl;
	    outFile << "Peak " << (j + 1) << " Integral: " << integral << std::endl;
	  }
	  outFile << std::endl;  // Separate histograms with a blank line
	  outFile.close();
	} else {
	  std::cerr << "Unable to open file for writing!" << std::endl;
	}
    }

    //file->Close();
    //delete file;
}
