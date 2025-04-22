#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

void performFitAndPlot(const std::string &histName, TH1 *hist, const std::vector<double> &peaks,
                       const std::string &outputSuffix, const std::string &outputFileName, double sigma = 0.08) {
    
  int nPeaks = peaks.size();
  std::string funcExpr = "pol1";
  for (int j = 0; j < nPeaks; ++j) funcExpr += "+gaus(" + std::to_string(3 * j + 2) + ")";

  TF1 *fitFunc = new TF1(("fitFunc_" + histName).c_str(), funcExpr.c_str(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  fitFunc->SetParameters(0, 0.0, 0.0); // Background intercept & slope
  fitFunc->SetParName(0, "Background_Intercept");
  fitFunc->SetParName(1, "Background_Slope");

  for (int j = 0; j < nPeaks; ++j) {
    fitFunc->SetParameter(3 * j + 2, hist->GetMaximum() * 0.5);
    fitFunc->SetParameter(3 * j + 3, peaks[j]);
    fitFunc->SetParameter(3 * j + 4, sigma);

    fitFunc->SetParLimits(3 * j + 3, peaks[j] - 0.1, peaks[j] + 0.1);
    fitFunc->SetParLimits(3 * j + 4, 0.05, sigma + 0.01);

    fitFunc->SetParName(3 * j + 2, ("Peak" + std::to_string(j + 1) + "_Amplitude").c_str());
    fitFunc->SetParName(3 * j + 3, ("Peak" + std::to_string(j + 1) + "_Position").c_str());
    fitFunc->SetParName(3 * j + 4, ("Peak" + std::to_string(j + 1) + "_Sigma").c_str());
  }

  fitFunc->SetNpx(1000);
  hist->Fit(fitFunc, "R");

  TCanvas *canvas = new TCanvas(("canvas_" + histName + outputSuffix).c_str(), histName.c_str(), 1600, 800);
  canvas->Divide(2, 1);

  canvas->cd(1);
  hist->SetStats(false);
  hist->Draw();
  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineStyle(2);
  fitFunc->Draw("SAME");

  for (int j = 0; j < nPeaks; ++j) {
    TF1 *peakFunc = new TF1(("peak_" + std::to_string(j)).c_str(),
			    "gaus(0)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    peakFunc->SetParameters(fitFunc->GetParameter(3 * j + 2),
			    fitFunc->GetParameter(3 * j + 3),
			    fitFunc->GetParameter(3 * j + 4));
    peakFunc->SetNpx(1000);
    peakFunc->SetLineColor(4);
    peakFunc->Draw("SAME");
  }

  canvas->cd(2);
  TPaveText *stats = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  stats->SetTextSize(0.04);
  stats->AddText(Form("Background Intercept: %.3f", fitFunc->GetParameter(0)));
  stats->AddText(Form("Background Slope: %.3f", fitFunc->GetParameter(1)));

  std::ofstream outFile(outputFileName, std::ios::app);
  outFile << "Histogram: " << histName << "\n";
  outFile << "Background Intercept: " << fitFunc->GetParameter(0) << "\n";
  outFile << "Background Slope: " << fitFunc->GetParameter(1) << "\n";

  for (int j = 0; j < nPeaks; ++j) {
    double amp = fitFunc->GetParameter(3 * j + 2);
    double pos = fitFunc->GetParameter(3 * j + 3);
    double sig = fitFunc->GetParameter(3 * j + 4);
    double integral = amp * TMath::Sqrt(2 * TMath::Pi()) * sig;

    stats->AddText(Form("Peak %d Amplitude: %.3f", j + 1, amp));
    stats->AddText(Form("Peak %d Position: %.3f", j + 1, pos));
    stats->AddText(Form("Peak %d Sigma: %.3f", j + 1, sig));
    stats->AddText(Form("Peak %d Integral: %.3f", j + 1, integral));

    outFile << "Peak " << j + 1 << " Amplitude: " << amp << "\n";
    outFile << "Peak " << j + 1 << " Position: " << pos << "\n";
    outFile << "Peak " << j + 1 << " Sigma: " << sig << "\n";
    outFile << "Peak " << j + 1 << " Integral: " << integral << "\n";
  }

  outFile << "\n";
  outFile.close();
  stats->Draw();
  canvas->Update();
  canvas->SaveAs(("plots_17O/fitting/" + histName + outputSuffix + ".png").c_str());
}

void fitting_17O() {
    TFile *file = TFile::Open("rings_17O.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::vector<std::string> histNames = {"Ex_d0", "Ex_d1", "Ex_d2", "Ex_d3", "Ex_d4", "Ex_d5"};

    std::vector<std::vector<double>> peakSets[] = {
        {
            {0.0},
            {0.0, 1.982},
            {1.982, 3.552, 3.63, 3.92},
            {1.982, 3.552, 3.63, 3.92, 5.255},
            {1.982, 3.552, 3.63, 3.92, 5.255, 6.2, 7.11},
            {1.982, 3.552, 3.63, 3.92, 5.255, 6.2, 7.11}
        },
        {
            {0.0},
            {0.0, 1.982},
            {1.982, 3.55, 3.92},
            {1.982, 3.55, 3.92, 5.255},
            {1.982, 3.55, 3.92, 5.255, 6.2, 7.11},
            {1.982, 3.55, 3.92, 5.255, 6.2, 7.11}
        }
    };

    std::string suffixes[] = {"", "_mergedPeaks"};
    std::string outputFiles[] = {
        "plots_17O/fitting/fit_parameters_17O.txt",
        "plots_17O/fitting/fit_parameters_17O_mergedPeaks.txt"
    };
    std::ofstream outFile1(outputFiles[0], std::ios::out | std::ios::trunc);
    std::ofstream outFile2(outputFiles[1], std::ios::out | std::ios::trunc);

    for (int k = 0; k < 2; ++k) {
        for (size_t i = 0; i < histNames.size(); ++i) {
            TH1 *hist = dynamic_cast<TH1*>(file->Get(histNames[i].c_str()));
            if (!hist) {
                std::cerr << "Histogram " << histNames[i] << " not found!" << std::endl;
                continue;
            }
            performFitAndPlot(histNames[i], hist, peakSets[k][i], suffixes[k], outputFiles[k]);
        }
    }

    file->Close();
}
