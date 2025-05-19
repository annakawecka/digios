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
    fitFunc->SetParLimits(3 * j + 4, 0.05, sigma + 0.02);

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
    double amp_err = fitFunc->GetParError(3 * j + 2);
    double pos = fitFunc->GetParameter(3 * j + 3);
    double pos_err = fitFunc->GetParError(3 * j + 3);
    double sig = fitFunc->GetParameter(3 * j + 4);
    double sig_err = fitFunc->GetParError(3 * j + 4);
    double integral = amp * TMath::Sqrt(2 * TMath::Pi()) * sig;
    double integral_err = sqrt(
				 pow(sig * sqrt(2 * TMath::Pi()) * amp_err, 2) +
				 pow(amp * sqrt(2 * TMath::Pi()) * sig_err, 2)
				 );

    stats->AddText(Form("Peak %d Amplitude: %.3f", j + 1, amp));
    stats->AddText(Form("Peak %d Position: %.3f", j + 1, pos));
    stats->AddText(Form("Peak %d Sigma: %.3f", j + 1, sig));
    stats->AddText(Form("Peak %d Integral: %.3f", j + 1, integral));

    outFile << "Peak " << j + 1 << " Amplitude: " << amp << "  (" << amp_err << ")" << "\n";
    outFile << "Peak " << j + 1 << " Position: " << pos << "  (" << pos_err << ")"<< "\n";
    outFile << "Peak " << j + 1 << " Sigma: " << sig << "  (" << sig_err << ")"<< "\n";
    outFile << "Peak " << j + 1 << " Integral: " << integral << "  (" << integral_err << ")"<< "\n";
  }

  outFile << "\n";
  outFile.close();
  stats->Draw();
  canvas->Update();
  canvas->SaveAs(("plots_17O/fitting_half_dets/" + histName + outputSuffix + ".png").c_str());
}

void fitting_17O_half_dets() {
    TFile *file = TFile::Open("rings_17O_half_dets.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::vector<std::string> histNames = {"Ex_d0_half_dets", "Ex_d1_half_dets", "Ex_d2_half_dets", "Ex_d3_half_dets", "Ex_d4_half_dets", "Ex_d5_half_dets",
					  "Ex_d6_half_dets", "Ex_d7_half_dets", "Ex_d8_half_dets", "Ex_d9_half_dets", "Ex_d10_half_dets", "Ex_d11_half_dets"};

    std::vector<std::vector<double>> peakSets[] = {
        {
            {0.0},
	    {0.0},
            {0.0, 1.982},
	    {0.0, 1.982},
            {0.0, 1.982, 3.552, 3.63, 3.92},
	    {0.0, 1.982, 3.552, 3.63, 3.92},
            {0.0, 1.982, 3.552, 3.63, 3.92, 5.255},
	    {0.0, 1.982, 3.552, 3.63, 3.92, 5.255},
            {1.982, 3.552, 3.63, 3.92, 5.255, 6.2, 7.11},
	    {1.982, 3.552, 3.63, 3.92, 5.255, 6.2, 7.11},
	    {1.982, 3.552, 3.63, 3.92, 5.255, 6.2, 7.11},
            {1.982, 3.552, 3.63, 3.92, 5.255, 6.2, 7.11}
        },
        {
            {0.0},
	    {0.0},
            {0.0, 1.982},
	    {0.0, 1.982},
            {0.0, 1.982, 3.55, 3.92},
	    {0.0, 1.982, 3.55, 3.92},
            {0.0, 1.982, 3.55, 3.92, 5.255},
	    {0.0, 1.982, 3.55, 3.92, 5.255},
            {1.982, 3.55, 3.92, 5.255, 6.2, 7.11},
	    {1.982, 3.55, 3.92, 5.255, 6.2, 7.11},
	    {1.982, 3.55, 3.92, 5.255, 6.2, 7.11},
            {1.982, 3.55, 3.92, 5.255, 6.2, 7.11}
        }
    };

    std::string suffixes[] = {"", "_mergedPeaks"};
    std::string outputFiles[] = {
        "plots_17O/fitting_half_dets/fit_parameters_17O_half_det.txt",
        "plots_17O/fitting_half_dets/fit_parameters_17O_mergedPeaks_half_det.txt"
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
