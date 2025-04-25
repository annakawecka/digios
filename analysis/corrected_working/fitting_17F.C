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
    //fitFunc->FixParameter(3 * j + 3, peaks[j]);
    fitFunc->SetParameter(3 * j + 4, sigma);

    fitFunc->SetParLimits(3 * j + 2, 0, 1000);
    fitFunc->SetParLimits(3 * j + 3, peaks[j] - 0.04, peaks[j] + 0.04);
    fitFunc->SetParLimits(3 * j + 4, 0.05, sigma + 0.02);

    fitFunc->SetParName(3 * j + 2, ("Peak" + std::to_string(j + 1) + "_Amplitude").c_str());
    fitFunc->SetParName(3 * j + 3, ("Peak" + std::to_string(j + 1) + "_Position").c_str());
    fitFunc->SetParName(3 * j + 4, ("Peak" + std::to_string(j + 1) + "_Sigma").c_str());
  }

  fitFunc->SetNpx(1000);
  hist->Fit(fitFunc, "R");

  TCanvas *canvas = new TCanvas(("canvas_" + histName + outputSuffix).c_str(), histName.c_str(), 1600, 1200);
  canvas->Divide(1, 1);

  canvas->cd(1);
  hist->SetStats(false);
  
  hist->Draw("E1");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(kRed);
  hist->SetMarkerSize(1);
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
    peakFunc->SetLineColor(kBlue);
    peakFunc->Draw("SAME");
  }

  /* canvas->cd(2); */
  /* TPaveText *stats = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC"); */
  /* stats->SetTextSize(0.04); */
  /* stats->AddText(Form("Background Intercept: %.3f", fitFunc->GetParameter(0))); */
  /* stats->AddText(Form("Background Slope: %.3f", fitFunc->GetParameter(1))); */

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

    /* stats->AddText(Form("Peak %d Amplitude: %.3f", j + 1, amp)); */
    /* stats->AddText(Form("Peak %d Position: %.3f", j + 1, pos)); */
    /* stats->AddText(Form("Peak %d Sigma: %.3f", j + 1, sig)); */
    /* stats->AddText(Form("Peak %d Integral: %.3f", j + 1, integral)); */

    outFile << "Peak " << j + 1 << " Amplitude: " << amp << "  (" << amp_err << ")" << "\n";
    outFile << "Peak " << j + 1 << " Position: " << pos << "  (" << pos_err << ")"<< "\n";
    outFile << "Peak " << j + 1 << " Sigma: " << sig << "  (" << sig_err << ")"<< "\n";
    outFile << "Peak " << j + 1 << " Integral: " << integral << "  (" << integral_err << ")"<< "\n";
  }

  outFile << "\n";
  outFile.close();
  /* stats->Draw(); */
  canvas->Update();
  canvas->SaveAs(("plots_17F/fitting/" + histName + outputSuffix + ".png").c_str());
  canvas->SaveAs(("plots_17F/fitting/" + histName + outputSuffix + ".root").c_str());
}

void fitting_17F() {
    TFile *file = TFile::Open("rings_17F.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::vector<std::string> histNames = {"Ex_d0", "Ex_d1", "Ex_d2", "Ex_d3", "Ex_d4", "Ex_d5"};

    std::vector<std::vector<double>> peakSets[] = {
        {
	  {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335},
	  {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258},
	  {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
	  {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
	  {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
	  {     0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976}
        }
    };

    std::string suffixes[] = {""};
    std::string outputFiles[] = {
        "plots_17F/fitting/fit_parameters_17F.txt"
    };
    std::ofstream outFile1(outputFiles[0], std::ios::out | std::ios::trunc);

    for (int k = 0; k < 1; ++k) {
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
