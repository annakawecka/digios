#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

void performFitAndPlot(const std::string &histName, TH1 *hist, const std::vector<double> &peaks,
                       const std::string &outputSuffix, const std::string &outputFileName,
                       double sigma = 0.08, bool fixDistances = false) {

    int nPeaks = peaks.size();
    std::string funcExpr = "pol1";
    std::vector<double> relDistances;

    if (fixDistances && nPeaks > 0) {
        for (int j = 1; j < nPeaks; ++j)
            relDistances.push_back(peaks[j] - peaks[0]);
        funcExpr += "+gaus(2)";
        for (int j = 1; j < nPeaks; ++j)
            funcExpr += "+gaus(" + std::to_string(3 * j + 2) + ")";
    } else {
        for (int j = 0; j < nPeaks; ++j)
            funcExpr += "+gaus(" + std::to_string(3 * j + 2) + ")";
    }

    TF1 *fitFunc = new TF1(("fitFunc_" + histName).c_str(), funcExpr.c_str(),
                           hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fitFunc->SetParameters(0, 0.0, 0.0);
    fitFunc->SetParName(0, "Background_Intercept");
    fitFunc->SetParName(1, "Background_Slope");

    if (fixDistances && nPeaks > 0) {
        fitFunc->SetParameter(2, hist->GetMaximum() * 0.5);
        fitFunc->SetParameter(3, peaks[0]);
        fitFunc->SetParameter(4, sigma);
        fitFunc->SetParLimits(2, 0, 1000);
        fitFunc->SetParLimits(3, peaks[0] - 0.03, peaks[0] + 0.03);
        fitFunc->SetParLimits(4, 0.05, 0.12);
        for (int j = 1; j < nPeaks; ++j) {
            int idx = 3 * j + 2;
            fitFunc->SetParameter(idx, hist->GetMaximum() * 0.5);
            //fitFunc->SetParameter(idx + 1, peaks[0] + relDistances[j - 1]);
            fitFunc->FixParameter(idx + 1, peaks[0] + relDistances[j - 1]);
            fitFunc->SetParameter(idx + 2, sigma);
            //fitFunc->FixParameter(idx + 2, sigma);
	    fitFunc->SetParLimits(idx + 2, 0.05, 0.12);
	    fitFunc->SetParLimits(idx, 0, 1000);
        }
    } else {
        for (int j = 0; j < nPeaks; ++j) {
            fitFunc->SetParameter(3 * j + 2, hist->GetMaximum() * 0.5);
            fitFunc->SetParameter(3 * j + 3, peaks[j]);
            fitFunc->SetParameter(3 * j + 4, sigma);
            fitFunc->SetParLimits(3 * j + 2, 0, 1000);
            fitFunc->SetParLimits(3 * j + 3, peaks[j] - 0.03, peaks[j] + 0.03);
            fitFunc->SetParLimits(3 * j + 4, 0.05, 0.12);
        }
    }

    fitFunc->SetNpx(2000);
    hist->Fit(fitFunc, "R");

    std::string plotDir = "plots_17F/fitting/" + outputSuffix;
    std::filesystem::create_directories(plotDir);
    TCanvas *canvas = new TCanvas(("canvas_" + histName + outputSuffix).c_str(), histName.c_str(), 1600, 1200);
    canvas->cd();
    hist->SetStats(false);
    
    hist->Draw("E1");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kRed);
    hist->SetMarkerSize(1);
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineStyle(2);
    fitFunc->Draw("SAME");

    for (int j = 0; j < nPeaks; ++j) {
        TF1 *peakFunc = new TF1(("peak_" + std::to_string(j)).c_str(), "gaus(0)",
                                hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        peakFunc->SetParameters(fitFunc->GetParameter(3 * j + 2),
                                fitFunc->GetParameter(3 * j + 3),
                                fitFunc->GetParameter(3 * j + 4));
        peakFunc->SetLineColor(kBlue);
        peakFunc->Draw("SAME");
    }

    canvas->SaveAs((plotDir + "/" + histName + ".png").c_str());
    canvas->SaveAs((plotDir + "/" + histName + ".root").c_str());

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
        double integral_err = sqrt(pow(sig * TMath::Sqrt(2 * TMath::Pi()) * amp_err, 2) +
                                   pow(amp * TMath::Sqrt(2 * TMath::Pi()) * sig_err, 2));
        outFile << "Peak " << j + 1 << " Amplitude: " << amp << " (" << amp_err << ")\n";
        outFile << "Peak " << j + 1 << " Position: " << pos << " (" << pos_err << ")\n";
        outFile << "Peak " << j + 1 << " Sigma: " << sig << " (" << sig_err << ")\n";
        outFile << "Peak " << j + 1 << " Integral: " << integral << " (" << integral_err << ")\n";
    }
    outFile << "\n";
    outFile.close();
}

void performFitAndPlotOld(const std::string &histName, TH1 *hist, const std::vector<double> &peaks,
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
    //fitFunc->FixParameter(3 * j + 4, sigma);

    fitFunc->SetParLimits(3 * j + 2, 0, 1000);
    fitFunc->SetParLimits(3 * j + 3, peaks[j] - 0.03, peaks[j] + 0.03);
    fitFunc->SetParLimits(3 * j + 4, 0.05, 0.12);

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
        },
	{
	  {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335},
	  {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159},
	  {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
	  {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
	  {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
	  {     0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976}
        }
    };

    //  0. (1+), 0.9372 (3+), 1.04155 (0+), 1.08054 (0-), 1.12136 (5+), 1.4, 1.6, 1.9, 1.70081 (1+), 2.10061 (2-), 2.52335 (2+), 3.06184 (2+), 3.13387 (1-), 3.3582 (3+), 3.72419 (1+), 3.79149 (3-), 3.83917 (2+), 4.1159 (3+), 4.2258 (2-), 4.36015 (1+), 4.3981 (4-), 4.652 (4+), 4.753 (0+), 4.8483 (5-), 4.860 (1-), 4.9636 (2+), 5.2976 (4+)

    double sigma_low = (0.0815 + 0.0903 + 0.0729 + 0.0843 + 0.0926) / 5.;
    double sigma_high = (0.124 + 0.097 + 0.117) / 3.;
    double sigma_mean = (sigma_low + sigma_high) / 2.;
    std::vector<double> sigmas = {sigma_low, sigma_mean, sigma_high};

    std::string suffixes[] = {"", "_positive"};
    std::string outputFiles[] = {
      "plots_17F/fitting/fit_parameters_17F.txt",
      "plots_17F/fitting/fit_parameters_17F_positive.txt"
    };
    std::ofstream outFile1(outputFiles[0], std::ios::out | std::ios::trunc);

    /*for (int k = 0; k < 2; ++k) {
        for (size_t i = 0; i < histNames.size(); ++i) {
            TH1 *hist = dynamic_cast<TH1*>(file->Get(histNames[i].c_str()));
            if (!hist) {
                std::cerr << "Histogram " << histNames[i] << " not found!" << std::endl;
                continue;
            }
            performFitAndPlot(histNames[i], hist, peakSets[k][i], suffixes[k], outputFiles[k], sigma_mean);
        }
	}*/

    for (int k = 0; k < 2; ++k) { // standard vs. "_positive"
        for (double sigma : sigmas) {
            for (bool fixDist : {false, true}) {
                std::string sigmaLabel = (std::abs(sigma - sigma_low) < 1e-5)   ? "_sigmaLow" :
                                         (std::abs(sigma - sigma_high) < 1e-5)  ? "_sigmaHigh" :
                                                                                 "_sigmaMean";
                std::string distLabel = fixDist ? "_fixedDist" : "_freeDist";
                std::string fullSuffix = suffixes[k] + sigmaLabel + distLabel;
                std::string fullOutputFile = "plots_17F/fitting/fit_parameters_17F" + fullSuffix + ".txt";

                std::ofstream outFile(fullOutputFile, std::ios::out | std::ios::trunc);

                for (size_t i = 0; i < histNames.size(); ++i) {
                    TH1 *hist = dynamic_cast<TH1*>(file->Get(histNames[i].c_str()));
                    if (!hist) {
                        std::cerr << "Histogram " << histNames[i] << " not found!" << std::endl;
                        continue;
                    }
                    performFitAndPlot(histNames[i], hist, peakSets[k][i], fullSuffix, fullOutputFile, sigma, fixDist);
                }
            }
        }
    }

    file->Close();
}
