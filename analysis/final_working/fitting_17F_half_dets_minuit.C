#include <TMinuit.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

// ---------- Globals used by the FCN ----------
static TH1 *gHist = nullptr;
static std::vector<double> gDeltas;   // relative offsets: peak_i = pos0 + gDeltas[i]
static int gNpeaks = 0;
static double gFitMin = -0.5;
static double gFitMax = 8.0;

// ---------- Utility: determine effective fit range from histogram ----------
void getEffectiveRange(TH1* hist, double &xmin, double &xmax, double margin = 0.1, double minWidth = 0.5) {
    int nbins = hist->GetNbinsX();
    int firstNonEmpty = 1;
    int lastNonEmpty = nbins;

    bool foundFirst = false;
    for (int i = 1; i <= nbins; ++i) {
        if (hist->GetBinContent(i) > 0) { firstNonEmpty = i; foundFirst = true; break; }
    }
    if (!foundFirst) {
        // no data -> fallback to full range
        xmin = hist->GetXaxis()->GetXmin();
        xmax = hist->GetXaxis()->GetXmax();
        return;
    }
    for (int i = nbins; i >= 1; --i) {
        if (hist->GetBinContent(i) > 0) { lastNonEmpty = i; break; }
    }

    xmin = hist->GetBinLowEdge(firstNonEmpty) - margin;
    xmax = hist->GetBinLowEdge(lastNonEmpty) + hist->GetBinWidth(lastNonEmpty) + margin;

    // enforce overall clamps
    if (xmin < -0.5) xmin = -0.5;
    if (xmax > 8.0)  xmax = 8.0;

    // ensure minimal width so minimizer has something to work with
    if ((xmax - xmin) < minWidth) {
        double center = 0.5 * (xmin + xmax);
        xmin = center - minWidth / 2.0;
        xmax = center + minWidth / 2.0;
        if (xmin < -0.5) xmin = -0.5;
        if (xmax > 8.0)  xmax = 8.0;
    }
}

// ---------- Automatic amplitude estimator -----------
// Uses a simple approach: take the maximum bin content near each expected peak
// and convert to amplitude for normalized Gaussian (area = height * sigma * sqrt(2*pi))
std::vector<double> estimateInitialAmplitudes(TH1* hist, const std::vector<double>& peakPositions, double sigma_guess = 0.08, int searchRadiusBins = 4) {
    std::vector<double> amps;
    double factor = sigma_guess * sqrt(2.0 * TMath::Pi()); // area = height * sigma * sqrt(2pi)
    for (double pos : peakPositions) {
        int bin = hist->FindBin(pos);
        double maxC = 0.0;
        for (int b = bin - searchRadiusBins; b <= bin + searchRadiusBins; ++b) {
            if (b < 1 || b > hist->GetNbinsX()) continue;
            double c = hist->GetBinContent(b);
            if (c > maxC) maxC = c;
        }
        double ampGuess = std::max(1.0, maxC * factor); // don't set extremely small initial amp
        amps.push_back(ampGuess);
    }
    return amps;
}

// ---------- Model: sum of normalized Gaussians (area = amplitude) + linear background ----------
// par layout for likelihoodFCN and TMinuit:
// par[0] = background intercept
// par[1] = background slope
// par[2] = sigma0
// par[3] = sigma1
// par[4] = pos0 (first peak position)
// par[5+i] = amplitude of peak i  (area under gaussian)
double modelFunction(double x, const double *par) {
    double bg = par[0] + par[1] * x;
    double sigma0 = par[2];
    double sigma1 = par[3];
    double pos0 = par[4];

    double val = bg;
    for (int i = 0; i < gNpeaks; ++i) {
        double amp = par[5 + i];
        double pos = pos0 + gDeltas[i];
        double sigma = sigma0 + sigma1 * (pos - pos0);
        if (sigma < 1e-3) sigma = 1e-3;
        // Use normalized Gaussian: TMath::Gaus(x, mean, sigma, kTRUE) => integral over x is 1
        val += amp * TMath::Gaus(x, pos, sigma, true);
    }
    return val;
}

// p[0] = background intercept
// p[1] = background slope
// p[2] = sigma0
// p[3] = sigma1
// p[4+i] = amplitude of peak i
// p[4+N+i] = position of peak i
double modelFreePositions(double x, double *par) {
  double val = par[0] + par[1]*x;
  double sigma0 = par[2];
  double sigma1 = par[3];

  for (int i = 0; i < gNpeaks; ++i) {
    double amp = par[4 + i];
    double pos = par[4 + gNpeaks + i];
    double sigma = sigma0 + sigma1 * pos;
    if (sigma < 1e-3) sigma = 1e-3;
    val += amp * TMath::Gaus(x, pos, sigma, true);
  }
  return val;
}

// ---------- Poisson log-likelihood FCN ----------
// fval = -2 * log L = 2 * sum( mu - n * log mu )  (drop ln(n!) const term)
void likelihoodFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t iflag) {
    // Set fit boundaries from global fit range
    int bin1 = gHist->FindBin(gFitMin);
    int bin2 = gHist->FindBin(gFitMax);
    double neg2logL = 0.0;

    for (int b = bin1; b <= bin2; ++b) {
        double x = gHist->GetBinCenter(b);
        double n = gHist->GetBinContent(b);
        double mu = modelFunction(x, par);
        if (mu <= 0) mu = 1e-12; // avoid log(0) and negative expectations

        if (n > 0) {
            neg2logL += 2.0 * (mu - n * std::log(mu));
        } else {
            // n == 0
            neg2logL += 2.0 * mu;
        }
    }

    fval = neg2logL;
}

void likelihoodFCNFreePos(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t iflag) {
    int bin1 = gHist->FindBin(gFitMin);
    int bin2 = gHist->FindBin(gFitMax);
    double neg2logL = 0.0;

    for (int b = bin1; b <= bin2; ++b) {
        double x = gHist->GetBinCenter(b);
        double n = gHist->GetBinContent(b);
        double mu = modelFreePositions(x, par);
        if (mu <= 0) mu = 1e-12;

        if (n > 0) {
            neg2logL += 2.0 * (mu - n * std::log(mu));
        } else {
            // n == 0
            neg2logL += 2.0 * mu;
        }
    }

    fval = neg2logL;
}

void fitHistogramLikelihoodTMinuitFreePos(TH1* hist,
					 const std::vector<double>& peakPositions,
					 TFile* outRootFile,
					 const std::string &outTxtFilename = "fit_results.txt",
					 bool savePlot = true) 
{

  if (!hist) {
    std::cerr << "Null histogram passed to fit routine.\n";
    return;
  }

  double fitMin, fitMax;
  getEffectiveRange(hist, fitMin, fitMax, 0.05, 0.4);
  gFitMin = fitMin;
  gFitMax = fitMax;
  
  gHist = hist;
  gNpeaks = peakPositions.size();

  std::vector<double> initAmps = estimateInitialAmplitudes(hist, peakPositions, 0.08, 2);

  int nPars = 4 + 2 * gNpeaks; // [bg0, bg1, sigma0, sigma1] + amplitudes + positions
  TMinuit minuit(nPars);
  //minuit.SetPrintLevel(1);
  minuit.SetFCN(likelihoodFCNFreePos);

  // Set initial parameter guesses
  double par, step, min, max;
  int ierrflg = 0;
    
  // Background intercept and slope
  minuit.mnparm(0, "BG_intercept",  hist->GetBinContent(hist->FindBin((fitMin+fitMax)/2.0)), 0.01, 0.0, 3., ierrflg);
  minuit.mnparm(1, "BG_slope",      0.0, 0.01, 0.0, 0.0, ierrflg);

  // sigma0, sigma1: choose reasonable limits so sigma near expected [0.075, 0.1]
  minuit.mnparm(2, "Sigma0", 0.085, 0.001, 0.07, 0.1, ierrflg);    // intercept
  minuit.mnparm(3, "Sigma1", 0.0,   0.001, 0.0, 0.05, ierrflg);  // slope (per unit x offset)


  // amplitudes
  for (int i = 0; i < gNpeaks; ++i) {
    double ainit = initAmps[i];
    double astep = std::max(1.0, 0.1 * ainit);
    double aupper = std::max(ainit * 20.0, ainit + 1000.0); // loose upper bound
    minuit.mnparm(4 + i, Form("Amp%d", i+1), ainit, astep, 0.0, aupper, ierrflg);
  }

  // peak positions (free)
  for (int i = 0; i < gNpeaks; ++i) {
    par = peakPositions[i];
    step = 0.001;
    min = par - 0.2;
    max = par + 0.2;
    minuit.DefineParameter(4 + gNpeaks + i, Form("pos_%d", i), par, step, min, max);
  }

  double arglist[10];
  arglist[0] = 10000;   // max calls
  arglist[1] = 1.0;    // tolerance (ERRDEF for likelihood ~1.0)
  minuit.mnexcm("MIGRAD", arglist, 2, ierrflg);

  // Recompute error matrix
  minuit.mnexcm("HESSE", arglist, 0, ierrflg);

  const char* fitStatusStr = minuit.fCstatu;
  int fitStatusCode = minuit.GetStatus();

  // --- Retrieve fitted parameters ---
  std::vector<double> parVal(nPars), parErr(nPars);
  for (int i = 0; i < nPars; ++i) {
    minuit.GetParameter(i, parVal[i], parErr[i]);
  }

  std::cout << "=== Fit summary for histogram: " << hist->GetName() << " ===\n";
  for (int i = 0; i < nPars; ++i) {
    TString parName;
    double val, err, bnd1, bnd2;
    int ivarbl;
    minuit.mnpout(i, parName, val, err, bnd1, bnd2, ivarbl);
    std::cout << Form("Param %2d : %12s = % .6g ± %.6g\n", i, parName.Data(), val, err);
  }

  // Write results to file
  std::ofstream out(outTxtFilename, std::ios::app);
  out << "Histogram: " << hist->GetName() << "\n";
  for (int i = 0; i < nPars; ++i) {
    TString parName;
    double val, err, bnd1, bnd2;
    int ivarbl;
    minuit.mnpout(i, parName, val, err, bnd1, bnd2, ivarbl);
    out << parName << " = " << val << " ± " << err << "\n";
  }
  out << "\n";
  out.close();

  // --- Draw fit and components ---
  TCanvas *c = new TCanvas(Form("c_%s", hist->GetName()), hist->GetName(), 1000, 700);
  gStyle->SetOptStat(0);
  hist->Draw("E");

  TF1 *fmodel = new TF1(Form("fmodel_%s", hist->GetName()), [&](double *x, double *p) {
      std::vector<double> par(nPars);
      for (int ii = 0; ii < nPars; ++ii) par[ii] = parVal[ii];
      return modelFreePositions(x[0], par.data());
    }, gFitMin, gFitMax, 0);
  fmodel->SetLineColor(kRed);
  fmodel->SetLineWidth(2);
  fmodel->SetNpx(500);
  fmodel->Draw("SAME");

  // Draw individual peaks with free positions
  for (int i = 0; i < gNpeaks; ++i) {
    double amp = parVal[4 + i];
    double pos = parVal[4 + gNpeaks + i]; // free peak position
    double sigma = parVal[2] + parVal[3] * pos;
    if (sigma < 1e-3) sigma = 1e-3;

    TF1 *fpeak = new TF1(Form("peak_%s_%d", hist->GetName(), i),
                         [](double *x, double *p) {
			   return p[0] * TMath::Gaus(x[0], p[1], p[2], true);
                         }, gFitMin, gFitMax, 3);

    fpeak->SetParameters(amp, pos, sigma);
    fpeak->SetLineColor(kBlue + i);
    fpeak->SetLineStyle(2);
    fpeak->SetNpx(500);
    fpeak->Draw("SAME");
  }

  // Draw background
  TF1 *fbkg = new TF1(Form("bkg_%s", hist->GetName()), [&](double *x, double*) {
      return parVal[0] + parVal[1] * x[0];
    }, gFitMin, gFitMax, 0);
  fbkg->SetLineColor(kGreen + 2);
  fbkg->SetLineStyle(3);
  fbkg->SetNpx(500);
  fbkg->Draw("SAME");

  // Show a small stats box
  TPaveText *pt = new TPaveText(0.55, 0.55, 0.9, 0.9, "NDC");
  pt->SetFillColor(0);
  pt->AddText(Form("Sigma0 = %.4g +- %.4g", parVal[2], parErr[2]));
  pt->AddText(Form("Sigma1 = %.4g +- %.4g", parVal[3], parErr[3]));
  for (int i = 0; i < gNpeaks; ++i) {
    pt->AddText(Form("Area%d = %.4g +- %.4g", i + 1, parVal[4 + i], parErr[4 + i]));
    pt->AddText(Form("Pos%d  = %.4g +- %.4g", i + 1, parVal[4 + gNpeaks + i], parErr[4 + gNpeaks + i]));
  }
  pt->Draw();

  if (savePlot) {
    c->SaveAs(Form("plots_17F/minuit_new/fit_%s_free.png", hist->GetName())); // PNGs
    if (outRootFile && outRootFile->IsOpen()) {
      outRootFile->cd();
      c->Write();
    }
  }

  std::ofstream outFile("plots_17F/minuit_new/fit_results_single_freePos.txt", std::ios::app);  // Append mode

  outFile << "Fit results for histogram: " << hist->GetName() << "\n";

  for (int i = 0; i < nPars; ++i) {
    TString parName;
    double val, err, b1, b2;
    int ivar;

    minuit.mnpout(i, parName, val, err, b1, b2, ivar);

    outFile << Form("%-15s = %.6g ± %.6g\n", parName.Data(), val, err);
  }

  outFile << "Fit status code: " << fitStatusCode << "\n";
  outFile << "Fit status message: " << fitStatusStr << "\n";

  outFile << "\n";

  for (int i = 0; i < peakPositions.size(); i++){
    TString parName;
    double val, err, b1, b2;
    int ivar;

    minuit.mnpout(4 + peakPositions.size() + i, parName, val, err, b1, b2, ivar);

    outFile << Form("%s - ref   = %.6g - %.6g = %.6g \n", parName.Data(), val, peakPositions[i], val - peakPositions[i]);
  }
  
  outFile << "-------------------------------------\n";

  outFile << "\n";

  outFile.close();
  
 
}


// ---------- Fit routine using TMinuit (Poisson likelihood) ----------
void fitHistogramLikelihoodTMinuit(TH1* hist,
                                   const std::vector<double>& peakPositions,
				   TFile* outRootFile,
                                   const std::string &outTxtFilename = "fit_results.txt",
                                   bool savePlot = true) {
  if (!hist) {
    std::cerr << "Null histogram passed to fit routine.\n";
    return;
  }

  // Compute fit range from histogram contents
  double fitMin, fitMax;
  getEffectiveRange(hist, fitMin, fitMax, 0.05, 0.4);
  gFitMin = fitMin;
  gFitMax = fitMax;

  // Setup globals
  gHist = hist;
  gNpeaks = peakPositions.size();
  gDeltas.resize(gNpeaks);
  for (int i = 0; i < gNpeaks; ++i) gDeltas[i] = peakPositions[i] - peakPositions[0];

  // Initial amplitude guesses (auto)
  std::vector<double> initAmps = estimateInitialAmplitudes(hist, peakPositions, 0.08, 2);

  // Number of parameters: bg0, bg1, sigma0, sigma1, pos0, amps(N) => 5 + N
  int nPars = 5 + gNpeaks;
  TMinuit minuit(nPars);
  minuit.SetFCN(likelihoodFCN);

  int ierflg = 0;

  // --- Define parameters (index, name, start, step, lower, upper, ierflg) ---
  // Background
  minuit.mnparm(0, "BG_intercept",  hist->GetBinContent(hist->FindBin((fitMin+fitMax)/2.0)), 0.01, 0., 5., ierflg);
  minuit.mnparm(1, "BG_slope",      0.0, 0.01, 0., 0.01, ierflg);

  // sigma0, sigma1: choose reasonable limits so sigma near expected [0.075, 0.1]
  minuit.mnparm(2, "Sigma0", 0.085, 0.001, 0.07, 0.1, ierflg);    // intercept, was :  0.085, 0.001, 0.07, 0.1, ierflg); 
  minuit.mnparm(3, "Sigma1", 0.0,   0.001, 0.0, 0.05, ierflg);  // slope (per unit x offset)

  // common group position (pos0)
  minuit.mnparm(4, "Pos0", peakPositions[0], 0.01, peakPositions[0] - 0.5, peakPositions[0] + 0.5, ierflg);

  // amplitudes (areas); enforce >= 0 and put an upper cap to avoid runaway
  for (int i = 0; i < gNpeaks; ++i) {
    double ainit = initAmps[i];
    double astep = std::max(1.0, 0.1 * ainit);
    double aupper = std::max(ainit * 20.0, ainit + 1000.0); // loose upper bound
    minuit.mnparm(5 + i, Form("Amp%d", i+1), ainit, astep, 0.0, aupper, ierflg);
  }

  // --- Run minimization ---
  double arglist[10];
  arglist[0] = 10000;   // max calls
  arglist[1] = 1.0;    // tolerance (ERRDEF for likelihood ~1.0)
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  // Recompute error matrix
  minuit.mnexcm("HESSE", arglist, 0, ierflg);

  const char* fitStatusStr = minuit.fCstatu;
  int fitStatusCode = minuit.GetStatus();

  // --- Retrieve fitted parameters ---
  std::vector<double> parVal(nPars), parErr(nPars);
  for (int i = 0; i < nPars; ++i) {
    minuit.GetParameter(i, parVal[i], parErr[i]);
  }

  std::cout << "=== Fit summary for histogram: " << hist->GetName() << " ===\n";
  for (int i = 0; i < nPars; ++i) {
    TString parName;
    double val, err, bnd1, bnd2;
    int ivarbl;
    minuit.mnpout(i, parName, val, err, bnd1, bnd2, ivarbl);
    std::cout << Form("Param %2d : %12s = % .6g ± %.6g\n", i, parName.Data(), val, err);
  }

  // Write results to file
  std::ofstream out(outTxtFilename, std::ios::app);
  out << "Histogram: " << hist->GetName() << "\n";
  for (int i = 0; i < nPars; ++i) {
    TString parName;
    double val, err, bnd1, bnd2;
    int ivarbl;
    minuit.mnpout(i, parName, val, err, bnd1, bnd2, ivarbl);
    out << parName << " = " << val << " ± " << err << "\n";
  }
  for (int i = 0; i < gNpeaks; ++i) out << Form("Pos%d = %.4g (%.4g)\n", i+1, parVal[4] + gDeltas[i], peakPositions[i]);
  out << "\n";
  out.close();

  // --- Draw fit and components ---
  TCanvas *c = new TCanvas(Form("c_%s", hist->GetName()), hist->GetName(), 1000, 700);
  gStyle->SetOptStat(0);
  hist->Draw("E");

  // Create TF1 to draw the full model (wraps modelFunction and fitted params)
  TF1 *fmodel = new TF1(Form("fmodel_%s", hist->GetName()), [&](double *x, double *p) {
      // build an array of parameters for modelFunction compatible layout
      std::vector<double> par(nPars);
      for (int ii = 0; ii < nPars; ++ii) par[ii] = parVal[ii];
      return modelFunction(x[0], par.data());
    }, gFitMin, gFitMax, 0);
  fmodel->SetLineColor(kRed);
  fmodel->SetLineWidth(2);
  fmodel->SetNpx(500);
  fmodel->Draw("SAME");

  // draw individual peaks
  for (int i = 0; i < gNpeaks; ++i) {
    double pos0 = parVal[4];
    double pos = pos0 + gDeltas[i];
    double sigma = parVal[2] + parVal[3] * (pos - pos0);
    if (sigma < 1e-3) sigma = 1e-3;
    double amp = parVal[5 + i];

    TF1 *fpeak = new TF1(Form("peak_%s_%d", hist->GetName(), i),
			 [](double *x, double *p) {
			   return p[0] * TMath::Gaus(x[0], p[1], p[2], true);
			 }, gFitMin, gFitMax, 3);

    fpeak->SetParameters(amp, pos, sigma);
    fpeak->SetLineColor(kBlue + i);
    fpeak->SetLineStyle(2);
    fpeak->SetNpx(500);
    fpeak->Draw("SAME");
  }

  // draw background
  TF1 *fbkg = new TF1(Form("bkg_%s", hist->GetName()), [&](double *x, double*) {
      return parVal[0] + parVal[1] * x[0];
    }, gFitMin, gFitMax, 0);
  fbkg->SetLineColor(kGreen+2);
  fbkg->SetLineStyle(3);
  fbkg->SetNpx(500);
  fbkg->Draw("SAME");

  // show a small stats box
  TPaveText *pt = new TPaveText(0.55, 0.55, 0.9, 0.9, "NDC");
  pt->SetFillColor(0);
  pt->AddText(Form("Sigma0 = %.4g +- %.4g", parVal[2], parErr[2]));
  pt->AddText(Form("Sigma1 = %.4g +- %.4g", parVal[3], parErr[3]));
  pt->AddText(Form("Pos0   = %.4g +- %.4g", parVal[4], parErr[4]));
  for (int i = 0; i < gNpeaks; ++i) pt->AddText(Form("Area%d = %.4g +- %.4g", i+1, parVal[5+i], parErr[5+i]));
  for (int i = 0; i < gNpeaks; ++i) pt->AddText(Form("Pos%d = %.4g", i+1, parVal[4] + gDeltas[i]));
  pt->Draw();

  if (savePlot) {
    c->SaveAs(Form("plots_17F/minuit_new/fit_%s.png", hist->GetName())); // PNGs
    //c->SaveAs(Form("plots_17F/minuit/fit_%s.png", hist->GetName())); // PNGs
    if (outRootFile && outRootFile->IsOpen()) {
      outRootFile->cd();
      c->Write();
    }
  }

  std::ofstream outFile("plots_17F/minuit_new/fit_results.txt", std::ios::app);  // Append mode
  //std::ofstream outFile("plots_17F/minuit/fit_results.txt", std::ios::app);  // Append mode

  outFile << "Fit results for histogram: " << hist->GetName() << "\n";

  for (int i = 0; i < nPars; ++i) {
    TString parName;
    double val, err, b1, b2;
    int ivar;

    minuit.mnpout(i, parName, val, err, b1, b2, ivar);

    outFile << Form("%-15s = %.6g ± %.6g\n", parName.Data(), val, err);
  }
  for (int i = 0; i < gNpeaks; ++i) outFile << Form("Pos%d = %.4g (%.4g)\n", i+1, parVal[4] + gDeltas[i], peakPositions[i]);

  outFile << "Fit status code: " << fitStatusCode << "\n";
  outFile << "Fit status message: " << fitStatusStr << "\n";
  outFile << "-------------------------------------\n";

  outFile << "\n";  // blank line between fits

  outFile.close();

}

// ---------- Example top-level function to demonstrate usage ----------
void fitting_17F_half_dets_minuit() {

  // fitting rings
  TFile *f = TFile::Open("rings_17F_half_dets.root", "READ");
  //TFile *f = TFile::Open("rings_17F_half_dets.root", "READ");
  if (!f || f->IsZombie()) { std::cerr << "Cannot open file\n"; return; }

  std::vector<std::string> histNames = {"Ex_d0_half_dets", "Ex_d1_half_dets", "Ex_d2_half_dets", "Ex_d3_half_dets", "Ex_d4_half_dets", "Ex_d5_half_dets",
					"Ex_d6_half_dets", "Ex_d7_half_dets", "Ex_d8_half_dets", "Ex_d9_half_dets", "Ex_d10_half_dets", "Ex_d11_half_dets"};

  // peak sets for each histogram (example)
  std::vector<std::vector<double>> peakSets = {
    /*{0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},
    {0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149, 3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753, 4.8483, 4.860, 4.9636, 5.2976},*/


    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335},
    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335},
    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159},
    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159},
    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {0.0, 0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {0.9372, 1.04155, 1.12136, 1.4, 1.6, 1.9, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {1.12136, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976},
    {1.12136, 1.70081, 2.52335, 3.06184, 3.3582, 3.72419, 3.83917, 4.1159, 4.36015, 4.652, 4.753, 4.9636, 5.2976}
  };

  // fitting individual detectors
  /*TFile *f = TFile::Open("rings_17F_single_dets.root", "READ");
  if (!f || f->IsZombie()) { std::cerr << "Cannot open file\n"; return; }

  std::vector<std::string> histNames = {"Ex_single_det0",
					"Ex_single_det1",
					"Ex_single_det2",
					"Ex_single_det3",
					"Ex_single_det4",
					"Ex_single_det5",
					"Ex_single_det6",
					"Ex_single_det7",
					"Ex_single_det8",
					"Ex_single_det9",
					"Ex_single_det10",
					"Ex_single_det11",
					"Ex_single_det12",
					"Ex_single_det13",
					"Ex_single_det14",
					"Ex_single_det15",
					"Ex_single_det16",
					"Ex_single_det17",
					"Ex_single_det18",
					"Ex_single_det19",
					"Ex_single_det20",
					"Ex_single_det21",
					"Ex_single_det22",
					"Ex_single_det23"
  };

  // peak sets for each histogram (example)
  std::vector<std::vector<double>> peakSets = {
    {0.0},
    {0.0, 1.98207},
    {1.98207, 3.55484, 3.63376, 3.92044},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {0.0}, // 6
    {0.0, 1.98207},
    {1.98207, 3.55484, 3.63376, 3.92044},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {0.0}, // 12
    {0.0, 1.98207},
    {1.98207, 3.55484, 3.63376, 3.92044},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {0.0},
    {0.0, 1.98207},
    {1.98207, 3.55484, 3.63376, 3.92044},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169},
    {1.98207, 3.55484, 3.63376, 3.92044, 5.2548},//, 6.19822, 7.1169}
    };*/

  TFile allFits("plots_17F/minuit_new/all_fits.root", "RECREATE");
  //TFile allFits("plots_17F/minuit/all_fits.root", "RECREATE");
  //TFile allFits("plots_17F/minuit/all_fits_single_dets_freePos.root", "RECREATE");

  for (size_t i = 0; i < histNames.size(); ++i) {
    TH1 *h = dynamic_cast<TH1*>(f->Get(histNames[i].c_str()));
    if (!h) { std::cerr << "Histogram " << histNames[i] << " not found\n"; continue; }

    std::cout << "Fitting " << histNames[i] << " ...\n";
    // automatic initial amplitudes estimated inside fit function, so just pass peakPositions
    fitHistogramLikelihoodTMinuit(h, peakSets[i], &allFits, "plots_17F/minuit_new/fit_results_all.txt", true);
    //fitHistogramLikelihoodTMinuit(h, peakSets[i], &allFits, "plots_17F/minuit/fit_results_all.txt", true);
    //fitHistogramLikelihoodTMinuitFreePos(h, peakSets[i], &allFits, "plots_17F/minuit/fit_results_all_single_freePos.txt", true);
  }

  f->Close();
}
