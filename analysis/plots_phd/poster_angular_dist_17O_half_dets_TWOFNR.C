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

// extracrting SF with TWOFNR instead of Ptolemy for the 1.982 state

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};

std::vector<TString> potentials = {
  "AK", "AV", "AM", "AG", "AP",
  "HK", "HV", "HM", "HG", "HP",
  "BK", "BV", "BM", "BG", "BP",
  "DK", "DV", "DM", "DG", "DP",
  "QK", "QV", "QM", "QG", "QP",
  "ZK", "ZV", "ZM", "ZG", "ZP",
  "LK", "LV", "LM", "LG", "LP"
};

bool checking_gs = false;

double combinedDWBA(double *x, double *par) {
  double result = 0;
  for (size_t i = 0; i < nr_of_functions; ++i) {
    double dwba_val = graphsDWBA[functions_ids[i]]->Eval(x[0]);
    result += par[i] * dwba_val;
  }
  return result;
}

void poster_angular_dist_17O_half_dets_TWOFNR() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex1982_half_dets = {
    {8.00 , 14.62, 6.62, 11.31, 1.2983, 28.4242, 1.42521}, // 2
    {14.86, 18.56, 3.70, 16.71, 1.0635, 17.2772, 1.11639}, // 3
    {19.89, 22.62, 2.73, 21.26, 0.9886, 16.2642, 1.07688}, // 4
    {22.76, 25.15, 2.39, 23.95, 0.9711, 12.7089, 0.99956}, // 5
    {26.10, 28.20, 2.10, 27.15, 0.9562, 10.2986, 0.86965}, // 6
    {28.30, 30.25, 1.94, 29.27, 0.9510, 10.5578, 0.88115}, // 7
    {31.07, 32.86, 1.78, 31.97, 0.9448, 6.15872 * 4./3., 0.559427}, // 8
    {32.95, 34.64, 1.69, 33.80, 0.9410, 4.24311 * 4./3., 0.590164}, // 9
    {35.32, 36.92, 1.60, 36.12, 0.9407, 3.42685 * 4./3., 0.524537}, //10
    {37.00, 38.53, 1.53, 37.77, 0.9362, 2.27111 * 4./3., 0.419665}, //11
  };

  std::vector<std::vector<double>> ex3920_half_dets = {
    { 8.00, 14.08, 6.08, 11.04, 1.1643, 26.7301, 1.43482}, // 4
    {14.33, 18.32, 4.00, 16.33, 1.1236, 22.9592, 1.34987}, // 5
    {19.72, 22.59, 2.87, 21.15, 1.0362, 14.5072, 1.03633}, // 6
    {22.73, 25.23, 2.50, 23.98, 1.0159, 15.1933, 1.05784}, // 7
    {26.26, 28.44, 2.18, 27.35, 1.0007, 8.45285 * 4./3., 0.67275}, // 8
    {28.55, 30.56, 2.02, 29.56, 0.9950, 9.04236 * 4./3., 0.8328}, // 9
    {31.36, 33.22, 1.85, 32.29, 0.9899, 7.50488 * 4./3., 0.7676}, //10
    {33.31, 35.07, 1.76, 34.19, 0.9877, 6.11097 * 4./3., 0.711327}, //11
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;
  std::vector<TString> titles;

  data = {ex1982_half_dets};
  Ex_values = {1.982};

  labels = {
    {"0 #font[42]{d}_{5/2}, 2^{+}", "1 #font[42]{s}_{1/2} 2^{+}"}
  };
  
  titles = {"1.982 MeV"};

  std::vector<TGraph*> graphs;

  TFile *file = TFile::Open("twofnr/output_17Odp18O.root");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Cannot open output_17Odp18O.root" << std::endl;
    return;
  }

  TGraph* g_d52 = (TGraph*)file->Get("0d52E1982");
  TGraph* g_s12 = (TGraph*)file->Get("1s12E1982");

  if (!g_d52 || !g_s12) {
    std::cerr << "ERROR: Required DWBA graphs not found!" << std::endl;
    return;
  }

  graphsDWBA.push_back(g_d52); // index 0
  graphsDWBA.push_back(g_s12); // index 1


  std::cout << "Loaded DWBA graphs:" << std::endl;
  std::cout << "  [0] " << g_d52->GetName() << std::endl;
  std::cout << "  [1] " << g_s12->GetName() << std::endl;


  const char* graphTitles[] = {
    "E_{x} = 1.982 MeV, 0 #font[42]{d}_{5/2}, 2^{+}"
  };

  file->Close();

  std::vector<int> color = {629, 596, 418, 801, 905, 8, 9, 1};
  int ncolor = 0;

  std::vector<double> theta_means;
  std::vector<double> int_corr_sins;
  std::vector<double> int_unc;
  std::vector<double> x_unc;

  for (const auto& det : ex1982_half_dets) {
    Tmin = det[0];
    Tmax = det[1];
    Dt = det[2];
    ThetaMean = det[3];
    sin_x_dx = det[4];
    integral_corr = det[5];
    integral_unc = det[6];
      
    theta_means.push_back(det[3]);
    int_corr_sins.push_back(det[5] / det[4] );
    //int_unc.push_back(TMath::Sqrt(integral_corr / sin_x_dx));
    int_unc.push_back(integral_unc / sin_x_dx);
    x_unc.push_back(0.0);
  }

  TGraphErrors* experimentGraph = new TGraphErrors(theta_means.size(), &theta_means[0], &int_corr_sins[0], &x_unc[0], &int_unc[0]);
  experimentGraph->SetTitle("");
  experimentGraph->SetMarkerStyle(20);
  experimentGraph->SetMarkerSize(2);
  experimentGraph->SetMarkerColor(kBlack);
  experimentGraph->SetLineColor(kBlack);
  experimentGraph->SetLineWidth(2);
  //experimentGraph->SetErrorSize(1.5);

  TCanvas* canvas = new TCanvas(Form("fit_canvas_%lu", 1), 
				Form("Ex = %.3f MeV", Ex_values[0]), 1400*2, 1400*2);

  experimentGraph->Draw("APE1 SAME");

  experimentGraph->GetHistogram()->GetXaxis()->SetRangeUser(0, 60);
  experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.01, 180);
  //if (mappingIndex == 2)
  //experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.4, 120);

  TAxis *axis = experimentGraph->GetXaxis();
  axis->SetLimits(0.,60.);
    
  experimentGraph->SetTitle(Form("%s;#theta_{CM} (deg);d#sigma/d#Omega (a. u.)", titles[0].Data()));

  gPad->SetLeftMargin(0.15);   // więcej miejsca na tytuł osi Y
  gPad->SetBottomMargin(0.15); // więcej miejsca na tytuł osi X
  gPad->SetRightMargin(0.05);  // wąski prawy margines
  gPad->SetTopMargin(0.1);

  experimentGraph->GetXaxis()->SetTitleSize(0.044);
  experimentGraph->GetYaxis()->SetTitleSize(0.044);

  experimentGraph->GetXaxis()->SetTitleOffset(1.2);
  experimentGraph->GetYaxis()->SetTitleOffset(1.5);

  experimentGraph->GetXaxis()->SetLabelSize(0.04);
  experimentGraph->GetYaxis()->SetLabelSize(0.04);

  experimentGraph->GetXaxis()->SetLabelOffset(0.01);
  experimentGraph->GetYaxis()->SetLabelOffset(0.01);

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(22);
  //latex.DrawLatexNDC(0.5, 0.92, Form("Ex = %.3f MeV", Ex_values[mappingIndex]));

  TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.89);
  legend->SetTextSize(0.04);
  legend->SetTextFont(42);
  legend->SetMargin(0.1);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  nr_of_functions = 2;
  functions_ids[0] = 0; // 0d5/2
  functions_ids[1] = 1; // 1s1/2

  TF1* fitFunction = new TF1(
			 "fitFunc_1982",
			 combinedDWBA,
			 0.0,
			 60.0,
			 2
			 );

  fitFunction->SetParameter(0, 1.0);
  fitFunction->SetParameter(1, 1.0);
  fitFunction->SetParLimits(0, 0.0, 100.0);
  fitFunction->SetParLimits(1, 0.0, 100.0);

  TFitResultPtr fitResult = experimentGraph->Fit(fitFunction, "RLS0");

  double chi2 = 1000000.;

  if (fitResult->IsValid()) {
    chi2 = fitResult->Chi2();
  }
  else {
    std::cout << "Invalid fit " << std::endl;
  }
      
  fitFunction->SetLineColor(color[ncolor]);
  fitFunction->SetLineWidth(4);
  if (nr_of_functions > 1)
    fitFunction->Draw("L SAME");

  if (nr_of_functions > 1)
    legend->AddEntry(fitFunction, "Total fit", "l");

  for (size_t comp = 0; comp < nr_of_functions; ++comp) {
    TF1* compFunc = new TF1(Form("fitComponent_%lu", comp),
			    combinedDWBA, 0, 60, 2);

    for (size_t i = 0; i < 2; ++i) {
      if (i == comp) compFunc->SetParameter(i, fitFunction->GetParameter(i));
      else compFunc->SetParameter(i, 0.0);
    }

    compFunc->SetLineColor(color[(ncolor + comp + 1) % color.size()]);
    compFunc->SetLineStyle(2);
    compFunc->SetLineWidth(4);
    compFunc->Draw("L SAME");

    legend->AddEntry(compFunc, Form("%s", labels[0][comp].Data()), "l");
  }

      

  std::cout << "======================================" << std::endl;
  std::cout << " E_x = 1.982 MeV " << std::endl;
  std::cout << " S(0d5/2) = " << fitFunction->GetParameter(0)
	    << " ± " << fitFunction->GetParError(0) << std::endl;
  std::cout << " S(1s1/2) = " << fitFunction->GetParameter(1)
	    << " ± " << fitFunction->GetParError(1) << std::endl;
  std::cout << " Chi2     = " << fitResult->Chi2() << std::endl;
  std::cout << " NDF      = " << fitResult->Ndf() << std::endl;
  std::cout << "======================================" << std::endl;

      

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.05);
  latex2.SetTextFont(42);          // 42 = Helvetica, 62 = bold Helvetica
  latex2.SetTextColor(kRed);
  latex2.SetTextColorAlpha(kGray+2, 0.4);
  latex2.DrawLatex(0.25, 0.22, "PRELIMINARY");

  //experimentGraph->Draw("APE1 SAME");

  canvas->SetLogy();
  //canvas->Update();

  legend->Draw();

  experimentGraph->Draw("PE SAME");

  gSystem->mkdir("plots_17O/Ex1982_TWOFNR", kTRUE);
  canvas->SaveAs("plots_17O/Ex1982_TWOFNR/fit_Ex1982.png");
  canvas->SaveAs("plots_17O/Ex1982_TWOFNR/fit_Ex1982.pdf");
}
