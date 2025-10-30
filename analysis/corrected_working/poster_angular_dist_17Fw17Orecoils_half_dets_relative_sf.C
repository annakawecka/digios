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
#include <TMath.h>

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};
std::vector<double> functions_relative_strengths;

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

double combinedDWBA_fixedRatio(double *x, double *par) {
  double result = 0;
  double A = par[0];
  for (size_t i = 0; i < nr_of_functions; ++i) {
    int gid = functions_ids[i];
    if ((size_t)gid >= graphsDWBA.size()) continue;
    double dwba_val = graphsDWBA[functions_ids[i]]->Eval(x[0]);
    result += A * functions_relative_strengths[i] * dwba_val;
  }
  return result;
}

double singleDWBA_component(double *x, double *par) {
  // par[0] = normalization A
  // par[1] = component index
  int idx = static_cast<int>(par[1]);
  if (idx < 0 || idx >= (int)nr_of_functions) return 0.0;
  int gid = functions_ids[idx];
  if ((size_t)gid >= graphsDWBA.size()) return 0.0;

  double dwba_val = graphsDWBA[gid]->Eval(x[0]);
  return par[0] * functions_relative_strengths[idx] * dwba_val;
}

std::vector<std::vector<int>> generateCombinations(const std::vector<int>& indices) {
  std::vector<std::vector<int>> combinations;
  size_t n = indices.size();
  for (size_t k = 1; k <= n; ++k) {
    std::vector<bool> select(n);
    std::fill(select.begin(), select.begin() + k, true);
    do {
      std::vector<int> subset;
      for (size_t i = 0; i < n; ++i) {
	if (select[i]) {
	  subset.push_back(indices[i]);
	}
      }
      combinations.push_back(subset);
    } while (std::prev_permutation(select.begin(), select.end()));
  }
  return combinations;
}

void poster_angular_dist_17Fw17Orecoils_half_dets_relative_sf() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex6136_6163_half_dets = {
    {12.90, 17.50, 4.61, 15.20, 1.2078, 4.01808 + 0.}, // 6
    {17.70, 21.02, 3.32, 19.36, 1.0991, 1.87111 + 0.}, // 7
    {22.30, 24.94, 2.64, 23.62, 1.0577, (3.28537 + 0.) * 4./3.}, // 8
    {25.07, 27.43, 2.36, 26.25, 1.0417, (4.51246 + 0.) * 4./3.}, // 9
    {28.35, 30.45, 2.10, 29.40, 1.0331, (5.23575 + 0.) * 4./3.}, //10
    {30.56, 32.52, 1.97, 31.54, 1.0286, (4.94403 + 0.) * 4./3.}, //11
  }; // thetaCM calculated for (6.136 + 6.163) / 2.0 = 6.1495
  
  std::vector<std::vector<double>> ex6633_6643_half_dets = {
    { 8.00, 14.53, 6.53, 11.26, 1.2751, 0.00013 + 14.4628}, // 6
    {14.80, 18.82, 4.02, 16.81, 1.1635, 0. + 9.81058}, // 7
    {20.29, 23.21, 2.93, 21.75, 1.0841, (0. + 3.84624) * 4./3.}, // 8
    {23.35, 25.91, 2.56, 24.63, 1.0656, (0. + 4.66187) * 4./3.}, // 9
    {26.90, 29.13, 2.24, 28.01, 1.0506, (0.00018577 + 3.37614) * 4./3.}, //10
    {29.24, 31.32, 2.07, 30.28, 1.0447, (0.0 + 3.12785)  * 4./3.}, //11
  }; // thetaCM calculated for (6.633 + 6.643) / 2. = 6.638

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;

  data = {ex6136_6163_half_dets, ex6633_6643_half_dets};
  Ex_values = {6.136, 6.633};
  labels = {
    {"\\ell = 2, 0d5/2\\;0^{+}, 6.136", "\\ell = 0, 1s1/2\\;3^{+}, 6.163"}, // 6136_6163 will probably be IAS of thr 5.34 (d) and 5.38 (s) states in O; there might be also addition of 6.108 state in F, I'm not sure what l is that, but also 6.108 is T=0 state so it's weaker
    {"\\ell = 3, 1f7/2\\;2^{-}, 6.163"} // 6.633 state - parity and T unknown, 6.643 - 2- state but T=1 so stronger? I don't know which state could it possibly correspond to in oxygen - maybe the 6.35 state? It's also 2-, but very weak
  };
  
  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("DWBA_17F_%s_he.root", pot.Data());
  TString outputDir = Form("plots_17Fw17Orecoils/minuit_extra5/ang_dist_fixedRatio/%s", pot.Data());
  gSystem->mkdir(outputDir, kTRUE);
  
  if (checking_gs) {
    //TFile *file = TFile::Open("DWBA_17O_pot_5255.root");
    file = TFile::Open("DWBA_17O_pot_gs.root");
  }
  else
    file = TFile::Open(filename);
  
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening ROOT file!" << std::endl;
    return;
  }

  TObjArray* objArray = (TObjArray*)file->Get("qList");
  if (!objArray) {
    std::cerr << "Error: Object array not found in file " << std::endl;
    file->Close();
    return;
  }

  for (int i = 0; i < objArray->GetEntries(); ++i) {
    TObject* obj = objArray->At(i);
    if (obj && obj->InheritsFrom(TGraph::Class())) {
      TGraph* graph = (TGraph*)obj;
      graphsDWBA.push_back(graph);

      std::cout << "Found TGraph: " << graph->GetName() << std::endl;
    }
  }

  std::vector<std::pair<std::vector<std::vector<double>>, std::vector<int>>> fitMappings;
  fitMappings = {
    {ex6136_6163_half_dets, {0, 1}},
    {ex6633_6643_half_dets, {2}}
  };

  std::vector<std::pair<std::vector<int>, std::vector<double>>> fitPairs;
  fitPairs =  {
    {{0, 1}, {0.16, 1.01}},
    {{2}, {1.0}}, 
  };

  std::vector<int> color = {629, 596, 418, 801, 905, 8, 9, 1, 49, 42, 40};
  int ncolor = 0;

  std::ofstream resultFile;

  TString outputFilename = Form("%s/DWBA_fit_results_17O.txt", outputDir.Data());
  resultFile.open(outputFilename);

  for (size_t mappingIndex = 0; mappingIndex < fitMappings.size(); ++mappingIndex) {
    const auto& dataSet = fitMappings[mappingIndex].first;
    const auto& fitIndices = fitMappings[mappingIndex].second;

    nr_of_functions = fitPairs[mappingIndex].first.size();
    functions_relative_strengths = fitPairs[mappingIndex].second;
    for (int i = 0; i < nr_of_functions; ++i)
      functions_ids[i] = fitPairs[mappingIndex].first[i];

    double bestChi2 = 1e9;
    std::vector<int> bestCombination;
    TF1* bestFitFunction = nullptr;
    //ncolor = 0;

    std::vector<double> theta_means;
    std::vector<double> int_corr_sins;
    std::vector<double> int_unc;
    std::vector<double> x_unc;

    for (const auto& det : dataSet) {
      Tmin = det[0];
      Tmax = det[1];
      Dt = det[2];
      ThetaMean = det[3];
      sin_x_dx = det[4];
      integral_corr = det[5];
      integral_unc = det[6];
      
      theta_means.push_back(det[3]);
      int_corr_sins.push_back(det[5] / det[4] );
      int_unc.push_back(TMath::Sqrt(integral_corr / sin_x_dx));
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

    TCanvas* canvas = new TCanvas(Form("fit_canvas_%lu", mappingIndex + 1), 
				  Form("Ex = %.3f MeV", Ex_values[mappingIndex]), 1400, 1400);

    experimentGraph->Draw("APE1");

    experimentGraph->GetHistogram()->GetXaxis()->SetRangeUser(0, 60);
    experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.01, 100);
    if (mappingIndex == 2)
      experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.4, 120);

    TAxis *axis = experimentGraph->GetXaxis();
    axis->SetLimits(0.,60.);
    
    experimentGraph->SetTitle(Form("Ex = %.3f MeV;#theta_{CM} (deg);d#sigma/d#Omega (a. u.)", Ex_values[mappingIndex]));

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
    legend->SetMargin(0.1);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);

    //TF1* fitFunction = new TF1("fitFunction", combinedDWBA, 0, 60, subset.size());
    TF1* fitFunction = new TF1(Form("fitFunction_%lu", mappingIndex), combinedDWBA_fixedRatio, 0, 60, 1);
    fitFunction->SetParameter(0, 1.0);
    fitFunction->SetParLimits(0, 0.0, 100.0);

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
    fitFunction->Draw("L SAME");

    if (nr_of_functions > 1)
      legend->AddEntry(fitFunction, "Total fit", "l");

    for (size_t comp = 0; comp < nr_of_functions; ++comp) {
      TF1* compFunc = new TF1(Form("fitComponent_%lu_%lu", mappingIndex, comp),
			      singleDWBA_component, 0, 60, 2);
      compFunc->SetParameter(0, fitFunction->GetParameter(0)); // same normalization
      compFunc->SetParameter(1, comp); // which DWBA component to draw
      compFunc->SetLineColor(color[(ncolor + comp + 1) % color.size()]);
      compFunc->SetLineStyle(2);
      compFunc->SetLineWidth(3);
      compFunc->Draw("L SAME");

      legend->AddEntry(compFunc, Form("%s", labels[mappingIndex][comp].Data()), "l");
    }

    resultFile << "========================= Fit for Ex = " << Ex_values[mappingIndex] << " MeV =========================\n";
    resultFile << "DWBA function indices: ";
    for (auto id : fitPairs[mappingIndex].first) resultFile << id << " ";
    resultFile << std::endl;
    for (size_t comp = 0; comp < nr_of_functions; ++comp) {
      int gid = functions_ids[comp];
      TString gname = graphsDWBA[gid]->GetName();  // or custom name if you have it
      resultFile << "  [" << comp << "] "
	   << gname.Data()
	   << "   SF_rel = " << functions_relative_strengths[comp]
	   << std::endl;
    }
    
    resultFile << "\nRelative SF ratios: ";
    for (auto r : fitPairs[mappingIndex].second) resultFile << r << " ";
    resultFile << "\nNormalization A = " << fitFunction->GetParameter(0)
               << " ± " << fitFunction->GetParError(0)
               << "\nChi2 = " << chi2 << "\n\n";

    ncolor++;
      

    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.05);
    latex2.SetTextFont(42);          // 42 = Helvetica, 62 = bold Helvetica
    latex2.SetTextColor(kRed);
    latex2.SetTextColorAlpha(kGray+2, 0.4);
    latex2.DrawLatex(0.25, 0.22, "PRELIMINARY");

    canvas->SetLogy();
    //canvas->Update();

    legend->Draw();

    TString outputFilename = Form("%s/fit_Ex_%lu.png", outputDir.Data(), mappingIndex + 1);
    canvas->SaveAs(outputFilename);
  }
  
  resultFile.close();
}
