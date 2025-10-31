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
  // par[1] = component index (as integer, but stored in double)
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

void poster_angular_dist_17O_half_dets_relative_sf() {

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

  std::vector<std::vector<double>> ex3552_3630_half_dets = {
    {10.92, 16.29, 5.37, 13.60, 1.2636, 63.0214 + 7.39285, 3.9}, // 4 , 3.9 + 346791 but looking at the fit I think it makes sense to take only the bigger one
    {16.49, 19.92, 3.42, 18.20, 1.0694, 66.6946, 2.21898}, // 5
    {21.17, 23.83, 2.66, 22.50, 1.0161, 48.8069, 1.86560}, // 6
    {23.96, 26.31, 2.35, 25.14, 0.9992, 37.5365, 1.63665}, // 7
    {27.29, 29.38, 2.08, 28.33, 0.9889, 20.0825 * 4./3., 1.0096}, // 8
    {29.48, 31.42, 1.94, 30.45, 0.9840, 15.6256 * 4./3., 1.0757}, // 9
    {32.20, 33.99, 1.80, 33.09, 0.9803, 11.3624 * 4./3., 0.9295}, //10
    {34.08, 35.79, 1.70, 34.94, 0.9760, 10.1256 * 4./3., 0.8921}, //11
  };

  std::vector<std::vector<double>> ex5255_5340_5375_half_dets = {
    {11.04, 16.48, 5.44, 13.76, 1.2938, 10.36 + 0 + 21.1155, sqrt( pow(1.47212, 2) + pow(1.7236, 2) )}, // 6
    {16.69, 20.20, 3.51, 18.45, 1.1096, 5.55457 + 0 + 3.17525, sqrt( pow(0.9439842, 2) + pow(1.2226, 2) )}, // 7
    {21.54, 24.26, 2.72, 22.90, 1.0595, (8.63332 + 0.725115 + 0.01896) * 4./3., 1.00871}, // 8
    {24.39, 26.81, 2.42, 25.60, 1.0449, (9.73567 + 4.43685 + 0) * 4./3., sqrt( pow(1.46941, 2) + pow(1.32538, 2) )}, // 9
    {27.75, 29.90, 2.14, 28.82, 1.0339, (12.7922 + 0 + 2.47708) * 4./3., sqrt( pow(1.45785, 2) + pow(1.18, 2) )}, //10
    {30.01, 32.00, 2.00, 31.00, 1.0281, (13.0306 + 2.34372 + 0) * 4./3., sqrt( pow(2.3225, 2) + pow(2.16732, 2) )}, //11
  }; // in Cleopatra/FindThetaCM energy 5.315 was used as the mean value of 5.255 and 5.375 states in the triplet

  std::vector<std::vector<double>> ex6200_half_dets = {
    { 8.00, 15.58, 7.58, 11.75, 1.5436, 2.200021, 0.455976}, // 7
    {17.42, 20.85, 3.43, 19.14, 1.1240, 1.77982 * 4./3., 0.324562}, // 8
    {21.01, 23.86, 2.85, 22.44, 1.0864, 1.23201 * 4./3., 0.3665}, // 9
    {24.94, 27.36, 2.42, 26.15, 1.0655, 1.04040 * 4./3., 0.354472}, //10
    {27.48, 29.69, 2.21, 28.59, 1.0593, 0.5536  * 4./3., 0.260893}, //11
  };

  std::vector<std::vector<double>> ex6930_half_dets = {
    {12.15, 17.19, 5.04, 14.67, 1.2773, 3.80869 * 4./3., 0.46119}, // 8
    {17.40, 20.89, 3.49, 19.14, 1.1460, 4.12727 * 4./3., 0.587297}, // 9
    {22.15, 24.91, 2.76, 23.53, 1.1025, 2.44996 * 4./3., 0.489051}, //10
    {25.05, 27.51, 2.46, 26.28, 1.0869, 2.41137 * 4./3., 0.471051}, //11
  };

  std::vector<std::vector<double>> ex7100_half_dets = {
    { 9.74, 16.15, 6.40, 12.95, 1.4348, 3.80869 * 4./3., 0.46119}, // 8
    {16.37, 20.10, 3.73, 18.23, 1.1665, 4.12727 * 4./3., 0.58729}, // 9
    {21.42, 24.28, 2.86, 22.85, 1.1122, 2.44996 * 4./3., 0.48905}, //10
    {24.42, 26.95, 2.53, 25.69, 1.0981, 2.41137 * 4./3., 0.47105}, //11
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;

  data = {ex1982_half_dets, ex3920_half_dets, ex3552_3630_half_dets, ex5255_5340_5375_half_dets, ex6200_half_dets, ex6930_half_dets, ex7100_half_dets};
  Ex_values = {1.982, 3.920, 3.552, 5.255, 6.200, 6.930, 6.930, 6.930, 6.930, 6.930, 7.100};
  labels = {
    {"\\ell = 0, 0d5/2", "\\ell = 2, 1s1/2"}, // 1.982
    {"\\ell = 0, 1s1/2", "\\ell = 2, 0d5/2", "\\ell = 3", "\\ell = 4"}, // 3.920
    {"\\ell = 2, 0d5/2\\;4^{+}, 3.552", "\\ell = 2, 0d5/2\\;0^{+}, 3.630",}, // 3.552 & 3.630
    {"\\ell = 2, 0d5/2\\;0^{+}, 5.340", "\\ell = 0, 1s1/2\\;2^{+}, 5.255", "\\ell = 0, 1s1/2\\;3^{+}, 5.375"}, // 5.255 & 5.340 & 5.375
    {"\\ell = 1, 1p3/2\\;1^{-}"}, // 6200
    {"\\ell = 1, 1p3/2\\;1^{-}"}, // 6930, 30
    {"\\ell = 3, 0f7/2\\;1^{-}"}, // 6930, 31
    {"\\ell = 2, 0d5/2\\;0^{+}"}, // 6930, 34
    {"\\ell = 0, 1s1/2\\;2^{+}"}, // 6930, 35
    {"\\ell = 2, 0d5/2\\;0^{+}", "\\ell = 0, 1s1/2\\;2^{+}"}, // 6930, 34, 35
    {"\\ell = 2, 0d5/2\\;4^{+}"}, // 7100, 35
  };

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("DWBA_17O_%s.root", pot.Data());
  TString outputDir = Form("plots_17O/minuit_extra5/ang_dist_fixedRatio_newunc/%s", pot.Data());
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
    {ex1982_half_dets, {0, 2}},//, 11, 12}},
    {ex3920_half_dets, {15, 13}},
    {ex3552_3630_half_dets, {36, 37}},
    {ex5255_5340_5375_half_dets, {26, 27, 28}},
    {ex6200_half_dets, {29}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex7100_half_dets, {38}}
  };

  std::vector<std::pair<std::vector<int>, std::vector<double>>> fitPairs;
  fitPairs =  {
    {{0, 2}, {0.21, 0.83}}, // 1.982, 1s1/2 and 0d5/2 relative SF
    {{15, 13}, {0.35, 0.66}}, // 3.920, 1s1/2 and 0d5/2 relative SF
    {{36, 37}, {1.57, 0.28}}, // 3.552 & 3.630, 4+ and 0+ relative SF
    {{26, 27, 28}, {0.16, 0.35, 1.01}}, // 5.255, 5.340, 5.375, 0+, 2+, 3+ relative SF, so 5.340, 5.255 and 5.375, respectively
    {{29}, {1.0}}, // 6.200
    {{30}, {1.0}}, // 6.930
    {{31}, {1.0}}, // 6.930
    {{34}, {1.0}}, // 6.930
    {{35}, {1.0}}, // 6.930
    {{34, 35}, {1.0, 1.0}}, // 6.930
    {{38}, {1.0}} // 7.100
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

    TCanvas* canvas = new TCanvas(Form("fit_canvas_%lu", mappingIndex + 1), 
				  Form("Ex = %.3f MeV", Ex_values[mappingIndex]), 1400, 1400);

    experimentGraph->Draw("APE1");

    experimentGraph->GetHistogram()->GetXaxis()->SetRangeUser(0, 60);
    experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.01, 100);
    if (mappingIndex == 2)
      experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.4, 120);

    TAxis *axis = experimentGraph->GetXaxis();
    axis->SetLimits(0.,60.);
    
    experimentGraph->SetTitle(Form("Ex = %.3f MeV (fixed SF);#theta_{CM} (deg);d#sigma/d#Omega (a. u.)", Ex_values[mappingIndex]));

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
