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

void poster_angular_dist_17F_half_dets_relative_sf() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex0937_1041_1121_half_dets = {
    { 8.40, 15.00, 6.60, 11.70, 1.3384, (13.2355 + 7.30814 + 5.04478)}, // 0
    {15.21, 18.71, 3.50, 16.96, 1.0216, (3.94081 + 11.0333 + 5.58712)}, // 1
    {19.99, 22.61, 2.62, 21.30, 0.9516, (1.69996 + 0 + 14.5829)}, // 2
    {22.74, 25.05, 2.31, 23.90, 0.9366, (1.42418 + 8.92131 + 8.23189)}, // 3
    {25.98, 28.01, 2.03, 27.00, 0.9218, (3.03473 + 0.173152 + 12.7338)}, // 4
    {28.11, 30.00, 1.89, 29.06, 0.9166, (3.5703 + 0. + 9.39902)}, // 5
    {30.77, 32.50, 1.73, 31.64, 0.9099, (3.27858 + 0. + 9.2654)}, // 6
    {32.59, 34.24, 1.65, 33.42, 0.9081, (1.31148 + 0. + 8.8856)}, // 7
    {34.95, 36.50, 1.55, 35.73, 0.9032, (0.449797 + 0.000258 + 2.93755) * 4./3.}, // 8
    {36.58, 38.07, 1.49, 37.32, 0.9044, (0.568829 + 0.01185 + 3.00139) * 4./3.}, // 9
    {38.67, 40.09, 1.42, 39.38, 0.9023, (3.15072) * 4./3.}, //10
    {40.17, 41.54, 1.38, 40.85, 0.8996, (0.677276) * 4./3.} //11
  }; // mean is 1.029

  std::vector<std::vector<double>> ex3061_half_dets = {
    {8.6  , 14.35, 5.75, 11.475, 1.1439,  5.76304}, // 2
    {14.60, 18.40, 3.80, 16.50,  1.0804,  5.7179}, // 3
    {19.76, 22.52, 2.77, 21.14,  0.9979,  4.93428}, // 4
    {22.66, 25.08, 2.42, 23.87,  0.9795,  4.64378}, // 5
    {26.04, 28.15, 2.11, 27.10,  0.9629,  3.41079}, // 6
    {28.26, 30.22, 1.96, 29.24,  0.9577,  3.23499}, // 7
    {31.05, 32.85, 1.80, 31.95,  0.9512,  1.07583 * 4./3.}, // 8
    {32.94, 34.65, 1.70, 33.80,  0.9476,  1.91486 * 4./3.}, // 9
    {35.33, 36.94, 1.61, 36.14,  0.9473,  0.76595 * 4./3.}, //10
    {37.02, 38.56, 1.54, 37.79,  0.9428,  1.63772 * 4./3.} //11
  };

  std::vector<std::vector<double>> ex3724_3839_half_dets = { // 3.79
    { 8.20, 14.74, 6.54, 11.47, 1.3005,  (0.487773 + 20.1041)}, // 3
    {16.54, 19.92, 3.38, 18.23, 1.0572,  (1.75304 + 3.55357)}, // 4
    {20.08, 22.84, 2.77, 21.46, 1.0117,  (3.2814 + 1.80978)}, // 5
    {23.92, 26.25, 2.33, 25.08, 0.9874,  (1.4294 + 3.70091)}, // 6
    {26.36, 28.49, 2.13, 27.43, 0.9791,  (1.14742 + 5.24187)}, // 7
    {29.38, 31.31, 1.92, 30.35, 0.9716,  (0. + 3.77986) * 4./3.}, // 8
    {31.41, 33.22, 1.81, 32.31, 0.9682,  (0.207111 + 2.60804) * 4./3.}, // 9
    {33.94, 35.63, 1.69, 34.79, 0.9627,  (0.8015 + 2.0848) * 4./3.}, //10
    {35.72, 37.33, 1.62, 36.53, 0.9619,  (0.753033 + 1.19958) * 4./3.} //11
  };

  std::vector<std::vector<double>> ex4115_half_dets = {
    {14.71, 18.56, 3.85, 16.63, 1.1012,  10.7236}, // 4
    {18.73, 21.72, 2.99, 20.22, 1.0329,  10.8906}, // 5
    {22.86, 25.31, 2.45, 24.09, 1.0004,  9.79403}, // 6
    {25.43, 27.65, 2.22, 26.54, 0.9901,  6.87214}, // 7
    {28.58, 30.57, 1.99, 29.57, 0.9808,  3.30765 * 4./3.}, // 8
    {30.67, 32.53, 1.86, 31.60, 0.9763,  3.78165 * 4./3.}, // 9
    {33.27, 35.01, 1.74, 34.14, 0.9739,  3.29471 * 4./3.}, //10
    {35.10, 36.75, 1.65, 35.92, 0.9695,  3.32874 * 4./3.} //11
  };

  std::vector<std::vector<double>> ex4360_half_dets = {
    {13.01, 17.42, 4.41, 15.22, 1.1582,  1.23553}, // 4
    {17.61, 20.81, 3.20, 19.21, 1.0534,  1.0575}, // 5
    {22.02, 24.57, 2.55, 23.30, 1.0097,  0.2632}, // 6
    {24.70, 26.99, 2.29, 25.85, 1.0004,  1.25035}, // 7
    {27.95, 29.99, 2.04, 28.97, 0.9885,  0.710929 * 4./3.}, // 8
    {30.09, 32.00, 1.91, 31.05, 0.9841,  0.807541 * 4./3.} // 9
  };

  /*std::vector<std::vector<double>> ex4652_4753_half_dets = { // theta values for 4.652
    {9.60 , 15.90, 6.30, 12.75, 1.3906, (21.0052 + 0.436834)}, // 4
    {16.12, 19.64, 3.52, 17.88, 1.0821, 20.5967}, // 5
    {20.94, 23.64, 2.70, 22.29, 1.0246, 15.1561}, // 6
    {23.77, 26.16, 2.39, 24.97, 1.0100, 13.0476}, // 7
    {27.16, 29.27, 2.11, 28.21, 0.9969, 6.08755 * 4./3.}, // 8
    {29.37, 31.34, 1.96, 30.35, 0.9921, 4.54548 * 4./3.}, // 9
    {32.12, 33.93, 1.81, 33.02, 0.9873, (3.86017) * 4./3.}, //10 - I don't trust this too much
    {34.02, 35.74, 1.72, 34.88, 0.9831, (3.65818 + 0.517828) * 4./3.} //11 - I don't trust this too much
    };*/

  std::vector<std::vector<double>> ex4652_4753_half_dets = { // theta values for 4.702 - avg value of energies
    { 8.20, 15.61, 7.41, 11.91, 1.5286, (21.0052 + 0.436834)}, // 4
    {15.84, 19.43, 3.59, 17.63, 1.0890, 20.5967}, // 5
    {20.74, 23.47, 2.73, 22.10, 1.0284, 15.1561}, // 6
    {23.61, 26.02, 2.42, 24.81, 1.0136, 13.0476}, // 7
    {27.02, 29.14, 2.12, 28.08, 0.9982, 6.08755 * 4./3.}, // 8
    {29.25, 31.22, 1.97, 30.23, 0.9936, 4.54548 * 4./3.}, // 9
    {32.01, 33.82, 1.82, 32.91, 0.9873, (3.86017) * 4./3.}, //10 - I don't trust this too much
    {33.92, 35.64, 1.73, 34.78, 0.9846, (3.65818 + 0.517828) * 4./3.} //11 - I don't trust this too much
  };

  std::vector<std::vector<double>> ex4964_half_dets = {
    {8.2  , 13.94, 5.74, 11.07, 1.1021, 8.57289}, // 4
    {14.20, 18.26, 4.06, 16.23, 1.1354, 7.30343}, // 5
    {19.67, 22.57, 2.90, 21.12, 1.0441, 4.63121}, // 6
    {22.71, 25.23, 2.52, 23.97, 1.0231, 4.8697}, // 7
    {26.27, 28.46, 2.19, 27.36, 1.0074, 3.66391 * 4./3.}, // 8
    {28.57, 30.60, 2.03, 29.58, 1.0016, 2.2864 * 4./3.}, // 9
    {31.40, 33.26, 1.86, 32.33, 0.9963, 3.12153 * 4./3.}, //10 - I don't trust this too much
    {33.36, 35.12, 1.77, 34.24, 0.9939, 1.73812 * 4./3.} //11 - I don't trust this too much
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;

  data = {ex4964_half_dets, ex3061_half_dets, ex4652_4753_half_dets, ex4360_half_dets, ex4115_half_dets, ex3724_3839_half_dets};
  Ex_values = {4.964, 3.062, 4.652, 4.360, 4.115, 3.790};
  labels = {
    {"\\ell = 0", "\\ell = 2"},
    {"\\ell = 0", "\\ell = 2"},
    {"\\ell = 2"},
    {"\\ell = 2"},
    {"\\ell = 0", "\\ell = 2"},
    {"\\ell = 0", "\\ell = 2"},
  };

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("DWBA_17F_%s.root", pot.Data());
  TString outputDir = Form("plots_17F/minuit_new/ang_dist_fixedRatio/%s", pot.Data());
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
    {ex4964_half_dets, {2, 3}},//, 11, 12}},
    {ex3061_half_dets, {5, 6}},
    {ex4652_4753_half_dets, {29, 30}},
    {ex4360_half_dets, {21}},
    {ex4115_half_dets, {24, 23}},
    {ex3724_3839_half_dets, {28, 26}}
  };

  std::vector<std::pair<std::vector<int>, std::vector<double>>> fitPairs;
  fitPairs =  {
    {{2, 3}, {0.35, 0.66}}, // 4.694, 1s1/2 and 0d5/2 SF ratio like in 18O
    {{5, 6}, {0.21, 0.83}}, // 3.061, 1s1/2 and 0d5/2 SF ratio like in 18O
    {{29, 30}, {1.57, 0.28}}, // 4.652 & 4.753, SF like in 18O for 3.552 & 3.630
    {{21}, {1.0}}, // 4.360
    {{24, 23}, {1.0, 1.0}}, // 4.115
    {{28, 26}, {1.0, 1.0}} // 3.724 & 3.839
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

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.85);
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
