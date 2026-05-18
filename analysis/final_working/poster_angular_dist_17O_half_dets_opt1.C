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
#include <iostream>

// this is the same as ../corrected_working/poster_angular_dist_17O_half_dets_opt1.C

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};
double min_y[13] = {0.2, 0.2, 2, 0.1, 0.1, 0.1, 1, 1, 1, 0.2, 0.1, 0.5, 0.2};
double max_y[13] = {110, 110, 110, 110, 10, 10, 10, 10, 10, 50, 30, 10, 10};

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

void poster_angular_dist_17O_half_dets_opt1() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex1982_half_dets = {
    {8.00 , 14.74, 6.74, 11.37, 1.3287, 28.4242, 1.42521}, // 2
    {14.74, 18.48, 3.73, 16.61, 1.0674, 17.2772, 1.11639}, // 3
    {19.97, 22.69, 2.71, 21.33, 0.9874, 16.2642, 1.07688}, // 4
    {22.69, 25.09, 2.40, 23.89, 0.9722, 12.7089, 0.99956}, // 5
    {26.16, 28.25, 2.09, 27.20, 0.9562, 10.2986, 0.86965}, // 6
    {28.25, 30.20, 1.95, 29.22, 0.9511, 10.5578, 0.88115}, // 7
    {31.12, 32.90, 1.78, 32.01, 0.9453, 6.15872 * 4./3., 0.559427 * 4./3.}, // 8
    {32.90, 34.60, 1.69, 33.75, 0.9413, 4.24311 * 4./3., 0.590164 * 4./3.}, // 9
    {35.37, 36.96, 1.59, 36.16, 0.9410, 3.42685 * 4./3., 0.524537 * 4./3.}, //10
    {36.96, 38.49, 1.53, 37.73, 0.9365, 2.27111 * 4./3., 0.419665 * 4./3.}, //11
  };

  std::vector<std::vector<double>> ex3920_half_dets = {
    { 8.00, 14.20, 6.20, 11.10, 1.1936, 26.7301, 1.43482}, // 4
    {14.20, 18.24, 4.03, 16.22, 1.1268, 22.9592, 1.34987}, // 5
    {19.80, 22.66, 2.86, 21.23, 1.0352, 14.5072, 1.03633}, // 6
    {22.66, 25.17, 2.51, 23.91, 1.0168, 15.1933, 1.05784}, // 7
    {26.32, 28.49, 2.17, 27.41, 1.0007, 8.45285 * 4./3., 0.67275 * 4./3.}, // 8
    {28.49, 30.51, 2.02, 29.50, 0.9951, 9.04236 * 4./3., 0.8328 * 4./3.}, // 9
    {31.41, 33.26, 1.85, 32.34, 0.9897, 7.50488 * 4./3., 0.7676 * 4./3.}, //10
    {33.26, 35.02, 1.76, 34.14, 0.9879, 6.11097 * 4./3., 0.711327 * 4./3.} //11
  };

  std::vector<std::vector<double>> ex3552_3630_half_dets = {
    {10.56, 16.16, 5.60, 13.36, 1.2937, 63.0214 + 7.39285, 3.9}, // 4 , 3.9 + 346791 but looking at the fit I think it makes sense to take only the bigger one
    {16.16, 19.65, 3.49, 17.91, 1.0734, 66.6946, 2.21898}, // 5
    {21.09, 23.75, 2.67, 22.42, 1.0168, 48.8069, 1.86560}, // 6
    {23.75, 26.13, 2.38, 24.94, 1.0034, 37.5365, 1.63665}, // 7
    {27.23, 29.32, 2.09, 28.28, 0.9900, 20.0825 * 4./3., 1.0096 * 4./3.}, // 8
    {29.32, 31.27, 1.95, 30.30, 0.9855, 15.6256 * 4./3., 1.0757 * 4./3.}, // 9
    {32.15, 33.95, 1.80, 33.05, 0.9809, 11.3624 * 4./3., 0.9295 * 4./3.}, //10
    {33.95, 35.66, 1.71, 34.80, 0.9765, 10.1256 * 4./3., 0.8921 * 4./3.} //11
  };

  std::vector<std::vector<double>> ex5255_5340_5375_half_dets = {
    {11.25, 16.59, 5.34, 13.92, 1.2936, 31.89, 1.49409}, // 6 manual fit with just 1 gaus
    {16.59, 20.12, 3.53, 18.35, 1.1122, 5.55457 + 0 + 3.17525, sqrt( pow(0.9439842, 2) + pow(1.2226, 2) )}, // 7
    {21.62, 24.33, 2.71, 22.97, 1.0596, (8.63332 + 0.725115 + 0.01896) * 4./3., 1.00871 * 4./3.}, // 8
    {24.33, 26.75, 2.42, 25.54, 1.0447, (9.73567 + 4.43685 + 0) * 4./3., (sqrt( pow(1.46941, 2) + pow(1.32538, 2) )) * 4./3.}, // 9
    {27.81, 29.95, 2.14, 28.88, 1.0339, (12.7922 + 0 + 2.47708) * 4./3., (sqrt( pow(1.45785, 2) + pow(1.18, 2) )) * 4./3.}, //10
    {29.95, 31.95, 2.00, 30.95, 1.0282, (13.0306 + 2.34372 + 0) * 4./3., (sqrt( pow(2.3225, 2) + pow(2.16732, 2) )) * 4./3.}, //11
  }; // in Cleopatra/FindThetaCM energy 5.315 was used as the mean value of 5.255 and 5.375 states in the triplet

  std::vector<std::vector<double>> ex6200_half_dets = {
    { 8.00, 15.46, 7.46, 11.73, 1.5166, 2.200021, 0.455976}, // 7
    {17.52, 20.93, 3.41, 19.23, 1.1233, 1.77982 * 4./3., 0.324562 * 4./3.}, // 8
    {20.93, 23.79, 2.86, 22.36, 1.0870, 1.23201 * 4./3., 0.3665 * 4./3.}, // 9
    {25.01, 27.42, 2.41, 26.21, 1.0648, 1.04040 * 4./3., 0.354472 * 4./3.}, //10
    {27.42, 29.64, 2.22, 28.53, 1.0593, 0.5536  * 4./3., 0.260893 * 4./3.}, //11
  };

  std::vector<std::vector<double>> ex6930_half_dets = {
    {12.15, 17.19, 5.04, 14.67, 1.2773, 3.80869 * 4./3., 0.46119}, // 8
    {17.40, 20.89, 3.49, 19.14, 1.1460, 4.12727 * 4./3., 0.587297}, // 9
    {22.15, 24.91, 2.76, 23.53, 1.1025, 2.44996 * 4./3., 0.489051}, //10
    {25.05, 27.51, 2.46, 26.28, 1.0869, 2.41137 * 4./3., 0.471051}, //11
  };

  std::vector<std::vector<double>> ex7100_half_dets = {
    {10.12, 16.26, 6.14, 13.19, 1.4011, 3.80869 * 4./3., 0.46119 * 4./3.}, // 8
    {16.26, 20.01, 3.75, 18.14, 1.1686, 4.12727 * 4./3., 0.58729 * 4./3.}, // 9
    {21.50, 24.35, 2.85, 22.92, 1.1116, 2.44996 * 4./3., 0.48905 * 4./3.}, //10
    {24.35, 26.89, 2.54, 25.62, 1.0978, 2.41137 * 4./3., 0.47105 * 4./3.}, //11
  };

  std::vector<std::vector<double>> ex0000_half_dets = {
    {  8.3, 15.02, 6.72, 11.66, 1.3581, 4.55924, 0.57738}, // 0
    {15.02, 18.54, 3.52, 16.78, 1.0173, 2.61834, 0.440614}, // 1
    {19.97, 22.58, 2.61, 21.28, 0.9464, 3.32898, 0.502553}, // 2 from manual fit 3.32898e+00   5.02553e-01
    {22.58, 24.90, 2.32, 23.74, 0.9323, 1.9949556, 0.60642250}, // 3 from manual fit
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;
  std::vector<TString> titles;

  data = {ex1982_half_dets, ex3920_half_dets, ex3552_3630_half_dets, ex5255_5340_5375_half_dets, ex6200_half_dets, ex6200_half_dets, ex6930_half_dets, ex7100_half_dets, ex0000_half_dets};
  Ex_values = {1.982, 3.920, 3.552, 5.255, 6.200, 6.200, 6.930, 6.930, 6.930, 6.930, 6.930, 7.100, 0.000};

  /*labels = {
    {"0 #font[42]{d}_{5/2}, 2^{+}", "1 #font[42]{s}_{1/2} 2^{+}"}, // 1.982
    {"1 #font[42]{s}_{1/2}, 2^{+}", "0 #font[42]{d}_{5/2} 2^{+}"}, // 3.920
    {"0 #font[42]{d}_{5/2}, 4^{+}, 3.552", "0 #font[42]{d}_{5/2}, 0^{+}, 3.630"}, // 3.552 & 3.630
    {"0 #font[42]{d}_{5/2}, 0^{+}, 5.340", "1 #font[42]{s}_{1/2}, 2^{+}, 5.255", "1 #font[42]{s}_{1/2}, 3^{+}, 5.375"}, // triplet
    {"1 #font[42]{p}_{3/2}, 1^{-}"}, // 6200
    {"0 #font[42]{f}_{7/2}, 1^{-}"}, // 6200
    {"1 #font[42]{p}_{3/2}, 1^{-}"}, // 6930
    {"0 #font[42]{f}_{7/2}, 1^{-}"}, // 6930
    {"0 #font[42]{d}_{5/2}, 0^{+}"}, // 6930
    {"1 #font[42]{s}_{1/2}, 2^{+}"}, // 6930
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 2^{+}"}, // combined
    {"0 #font[42]{d}_{5/2}, 4^{+}"}, // 7100
    {"0 #font[42]{d}_{5/2}, 0^{+}"}  // g.s.
    };*/

  labels = {
    {"1 #font[42]{s}_{1/2} 2^{+}", "0 #font[42]{d}_{5/2}, 2^{+}"}, // 1.982
    {"1 #font[42]{s}_{1/2}, 2^{+}", "0 #font[42]{d}_{5/2} 2^{+}"}, // 3.920
    {"0 #font[42]{d}_{5/2}, 4^{+}"}, // 3.552 & 3.630
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}"}, // 5.255 & 5.340 & 5.375
    {"1 #font[42]{p}_{3/2}, 1^{-}"}, // 6200
    {"0 #font[42]{f}_{7/2}, 1^{-}"}, // 6200
    {"1 #font[42]{p}_{3/2}, 1^{-}"}, // 6930
    {"0 #font[42]{f}_{7/2}, 1^{-}"}, // 6930
    {"0 #font[42]{d}_{5/2}, 0^{+}"}, // 6930
    {"1 #font[42]{s}_{1/2}, 2^{+}"}, // 6930
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 2^{+}"}, // combined
    {"0 #font[42]{d}_{5/2}, 4^{+}"}, // 7100
    {"0 #font[42]{d}_{5/2}, 0^{+}"}  // g.s.
  };
  
  titles = {"1.982 MeV", "3.920 MeV", "3.552 & 3.630 MeV", "5.255 & 5.340 & 5.375 MeV", "6.200 MeV", "6.200 MeV", "6.930 MeV", "6.930 MeV", "6.930 MeV", "6.930 MeV", "6.930 MeV", "7.100 MeV", "0.000 MeV"};

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("../corrected_working/DWBA_17O_%s.root", pot.Data());
  TString outputDir = Form("plots_17O/minuit_extra5/ang_dist_newunc/%s", pot.Data());
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

  std::vector<std::vector<std::vector<int>>> fitPairs;

  std::cout << "Mappings " << std::endl;
  
  fitMappings = {
    {ex1982_half_dets, {0, 2}},//, 11, 12}},
    {ex3920_half_dets, {15, 13}},
    {ex3552_3630_half_dets, {36, 37}},
    {ex5255_5340_5375_half_dets, {26, 27, 28}},
    {ex6200_half_dets, {29}},
    {ex6200_half_dets, {40}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex7100_half_dets, {38}},
    {ex0000_half_dets, {39}},
  };
  fitPairs =  {
    {{0, 2}},
    {{15, 13}},
    {{36}},
    {{26, 28}},
    {{29}},
    {{40}},
    {{30}},
    {{31}},
    {{34}},
    {{35}},
    {{34, 35}},
    {{38}},
    {{39}}
  };

  std::vector<int> color = {629, 596, 418, 801, 905, 8, 9, 1};
  int ncolor = 0;

  std::ofstream resultFile;

  TString outputFilename = Form("%s/DWBA_fit_results_17O.txt", outputDir.Data());
  resultFile.open(outputFilename);

  for (size_t mappingIndex = 0; mappingIndex < fitMappings.size(); ++mappingIndex) {
    const auto& dataSet = fitMappings[mappingIndex].first;
    const auto& fitIndices = fitMappings[mappingIndex].second;

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
      std::cout << (integral_corr-integral_unc)/integral_unc << std::endl;
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
				  Form("Ex = %.3f MeV", Ex_values[mappingIndex]), 1400*2, 1400*2);

    experimentGraph->Draw("APE1 SAME");

    experimentGraph->GetHistogram()->GetXaxis()->SetRangeUser(0, 45);
    experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(min_y[mappingIndex], max_y[mappingIndex]);
    //if (mappingIndex == 2)
      //experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.4, 120);

    TAxis *axis = experimentGraph->GetXaxis();
    axis->SetLimits(0.,45.);
    
    experimentGraph->SetTitle(Form("%s;#theta_{CM} (deg);d#sigma/d#Omega (arb. units)", titles[mappingIndex].Data()));

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

    std::vector<int> bestSubset;

    std::vector<TF1*> allFitFunctions;

    for (const auto& subset : fitPairs[mappingIndex]) {
      //TF1* fitFunction = new TF1("fitFunction", combinedDWBA, 0, 60, subset.size());
      TF1* fitFunction = new TF1(Form("fitFunction_%lu", subset.size()), combinedDWBA, 0, 60, subset.size());
      for (size_t i = 0; i < subset.size(); ++i) {
	fitFunction->SetParameter(i, 0.0);
	fitFunction->SetParLimits(i, 0.0, 100.0);
      }

      nr_of_functions = subset.size();
      for(int i=0; i<nr_of_functions; i++){
        functions_ids[i] = subset[i];
      }

      TFitResultPtr fitResult = experimentGraph->Fit(fitFunction, "RLS0");

      double chi2 = 1000000.;
      double ndf = 1.0;

      if (fitResult->IsValid()) {
	chi2 = fitResult->Chi2();
	ndf = fitResult->Ndf();
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

      for (size_t comp = 0; comp < subset.size(); ++comp) {
	TF1* compFunc = new TF1(Form("fitComponent_%lu_%lu", mappingIndex, comp),
				combinedDWBA, 0, 60, subset.size());

	for (size_t i = 0; i < subset.size(); ++i) {
	  if (i == comp) compFunc->SetParameter(i, fitFunction->GetParameter(i));
	  else compFunc->SetParameter(i, 0.0);
	}

	compFunc->SetLineColor(color[(ncolor + comp + 1) % color.size()]);
	compFunc->SetLineStyle(2);
	compFunc->SetLineWidth(4);
	compFunc->Draw("L SAME");

	legend->AddEntry(compFunc, Form("%s", labels[mappingIndex][comp].Data()), "l");
      }

      

      resultFile << "========================= Fit for Ex = " << Ex_values[mappingIndex] << " MeV =========================" << std::endl;
      resultFile << "Fitted Function Indices: ";
      for (const auto& element : subset) {
	resultFile << element << " ";
      }
      resultFile << std::endl;

      resultFile << "Function Names: "  << std::endl;
      for (const auto& element : subset) {
	resultFile << graphsDWBA[element]->GetName() << std::endl;
      }
      resultFile << std::endl;

      resultFile << "Fit Parameters: ";
      for (int i = 0; i < fitFunction->GetNpar(); ++i) {
	resultFile << fitFunction->GetParameter(i) << " (" << fitFunction->GetParError(i) << ") ";
      }
      resultFile << std::endl;

      resultFile << "Chi2 (from TF1): " << chi2 << std::endl;
      resultFile << "NDF: " << ndf << std::endl;
      resultFile << "Chi2/ndf: " << (double)chi2/ndf << std::endl;
      resultFile << "------------------------------------------------------------" << std::endl << std::endl;

      ncolor++;
    }
      

    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.05);
    latex2.SetTextFont(42);          // 42 = Helvetica, 62 = bold Helvetica
    latex2.SetTextColor(kRed);
    latex2.SetTextColorAlpha(kGray+2, 0.4);
    //latex2.DrawLatex(0.25, 0.22, "PRELIMINARY");

    //experimentGraph->Draw("APE1 SAME");

    canvas->SetLogy();
    //canvas->Update();

    legend->Draw();

    experimentGraph->Draw("PE SAME");

    TString outputFilename = Form("%s/fit_Ex_%lu.png", outputDir.Data(), mappingIndex + 1);
    canvas->SaveAs(outputFilename);
    outputFilename = Form("%s/fit_Ex_%lu.pdf", outputDir.Data(), mappingIndex + 1);
    canvas->SaveAs(outputFilename);
    outputFilename = Form("%s/fit_Ex_%lu.root", outputDir.Data(), mappingIndex + 1);
    canvas->SaveAs(outputFilename);
  }
  
  resultFile.close();
}
