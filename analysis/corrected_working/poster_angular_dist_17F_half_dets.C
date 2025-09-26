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

void poster_angular_dist_17O_half_dets() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex0937_1041_1121_half_dets = {
    {8.00 , 14.62, 6.62, 11.31, 1.2983, (13.2355 + 7.30814 + 5.04478)}, // 0
    {8.00 , 14.62, 6.62, 11.31, 1.2983, (3.94081 + 11.0333 + 5.58712)}, // 1
    {8.00 , 14.62, 6.62, 11.31, 1.2983, (1.69996 + 0 + 14.5829)}, // 2
    {14.86, 18.56, 3.70, 16.71, 1.0635, (1.42418 + 8.92131 + 8.23189)}, // 3
    {19.89, 22.62, 2.73, 21.26, 0.9886, (3.03473 + 0.173152 + 12.7338)}, // 4
    {22.76, 25.15, 2.39, 23.95, 0.9711, (3.5703 + 0. + 9.39902)}, // 5
    {26.10, 28.20, 2.10, 27.15, 0.9562, (3.27858 + 0. + 9.2654)}, // 6
    {28.30, 30.25, 1.94, 29.27, 0.9510, (1.31148 + 0. + 8.8856)}, // 7
    {31.07, 32.86, 1.78, 31.97, 0.9448, (0.449797 + 0.000258 + 2.93755) * 4./3.}, // 8
    {32.95, 34.64, 1.69, 33.80, 0.9410, (0.568829 + 0.01185 + 3.00139) * 4./3.}, // 9
    {35.32, 36.92, 1.60, 36.12, 0.9407, (3.15072) * 4./3.}, //10
    {37.00, 38.53, 1.53, 37.77, 0.9362, (0.677276) * 4./3.}, //11
  }; // bad angles

  std::vector<std::vector<double>> ex3061_half_dets = {
    {8.0  , 14.35, 6.35, 11.175, 1.23067, 5.76304}, // 2
    {14.60, 18.40, 3.80, 16.50,  1.0804,  5.7179}, // 3
    {19.76, 22.52, 2.77, 21.14,  0.9979,  4.93428}, // 4
    {22.66, 25.08, 2.42, 23.87,  0.9795,  4.64378}, // 5
    {26.04, 28.15, 2.11, 27.10,  0.9629,  3.41079}, // 6
    {28.26, 30.22, 1.96, 29.24,  0.9577,  3.23499}, // 7
    {31.05, 32.85, 1.80, 31.95,  0.9512,  1.07583 * 4./3.}, // 8
    {32.94, 34.65, 1.70, 33.80,  0.9476,  1.91486 * 4./3.}, // 9
    {35.33, 36.94, 1.61, 36.14,  0.9473,  0.76595 * 4./3.}, //10
    {37.02, 38.56, 1.54, 37.79,  0.9428,  1.63772 * 4./3.}, //11
  };

  std::vector<std::vector<double>> ex3741_3839_half_dets = {
    {14.60, 18.40, 3.80, 16.50,  1.0804,  (0.487773 + 20.1041)}, // 3
    {19.76, 22.52, 2.77, 21.14,  0.9979,  (1.75304 + 3.55357)}, // 4
    {22.66, 25.08, 2.42, 23.87,  0.9795,  (3.2814 + 1.80978)}, // 5
    {26.04, 28.15, 2.11, 27.10,  0.9629,  (1.4294 + 3.70091)}, // 6
    {28.26, 30.22, 1.96, 29.24,  0.9577,  (1.14742 + 5.24187)}, // 7
    {31.05, 32.85, 1.80, 31.95,  0.9512,  (0. + 3.77986) * 4./3.}, // 8
    {32.94, 34.65, 1.70, 33.80,  0.9476,  (0.207111 + 2.60804) * 4./3.}, // 9
    {35.33, 36.94, 1.61, 36.14,  0.9473,  (0.8015 + 2.0848) * 4./3.}, //10
    {37.02, 38.56, 1.54, 37.79,  0.9428,  (0.753033 + 1.19958) * 4./3.}, //11
  };  // bad angles

  std::vector<std::vector<double>> ex4115_half_dets = {
    {19.76, 22.52, 2.77, 21.14,  0.9979,  10.7236}, // 4
    {22.66, 25.08, 2.42, 23.87,  0.9795,  10.8906}, // 5
    {26.04, 28.15, 2.11, 27.10,  0.9629,  9.79403}, // 6
    {28.26, 30.22, 1.96, 29.24,  0.9577,  6.87214}, // 7
    {31.05, 32.85, 1.80, 31.95,  0.9512,  3.30765 * 4./3.}, // 8
    {32.94, 34.65, 1.70, 33.80,  0.9476,  3.78165 * 4./3.}, // 9
    {35.33, 36.94, 1.61, 36.14,  0.9473,  3.29471 * 4./3.}, //10
    {37.02, 38.56, 1.54, 37.79,  0.9428,  3.32874 * 4./3.}, //11
  };  // bad angles

  std::vector<std::vector<double>> ex4360_half_dets = {
    {19.76, 22.52, 2.77, 21.14,  0.9979,  1.23553}, // 4
    {22.66, 25.08, 2.42, 23.87,  0.9795,  1.0575}, // 5
    {26.04, 28.15, 2.11, 27.10,  0.9629,  0.2632}, // 6
    {28.26, 30.22, 1.96, 29.24,  0.9577,  1.25035}, // 7
    {31.05, 32.85, 1.80, 31.95,  0.9512,  0.710929 * 4./3.}, // 8
    {32.94, 34.65, 1.70, 33.80,  0.9476,  0.807541 * 4./3.}, // 9
  };  // bad angles

  std::vector<std::vector<double>> ex4652_4753_half_dets = {
    {9.60 , 15.90, 6.30, 12.75, 1.3906, (21.0052 + 0.436834)}, // 4
    {16.12, 19.64, 3.52, 17.88, 1.0821, 20.5967}, // 5
    {20.94, 23.64, 2.70, 22.29, 1.0246, 15.1561}, // 6
    {23.77, 26.16, 2.39, 24.97, 1.0100, 13.0476}, // 7
    {27.16, 29.27, 2.11, 28.21, 0.9969, 6.08755 * 4./3.}, // 8
    {29.37, 31.34, 1.96, 30.35, 0.9921, 4.54548 * 4./3.}, // 9
    {32.12, 33.93, 1.81, 33.02, 0.9873, (3.86017) * 4./3.}, //10 - I don't trust this too much
    {34.02, 35.74, 1.72, 34.88, 0.9831, (3.65818 + 0.517828) * 4./3.} //11 - I don't trust this too much
  };

  std::vector<std::vector<double>> ex4964_half_dets = {
    {8.0  , 13.94, 5.94, 10.97, 1.1304, 8.57289}, // 4
    {14.20, 18.26, 4.06, 16.23, 1.1354, 7.30343}, // 5
    {19.67, 22.57, 2.90, 21.12, 1.0441, 4.63121}, // 6
    {22.71, 25.23, 2.52, 23.97, 1.0231, 4.8697}, // 7
    {26.27, 28.46, 2.19, 27.36, 1.0074, 3.66391 * 4./3.}, // 8
    {28.57, 30.60, 2.03, 29.58, 1.0016, 2.2864 * 4./3.}, // 9
    {31.40, 33.26, 1.86, 32.33, 0.9963, 3.12153 * 4./3.} //10 - I don't trust this too much
    {33.36, 35.12, 1.77, 34.24, 0.9939, 1.73812 * 4./3.}, //11 - I don't trust this too much
  };
  
  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;

  data = {ex1982_half_dets, ex3920_half_dets, ex3552_3630_half_dets, ex5255_5340_5375_half_dets, ex6200_half_dets, ex6930_half_dets};
  Ex_values = {1.982, 3.920, 3.552, 5.255, 6.200, 6.930, 6.930, 6.930, 6.930, 6.930};
  labels = {
    {"\\ell = 0", "\\ell = 2", "\\ell = 3", "\\ell = 4"},
    {"\\ell = 0", "\\ell = 2", "\\ell = 3", "\\ell = 4"},
    {"\\ell = 2", "\\ell = 2", "\\ell = 3", "\\ell = 4"},
    {"\\ell = 2", "\\ell = 0", "\\ell = 0", "\\ell = 4"}, // 5255
    {"\\ell = 1"}, // 6200
    {"\\ell = 1"}, // 6930, 30
    {"\\ell = 3"}, // 6930, 31
    {"\\ell = 2"}, // 6930, 34
    {"\\ell = 0"}, // 6930, 35
    {"\\ell = 2", "\\ell = 0"}, // 6930, 34, 35
  };

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("DWBA_17O_%s.root", pot.Data());
  TString outputDir = Form("plots_17O/minuit_extra5/ang_dist/%s", pot.Data());
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
    {ex3552_3630_half_dets, {6}},
    {ex5255_5340_5375_half_dets, {26, 27}},
    {ex6200_half_dets, {29}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}},
    {ex6930_half_dets, {30, 31, 34, 35}}
  };
  fitPairs =  {
    {{0, 2}},
    {{15, 13}},
    {{6}},
    {{26, 27}},
    {{29}},
    {{30}},
    {{31}},
    {{34}},
    {{35}},
    {{34, 35}}
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
      resultFile << "------------------------------------------------------------" << std::endl << std::endl;

      ncolor++;
    }
      

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
