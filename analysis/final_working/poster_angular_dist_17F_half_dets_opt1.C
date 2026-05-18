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
double min_y[13] = {0.01, 0.1, 0.02, 0.9, 0.01, 0.1, 0.2, 0.3, 1, 0.2, 0.1, 0.5, 0.2};
double max_y[13] = {180,  80,  80,   80,  10,   110, 110, 8, 10, 50, 30, 10, 10};

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

void poster_angular_dist_17F_half_dets_opt1() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex0937_1041_1121_half_dets = {
    { 8.40, 15.11, 6.71, 11.75, 1.3664, (13.2355 + 7.30814 + 5.04478),  sqrt( pow( 2.9457, 2) + pow( 3.87, 2) + pow( 2.549, 2) ) }, // 0
    {15.11, 18.63, 3.52, 16.87, 1.0229, (3.94081 + 11.0333 + 5.58712),  sqrt( pow( 2.04, 2) + pow( 3.84, 2) + pow( 4.668, 2) )}, // 1
    {20.06, 22.68, 2.61, 21.37, 0.9516, (1.69996 + 0 + 14.5829),        sqrt( pow( 1.4475, 2) + pow(2.12, 2) + pow( 1.77443, 2) )}, // 2
    {22.68, 25.00, 2.32, 23.84, 0.9376, (1.42418 + 8.92131 + 8.23189),  sqrt( pow( 1.37, 2) + pow( 3.689, 2) + pow( 3.3501, 2) )}, // 3
    {26.04, 28.06, 2.03, 27.05, 0.9214, (3.03473 + 0.173152 + 12.7338), sqrt( pow( 0.92, 2) + pow( 1.7, 2) )}, // 4
    {28.06, 29.95, 1.89, 29.01, 0.9161, (3.5703 + 0. + 9.39902),        sqrt( pow( 0.72, 2) + pow( 0.98, 2) )}, // 5
    {30.82, 32.55, 1.73, 31.68, 0.9096, (3.27858 + 0. + 9.2654),        sqrt( pow( 0.6087, 2) + pow( 0.9217, 2) )}, // 6
    {32.55, 34.20, 1.65, 33.37, 0.9083, (1.31148 + 0. + 8.8856),        sqrt( pow( 0.477, 2) + pow( 0.9083, 2) )}, // 7
    {34.99, 36.54, 1.54, 35.77, 0.9028, (0.449797 + 0.000258 + 2.93755) * 4./3., (sqrt( pow( 0.267, 2) + pow( 1.158, 2) + pow( 0.517, 2))) * 4./3.}, // 8
    {36.54, 38.03, 1.49, 37.29, 0.9047, (0.568829 + 0.01185 + 3.00139) * 4./3.,  (sqrt( pow( 0.34, 2) + pow( 0.617, 2) )) * 4./3.}, // 9
    {38.71, 40.13, 1.42, 39.42, 0.9020, (3.15072) * 4./3.,  (sqrt( pow( 0.49, 2) + pow( 0., 2) )) * 4./3.}, //10
    {40.13, 41.51, 1.38, 40.82, 0.8995, (0.677276) * 4./3., (sqrt( pow( 0.23, 2) + pow( 0., 2) )) * 4./3.} //11
  }; // mean is 1.029

  std::vector<std::vector<double>> ex0937_1041_1121_half_dets_manual_fit = {
    { 8.40, 15.11, 6.71, 11.75, 1.3664, 25.2176, 1.33959}, // 0
    {15.11, 18.63, 3.52, 16.87, 1.0229, 20.3099, 1.19319}, // 1
    {20.06, 22.68, 2.61, 21.37, 0.9516, 15.7182, 1.05023}, // 2
    {22.68, 25.00, 2.32, 23.84, 0.9376, 18.5315, 1.13960}, // 3
    {26.04, 28.06, 2.03, 27.05, 0.9214, 15.9595, 1.05709}, // 4
    {28.06, 29.95, 1.89, 29.01, 0.9161, 13.3011, 0.96627}, // 5
    {30.82, 32.55, 1.73, 31.68, 0.9096, 14.4273, 1,10382}, // 6
    {32.55, 34.20, 1.65, 33.37, 0.9083, 10.2674, 0.849735}, // 7
    {34.99, 36.54, 1.54, 35.77, 0.9028, 5.01511 * 4./3., 0.626299 * 4./3.}, // 8
    {36.54, 38.03, 1.49, 37.29, 0.9047, 4.16025 * 4./3., 0.572868 * 4./3.}, // 9
    {38.71, 40.13, 1.42, 39.42, 0.9020, (3.15072) * 4./3.,  (sqrt( pow( 0.49, 2) + pow( 0., 2) )) * 4./3.}, //10
    {40.13, 41.51, 1.38, 40.82, 0.8995, (0.677276) * 4./3., (sqrt( pow( 0.23, 2) + pow( 0., 2) )) * 4./3.} //11
  }; 

  std::vector<std::vector<double>> ex3061_half_dets = {
    {8.6  , 14.48, 5.88, 11.54, 1.1763,  5.83176, 0.658813}, // 2
    {14.48, 18.32, 3.84, 16.40, 1.0839,  6.27474, 0.702328}, // 3
    {19.84, 22.59, 2.75, 21.22, 0.9967,  4.93428, 0.60}, // 4
    {22.59, 25.02, 2.43, 23.81, 0.9804,  4.64378, 0.608}, // 5
    {26.10, 28.21, 2.11, 27.15, 0.9629,  3.41079, 0.51562}, // 6
    {28.21, 30.17, 1.96, 29.19, 0.9578,  3.23499, 0.497581}, // 7
    {31.10, 32.90, 1.80, 32.00, 0.9517,  1.07583 * 4./3., 0.309092 * 4./3.}, // 8
    {32.90, 34.61, 1.71, 33.75, 0.9478,  1.91486 * 4./3., 0.382529 * 4./3.}, // 9
    {35.38, 36.98, 1.61, 36.18, 0.9476,  0.76595 * 4./3., 0.272617 * 4./3.}, //10
    {36.98, 38.52, 1.54, 37.75, 0.9427,  1.63772 * 4./3., 0.365782 * 4./3.} //11
  }; 

  std::vector<std::vector<double>> ex3724_3839_half_dets = { // 3.79
    { 8.20, 14.62, 6.42, 11.41, 1.2701,  (0.487773 + 20.1041), sqrt( pow(5.06347, 2) + pow(3.3762, 2) )}, // 3
    {16.64, 20.00, 3.36, 18.32, 1.0560,  (1.75304 + 3.55357),  sqrt( pow(0.561044, 2) + pow(0.691994, 2) )}, // 4
    {20.00, 22.77, 2.77, 21.39, 1.0112,  (3.2814 + 1.80978),   sqrt( pow(0.8212, 2) + pow(0.823344, 2) )}, // 5
    {23.98, 26.31, 2.32, 25.15, 0.9863,  (1.4294 + 3.70091),   sqrt( pow(0.556658, 2) + pow(0.705114, 2) )}, // 6
    {26.31, 28.44, 2.13, 27.37, 0.9791,  (1.14742 + 5.24187),  sqrt( pow(0.6758, 2) + pow(0.885512, 2) )}, // 7
    {29.44, 31.36, 1.92, 30.40, 0.9715,  (0. + 3.77986) * 4./3.,       (sqrt( pow(0.5354, 2) + pow(0., 2) )) * 4./3.}, // 8
    {31.36, 33.17, 1.81, 32.26, 0.9684,  (0.207111 + 2.60804) * 4./3., (sqrt( pow(0.310743, 2) + pow(0.51939, 2) )) * 4./3.}, // 9
    {33.99, 35.67, 1.68, 34.83, 0.9624,  (0.8015 + 2.0848) * 4./3.,    (sqrt( pow(0.526342, 2) + pow(0.620528, 2) )) * 4./3.}, //10
    {35.67, 37.29, 1.62, 36.48, 0.9622,  (0.753033 + 1.19958) * 4./3., (sqrt( pow(0.58534, 2) + pow(0.668129, 2) )) * 4./3.} //11
  };  

  std::vector<std::vector<double>> ex4115_half_dets = {
    {14.84, 18.65, 3.81, 16.74, 1.0972,  10.7236, 0.908082}, // 4
    {18.65, 21.65, 3.00, 20.15, 1.0339,  10.8906, 0.933796}, // 5
    {22.93, 25.38, 2.44, 24.16, 0.9992,  9.79403, 0.865031}, // 6
    {25.38, 27.60, 2.22, 26.49, 0.9901,  6.87214, 0.761571}, // 7
    {28.63, 30.62, 1.98, 29.63, 0.9807,  3.30765 * 4./3., 0.522647 * 4./3.}, // 8
    {30.62, 32.48, 1.87, 31.55, 0.9765,  3.78165 * 4./3., 0.550127 * 4./3.}, // 9
    {33.32, 35.06, 1.73, 34.19, 0.9737,  3.29471 * 4./3., 0.518427 * 4./3.}, //10
    {35.06, 36.71, 1.65, 35.88, 0.9692,  3.32874 * 4./3., 0.690289 * 4./3.} //11
  }; 

  std::vector<std::vector<double>> ex4360_half_dets = {
    {13.15, 17.52, 4.36, 15.34, 1.1544,  1.23553, 0.397278}, // 4
    {17.52, 20.74, 3.22, 19.13, 1.0544,  1.0575,  0.5129}, // 5
    {22.09, 24.64, 2.55, 23.36, 1.0100,  0.2632,  0.259}, // 6
    {24.64, 26.94, 2.30, 25.79, 1.0003,  1.25035, 0.431587}, // 7
    {28.01, 30.04, 2.04, 29.02, 0.9877,  0.710929 * 4./3., 0.29516 * 4./3.}, // 8
    {30.04, 31.95, 1.91, 31.00, 0.9835,  0.807541 * 4./3., 0.3157 * 4./3.} // 9
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
    { 8.20, 15.72, 7.52, 11.96, 1.5583, (21.0052 + 0.436834), 1.6249}, // 4
    {15.72, 19.35, 3.63, 17.54, 1.0922, 20.5967, 1.27427}, // 5
    {20.82, 23.54, 2.72, 22.18, 1.0274, 15.1561, 1.05519}, // 6
    {23.54, 25.96, 2.42, 24.75, 1.0139, 13.0476, 1.00247}, // 7
    {27.08, 29.19, 2.12, 28.13, 0.9982, 6.08755 * 4./3., 0.673348 * 4./3.}, // 8
    {29.19, 31.17, 1.98, 30.18, 0.9937, 4.54548 * 4./3., 0.596529 * 4./3.}, // 9
    {32.05, 33.87, 1.82, 32.96, 0.9879, (3.86017) * 4./3., 0.562701 * 4./3.}, //10 - I don't trust this too much
    {33.87, 35.60, 1.73, 34.73, 0.9849, (3.65818 + 0.517828) * 4./3., (sqrt( pow(1.27, 2) + pow(1.17135, 2) )) * 4./3.} //11 - I don't trust this too much
  };

  std::vector<std::vector<double>> ex4964_half_dets = {
    {8.2  , 14.07, 5.87, 11.13, 1.1336, 8.57289, 0.847075}, // 4
    {14.07, 18.17, 4.10, 16.12, 1.1386, 7.30343, 0.783341}, // 5
    {19.75, 22.64, 2.89, 21.20, 1.0432, 4.63121, 0.603035}, // 6
    {22.64, 25.17, 2.53, 23.90, 1.0240, 4.8697,  0.637558}, // 7
    {26.32, 28.51, 2.19, 27.42, 1.0075, 3.66391 * 4./3., 0.528969 * 4./3.}, // 8
    {28.51, 30.54, 2.03, 29.53, 1.0017, 2.2864 * 4./3.,  0.430266 * 4./3.}, // 9
    {31.45, 33.31, 1.86, 32.38, 0.9961, 3.12153 * 4./3., 0.515049 * 4./3.}, //10 - I don't trust this too much
    {33.31, 35.08, 1.77, 34.19, 0.9942, 1.73812 * 4./3., 0.51603 * 4./3.} //11 - I don't trust this too much
  };

  std::vector<std::vector<double>> ex000_full_dets = { // fitted manually in TBrowser for full dets
    {16.21, 22.07, 5.86, 19.14, 1.9224, 4.668,     0.6991701}, // 0
    {23.23, 27.47, 4.24, 25.35, 1.8172, 3.1347672, 0.46539170}, // 1
    {28.40, 31.96, 3.56, 30.18, 1.7898, 2.3449830, 0.41363009}, // 2
    {32.76, 35.91, 3.15, 34.33, 1.7752, 1.8110064, 0.36372493}, // 3
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;
  std::vector<TString> titles;

  data = {ex0937_1041_1121_half_dets, ex4964_half_dets, ex3061_half_dets, ex4652_4753_half_dets, ex4360_half_dets, ex4115_half_dets, ex3724_3839_half_dets, ex000_full_dets};
  Ex_values = {1.041, 4.964, 3.062, 4.652, 4.360, 4.115, 3.790, 0.000};
  labels = {
    {"1 #font[42]{s}_{1/2}, 3^{+}", "0 #font[42]{d}_{5/2}, 5^{+}"}, // 0.937, 1.042, 1.121
    {"1 #font[42]{s}_{1/2}, 2^{+}", "0 #font[42]{d}_{5/2}, 2^{+}"}, // 4.964
    {"1 #font[42]{s}_{1/2}, 2^{+}", "0 #font[42]{d}_{5/2}, 2^{+}"}, // 3.061
    {"0 #font[42]{d}_{5/2}, 4^{+}"}, // 4.652, 4.753
    {"0 #font[42]{d}_{5/2}, 1^{+}"}, // 4.360
    {"1 #font[42]{s}_{1/2}, 3^{+}", "0 #font[42]{d}_{5/2}, 2^{+}"}, // 4.115
    {"1 #font[42]{s}_{1/2}, 2^{+}", "0 #font[42]{d}_{5/2}, 2^{+}"}, // 3.724, 3.839
    {"0 #font[42]{d}_{3/2}, 1^{+}"}, // 0.000
  };
  

  titles = {"0.937 & 1.042 & 1.121 MeV", "4.964 MeV", "3.062 MeV", "4.652 & 4.753 MeV", "4.360 MeV", "4.115 MeV", "3.724 & 3.839 MeV", "0.000 MeV"};

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("../corrected_working/DWBA_17F_%s.root", pot.Data());
  TString outputDir = Form("plots_17F/minuit_new/ang_dist_newunc/%s", pot.Data());
  gSystem->mkdir(outputDir, kTRUE);
  
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
    {ex0937_1041_1121_half_dets, {26, 28}},
    {ex4964_half_dets, {2, 3}},//, 11, 12}},
    {ex3061_half_dets, {5, 6}},
    {ex4652_4753_half_dets, {29}},
    {ex4360_half_dets, {21}},
    {ex4115_half_dets, {20, 19}},
    {ex3724_3839_half_dets, {24, 23}},
    {ex000_full_dets, {39}}
  };
  fitPairs =  {
    {{26, 28}},
    {{2, 3}},
    {{5, 6}},
    {{29}},
    {{21}},
    {{20, 19}},
    {{24, 23}},
    {{39}}
  };

  std::vector<int> color = {629, 596, 418, 801, 905, 8, 9, 1};
  int ncolor = 0;

  std::ofstream resultFile;

  TString outputFilename = Form("%s/DWBA_fit_results_17F.txt", outputDir.Data());
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
