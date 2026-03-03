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
#include <TF1.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TSystem.h>

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};

std::vector<TString> potentials = {
  "AK", "AV", "AM", "AG", "AP",
  "HK", "HV", "HM", "HG", "HP",
  "BK", "BV", "BM", "BG", "BP",
  "DK", "DV", "DM", "DG", "DP",
  "LK", "LV", "LM", "LG", "LP",
  "QK", "QV", "QM", "QG", "QP",
  "ZK", "ZV", "ZM", "ZG", "ZP"
};

bool checking_gs = true;

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

void comp_pot_17O_2() {

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

  std::vector<std::vector<double>> ex0000_half_dets = {
    {  8.3, 18.54, 10.24, 13.42, 2.3766, 4.6937190, 0.99007398}, // 0
    {19.97, 24.90, 4.92 , 22.43, 1.8792, 5.2084505, 0.96378904}, // 1
    {25.94, 29.85, 3.91 , 27.89, 1.8279, 3.5981113, 0.98218343}, // 2
    {30.71, 34.09, 3.38 , 32.40, 1.8100, 1.9949556, 0.60642250}, // 3
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;
  std::vector<TString> titles;

  data = {ex0000_half_dets, ex1982_half_dets, ex3920_half_dets};
  Ex_values = {0.000, 1.982, 3.920};

  labels = {
    {"AK"}, {"AV"}, {"AM"}, {"AG"}, {"AP"}, {"HK"}, {"HV"}, {"HM"}, {"HG"}, {"HP"}, {"BK"}, {"BV"}, {"BM"}, {"BG"}, {"BP"}, {"DK"}, {"DV"}, {"DM"}, {"DG"}, {"DP"}, {"LK"}, {"LV"}, {"LM"}, {"LG"}, {"LP"}, {"QK"}, {"QV"}, {"QM"}, {"QG"}, {"QP"}, {"ZK"}, {"ZV"}, {"ZM"}, {"ZG"}, {"ZP"}, // 0.000
    {"AK"}, {"AV"}, {"AM"}, {"AG"}, {"AP"}, {"HK"}, {"HV"}, {"HM"}, {"HG"}, {"HP"}, {"BK"}, {"BV"}, {"BM"}, {"BG"}, {"BP"}, {"DK"}, {"DV"}, {"DM"}, {"DG"}, {"DP"}, {"LK"}, {"LV"}, {"LM"}, {"LG"}, {"LP"}, {"QK"}, {"QV"}, {"QM"}, {"QG"}, {"QP"}, {"ZK"}, {"ZV"}, {"ZM"}, {"ZG"}, {"ZP"}, // 1.982
    {"AK"}, {"AV"}, {"AM"}, {"AG"}, {"AP"}, {"HK"}, {"HV"}, {"HM"}, {"HG"}, {"HP"}, {"BK"}, {"BV"}, {"BM"}, {"BG"}, {"BP"}, {"DK"}, {"DV"}, {"DM"}, {"DG"}, {"DP"}, {"LK"}, {"LV"}, {"LM"}, {"LG"}, {"LP"}, {"QK"}, {"QV"}, {"QM"}, {"QG"}, {"QP"}, {"ZK"}, {"ZV"}, {"ZM"}, {"ZG"}, {"ZP"} // 3.920
  };
  
  titles = {"0.000 MeV", "1.982 MeV", "3.920 MeV"};

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("../corrected_working/DWBA_17O_%s.root", pot.Data());
  TString outputDir = Form("plots_17O/minuit_extra5/ang_dist_newunc_potentials/", pot.Data());
  gSystem->mkdir(outputDir, kTRUE);
  
  if (checking_gs) {
    file = TFile::Open("DWBA_17O_pot.root");
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
      graphsDWBA.push_back((TGraph*)graph->Clone());

      std::cout << "Found TGraph: " << graph->GetName() << std::endl;
    }
  }

  std::vector<std::vector<int>> fitCombos_0000 = {
    {0}, {1}, {2}, {3}, {4},
    {5}, {6}, {7}, {8}, {9},
    {10}, {11}, {12}, {13}, {14},
    {15}, {16}, {17}, {18}, {19},
    {20}, {21}, {22}, {23}, {24},
    {25}, {26}, {27}, {28}, {29},
    {30}, {31}, {32}, {33}, {34}
  };

  std::vector<std::vector<int>> fitCombos_1982 = {
    {35, 36}, {37, 38}, {39, 40}, {41, 42}, {43, 44},
    {45, 46}, {47, 48}, {49, 50}, {51, 52}, {53, 54},
    {55, 56}, {57, 58}, {59, 60}, {61, 62}, {63, 64},
    {65, 66}, {67, 68}, {69, 70}, {71, 72}, {73, 74},
    {75, 76}, {77, 78}, {79, 80}, {81, 82}, {83, 84},
    {85, 86}, {87, 88}, {89, 90}, {91, 92}, {93, 94}
  };

  std::vector<std::vector<int>> fitCombos_3920 = {
    {95, 96}, {97, 98}, {99, 100}, {101, 102}, {103, 104},
    {105, 106}, {107, 108}, {109, 110}, {111, 112}, {113, 114},
    {115, 116}, {117, 118}, {119, 120}, {121, 122}, {123, 124},
    {125, 126}, {127, 128}, {129, 130}, {131, 132}, {133, 134},
    {135, 136}, {137, 138}, {139, 140}, {141, 142}, {143, 144},
    {145, 146}, {147, 148}, {149, 150}, {151, 152}, {153, 154},
    {155, 156}, {157, 158}, {159, 160}, {161, 162}, {163, 164}
  };

  //std::vector<int> color = {629, 596, 418, 801, 905, 8, 9, 1};
  std::vector<int> color = {600, 600-4, 600-7, 600-6, 600-2,
			    416+1, 416+2, 416-2, 416-3, 416-6,
			    632+1, 632-4, 632-7, 632-3, 632+2,
			    880+1, 880+2, 880+3, 880+4, 880+5,
			    860, 860+1, 860+2, 860+3, 860+4,
			    800, 800+1, 800+2, 800+3, 800+4,
			    900, 900+1, 900+2, 900+3, 900+4};
  int ncolor = 0;

  std::ofstream resultFile;

  TString outputFilename = Form("%s/DWBA_fit_results_17O.txt", outputDir.Data());
  resultFile.open(outputFilename);

  TCanvas* canvas = new TCanvas(
				"c_dwba",
				"DWBA fits",
				1600,
				1400
				);
  canvas->Divide(2, 2);
  
  for (size_t s = 0; s < 3; ++s) {
    int cntr = 0;

    canvas->cd(s + 1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    const std::vector<std::vector<double>> *currentData = nullptr;
    const std::vector<std::vector<int>>    *currentCombos = nullptr;
    TString title = titles[s];;

    const auto &detectors = data[s];
    const auto &combos =
      (s==0) ? fitCombos_0000 :
      (s==1) ? fitCombos_1982 :
               fitCombos_3920;

    double bestChi2 = 1e9;
    std::vector<int> bestCombination;
    TF1* bestFitFunction = nullptr;
    //ncolor = 0;

    std::vector<double> theta_means;
    std::vector<double> int_corr_sins;
    std::vector<double> int_unc;
    std::vector<double> x_unc;

    for (const auto& det : data[s]) {
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
    experimentGraph->SetMarkerSize(1);
    experimentGraph->SetMarkerColor(kBlack);
    experimentGraph->SetLineColor(kBlack);
    experimentGraph->SetLineWidth(2);

    experimentGraph->Draw("AP");

    experimentGraph->GetHistogram()->GetXaxis()->SetRangeUser(0, 60);
    experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.01, 180);
    //if (mappingIndex == 2)
      //experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.4, 120);

    TAxis *axis = experimentGraph->GetXaxis();
    axis->SetLimits(0.,60.);
    
    experimentGraph->SetTitle(Form("%s;#theta_{CM} (deg);d#sigma/d#Omega (a. u.)", titles[s].Data()));

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.1);

    experimentGraph->GetXaxis()->SetTitleSize(0.044);
    experimentGraph->GetYaxis()->SetTitleSize(0.044);

    experimentGraph->GetXaxis()->SetTitleOffset(1.2);
    experimentGraph->GetYaxis()->SetTitleOffset(1.5);

    experimentGraph->GetXaxis()->SetLabelSize(0.04);
    experimentGraph->GetYaxis()->SetLabelSize(0.04);

    experimentGraph->GetXaxis()->SetLabelOffset(0.01);
    experimentGraph->GetYaxis()->SetLabelOffset(0.01);

    std::vector<int> bestSubset;

    std::vector<TF1*> allFitFunctions;

    int col = 0;

    for (const auto& subset : combos) {

      nr_of_functions = subset.size();
      for (int i=0;i<nr_of_functions;i++)
        functions_ids[i]=subset[i];

      std::cout << "After nr_of func " << nr_of_functions << std::endl;

      TF1 *f =
        new TF1("f", combinedDWBA, 0, 60, subset.size());

      std::cout << "After f " << std::endl;

      for (int p = 0; p < nr_of_functions; ++p) {
        f->SetParameter(p, 0.0);
        f->SetParLimits(p, 0, 100);
      }

      std::cout << "After pars " << std::endl;

      TFitResultPtr fitResult = experimentGraph->Fit(f, "RLS0");

       std::cout << "After fit " << std::endl;

      double chi2 = 1000000.;

      if (fitResult->IsValid()) {
	chi2 = fitResult->Chi2();
      }
      else {
	std::cout << "Invalid fit " << std::endl;
      }
      
      f->SetLineColor(color[col % color.size()]);
      f->SetLineWidth(1);
      //f->Draw("L SAME");

      const int Npoints = 200;  // resolution of the combined curve
      TGraph* combinedGraph = new TGraph(Npoints);

      for (int i = 0; i < Npoints; ++i) {
	double theta = i * 60.0 / (Npoints - 1);
	double y = f->Eval(theta);
	combinedGraph->SetPoint(i, theta, y);
      }

      if (s == 1 && cntr == 19)
	col = col + 5;
      
      combinedGraph->SetLineColor(color[col % color.size()]);
      combinedGraph->SetLineWidth(1);
      combinedGraph->Draw("L SAME");

      experimentGraph->Draw("P SAME");

      resultFile << "========================= Fit for " << title << "  =========================" << std::endl;
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
      for (int i = 0; i < f->GetNpar(); ++i) {
	resultFile << f->GetParameter(i) << " (" << f->GetParError(i) << ") ";
      }
      resultFile << std::endl;

      resultFile << "Chi2 (from TF1): " << chi2 << std::endl;
      resultFile << "------------------------------------------------------------" << std::endl << std::endl;

      col++;

     

      TLatex latex2;
      latex2.SetNDC();
      latex2.SetTextSize(0.05);
      latex2.SetTextFont(42);
      latex2.SetTextColor(kBlack);
      if (s == 0)
	latex2.DrawLatex(0.25, 0.22, "#it{l} = 2");
      if (s == 1)
	latex2.DrawLatex(0.25, 0.22, "#it{l} = 0 + 2");
      if (s == 2)
	latex2.DrawLatex(0.25, 0.22, "#it{l} = 0 + 2");
    }
      

    //TLatex lat;
    //lat.SetNDC();
    //lat.SetTextSize(0.045);
    //lat.DrawLatex(0.20,0.85,title);

    cntr++;
  }
  canvas->SetLogy();

  canvas->cd(4);
  gPad->Clear();

  TLegend *leg = new TLegend(0.15,0.1,0.95,0.9);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);

  for (size_t i=0;i<potentials.size();++i) {
    TGraph *dum = new TGraph();
    dum->SetLineWidth(3);
    dum->SetLineColor(color[i]);
    leg->AddEntry(dum,potentials[i],"l");
  }
  leg->SetNColumns(5);

  leg->Draw();

  /* ================= SAVE ================= */

  canvas->SaveAs("DWBA_17O_2x2.png");
  canvas->SaveAs("DWBA_17O_2x2.pdf");
  
  resultFile.close();
}
