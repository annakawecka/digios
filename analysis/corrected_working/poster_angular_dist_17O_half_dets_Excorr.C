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
  "AV", "HV", "HM", "HG", "HP", "BK", "BV", "BM", "BG", "BP",
  "DK", "DV", "DM", "DG", "DP", "QK", "QV", "QM", "QG", "QP",
  "ZK", "ZV", "ZM", "ZG", "ZP"
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

void poster_angular_dist_17O_half_dets_Excorr() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex3061_half_dets = {
    {8.0  , 14.35, 6.35, 11.175, 1.23067, 5.81208}, // 2
    {14.60, 18.40, 3.80, 16.50,  1.0804,  5.27457}, // 3
    {19.76, 22.52, 2.77, 21.14,  0.9979,  4.93428}, // 4
    {22.66, 25.08, 2.42, 23.87,  0.9795,  4.64378}, // 5
    {26.04, 28.15, 2.11, 27.10,  0.9629,  3.48367}, // 6
    {28.26, 30.22, 1.96, 29.24,  0.9577,  3.27969}, // 7
    {31.05, 32.85, 1.80, 31.95,  0.9512,  1.09372 * 4./3.}, // 8
    {32.94, 34.65, 1.70, 33.80,  0.9476,  1.8835  * 4./3.}, // 9
    {35.33, 36.94, 1.61, 36.14,  0.9473,  1.32144 * 4./3.}, //10
    {37.02, 38.56, 1.54, 37.79,  0.9428,  2.04376 * 4./3.}, //11
  };

  std::vector<std::vector<double>> ex4964_half_dets = {
    {8.0  , 13.94, 5.94, 10.97, 1.1304, 8.57289}, // 4
    {14.20, 18.26, 4.06, 16.23, 1.1354, 7.30343}, // 5
    {19.67, 22.57, 2.90, 21.12, 1.0441, 4.46986}, // 6
    {22.71, 25.23, 2.52, 23.97, 1.0231, 4.04592}, // 7
    {26.27, 28.46, 2.19, 27.36, 1.0074, 3.61657 * 4./3.}, // 8
    {28.57, 30.60, 2.03, 29.58, 1.0016, 2.27792 * 4./3.}, // 9
    {31.40, 33.26, 1.86, 32.33, 0.9963, 2.58783 * 4./3.} //10 - I don't trust this too much
    //{33.36, 35.12, 1.77, 34.24, 0.9939, 0.770579 * 4./3.}, //11 - I don't trust this too much
  };

  std::vector<std::vector<double>> ex4652_4753_half_dets = {
    {9.60 , 15.90, 6.30, 12.75, 1.3906, (21.0052 + 0.436834)}, // 4
    {16.12, 19.64, 3.52, 17.88, 1.0821, 20.5967}, // 5
    {20.94, 23.64, 2.70, 22.29, 1.0246, 15.3903}, // 6
    {23.77, 26.16, 2.39, 24.97, 1.0100, 13.4816}, // 7
    {27.16, 29.27, 2.11, 28.21, 0.9969, 5.83519 * 4./3.}, // 8
    {29.37, 31.34, 1.96, 30.35, 0.9921, 4.5355 * 4./3.}, // 9
    {32.12, 33.93, 1.81, 33.02, 0.9873, (3.76639 + 0.727264) * 4./3.}, //10 - I don't trust this too much
    //{34.02, 35.74, 1.72, 34.88, 0.9831, 5.64157 * 4./3.} //11 - I don't trust this too much
  };

  std::vector<std::vector<double>> ex1982_half_dets = {
    {8.00 , 14.62, 6.62, 11.31, 1.2983, 28.3942}, // 2
    {14.86, 18.56, 3.70, 16.71, 1.0635, 17.1773}, // 3
    {19.89, 22.62, 2.73, 21.26, 0.9886, 16.2642}, // 4
    {22.76, 25.15, 2.39, 23.95, 0.9711, 12.7089}, // 5
    {26.10, 28.20, 2.10, 27.15, 0.9562, 10.4}, // 6
    {28.30, 30.25, 1.94, 29.27, 0.9510, 10.5249}, // 7
    {31.07, 32.86, 1.78, 31.97, 0.9448, 5.96941 * 4./3.}, // 8
    {32.95, 34.64, 1.69, 33.80, 0.9410, 4.00472 * 4./3.}, // 9
    {35.32, 36.92, 1.60, 36.12, 0.9407, 3.44069 * 4./3.}, //10
    {37.00, 38.53, 1.53, 37.77, 0.9362, 2.11812 * 4./3.}, //11
  };

  std::vector<std::vector<double>> ex3920_half_dets = {
    { 8.00, 14.08, 6.08, 11.04, 1.1643, 26.7301}, // 4
    {14.33, 18.32, 4.00, 16.33, 1.1236, 22.9592}, // 5
    {19.72, 22.59, 2.87, 21.15, 1.0362, 14.8715}, // 6 // from automatic fitting it was 13.9381, from manual fitting 14.871480
    {22.73, 25.23, 2.50, 23.98, 1.0159, 15.9390}, // 7 // from automatic fitting it was 14.6462, from manual 15.938928
    {26.26, 28.44, 2.18, 27.35, 1.0007, 9.23547 * 4./3.}, // 8 // from automatic fitting it was 8.40969, from manual 
    {28.55, 30.56, 2.02, 29.56, 0.9950, 9.04759 * 4./3.}, // 9 // automatic fitting looked ok
    {31.36, 33.22, 1.85, 32.29, 0.9899, 7.46738 * 4./3.}, //10 // not to be super trusted
    {33.31, 35.07, 1.76, 34.19, 0.9877, 6.32486 * 4./3.}, //11 // not to be super trusted
  };

  std::vector<std::vector<double>> ex3552_3630_half_dets = {
    {10.92, 16.29, 5.37, 13.60, 1.2636, 63.0214 + 7.39285}, // 4
    {16.49, 19.92, 3.42, 18.20, 1.0694, 66.6946}, // 5
    {21.17, 23.83, 2.66, 22.50, 1.0161, 48.6652}, // 6 // from automatic fitting 49.5023, from manual 
    {23.96, 26.31, 2.35, 25.14, 0.9992, 37.9169}, // 7
    {27.29, 29.38, 2.08, 28.33, 0.9889, 20.050 * 4./3.}, // 8
    {29.48, 31.42, 1.94, 30.45, 0.9840, 15.5393 * 4./3.}, // 9
    {32.20, 33.99, 1.80, 33.09, 0.9803, 8.8607956 * 4./3.}, //10 // not to be super trusted, from automatic fitting 11.5826, from manual 8.8607956
    //{34.08, 35.79, 1.70, 34.94, 0.9760,  * 4./3.}, //11 // not to be super trusted
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;

  data = {ex1982_half_dets, ex3920_half_dets, ex3552_3630_half_dets};
  Ex_values = {1.982, 3.920, 3.552};
  labels = {
    {"\\ell = 0", "\\ell = 2", "\\ell = 3", "\\ell = 4"},
    {"\\ell = 0", "\\ell = 2", "\\ell = 3", "\\ell = 4"},
    {"\\ell = 2", "\\ell = 2", "\\ell = 3", "\\ell = 4"}
  };

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("DWBA_17O_%s.root", pot.Data());
  TString outputDir = Form("plots_17O/minuit_Excorr/ang_dist/%s", pot.Data());
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
    {ex3920_half_dets, {13, 15}},
    {ex3552_3630_half_dets, {10}},
  };
  fitPairs =  {
    {{0, 2}},
    {{13, 15}},
    {{10}},
  };

  std::vector<int> color = { 1, 629, 596, 418, 801, 905, 8, 9};
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

    /*for (size_t j = 0; j < fitIndices.size(); ++j) {
      int fitIndex = fitIndices[j];
      TGraph* fitGraph = graphsDWBA[fitIndex];
      fitGraph->SetLineColor(color[ncolor]);
      fitGraph->SetLineWidth(4);
      fitGraph->Draw("L");

      ncolor++;

      legend->AddEntry(fitGraph, Form("%s", labels[mappingIndex][j].Data()), "l");
      }*/

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

      fitFunction->SetLineColor(color[ncolor]);
      fitFunction->SetLineWidth(4);
      fitFunction->Draw("L SAME");

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
