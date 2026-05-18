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

void poster_angular_dist_17Fw17Orecoils_half_dets_opt1_sd() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

   std::vector<std::vector<double>> ex6136_6163_half_dets = {
    {13.06, 17.60, 4.54, 15.33, 1.2010, 4.01808 + 0., sqrt( pow(0.5717, 2) + pow(0.6187, 2) )}, // 6
    {17.60, 20.94, 3.33, 19.27, 1.1003, 1.87111 + 0., sqrt( pow(0.430758, 2) + pow(0.21, 2) )}, // 7
    {22.37, 25.01, 2.63, 23.69, 1.0579, (3.28537 + 0.) * 4./3., sqrt( pow(0.5086, 2) + pow(0.4576, 2) )* 4./3.}, // 8
    {25.01, 27.37, 2.36, 26.19, 1.0415, (4.51246 + 0.) * 4./3., sqrt( pow(0.592214, 2) + pow(0.4192, 2) )* 4./3.}, // 9
    {28.40, 30.50, 2.10, 29.45, 1.0331, (5.23575 + 0.) * 4./3., sqrt( pow(0.6359, 2) + pow(1.95704, 2) )* 4./3.}, //10
    {30.50, 32.47, 1.97, 31.49, 1.0287, (4.94403 + 0.) * 4./3., sqrt( pow(0.6184, 2) + pow(0.845107, 2) )* 4./3.}, //11
  }; // thetaCM calculated for (6.136 + 6.163) / 2.0 = 6.1495
  
  std::vector<std::vector<double>> ex6633_6643_half_dets = {
    { 8.00, 14.66, 6.66, 11.33, 1.3084, 0.00013 + 14.4628,              sqrt( pow(1.03, 2) + pow(5.09, 2) )}, // 6
    {14.66, 18.73, 4.06, 16.70, 1.1670, 0. + 9.81058,                   sqrt( pow(2.13, 2) + pow(0.85, 2) )}, // 7
    {20.37, 23.28, 2.91, 21.83, 1.0834, (0. + 3.84624) * 4./3.,         sqrt( pow(1.5904, 2) +  pow(0.5425,  2))* 4./3.}, // 8
    {23.28, 25.85, 2.56, 24.57, 1.0653, (0. + 4.66187) * 4./3.,         sqrt( pow(1.66424, 2) + pow(0.6285,  2))* 4./3.}, // 9
    {26.96, 29.19, 2.23, 28.07, 1.0497, (0.00018577 + 3.37614) * 4./3., sqrt( pow(2.21755, 2) + pow(0.53759, 2))* 4./3.}, //10
    {29.19, 31.26, 2.08, 30.23, 1.0447, (0.0 + 3.12785)  * 4./3.,       sqrt( pow(1.23134, 2) + pow(0.512909,2))* 4./3.}, //11
  }; // thetaCM calculated for (6.633 + 6.643) / 2. = 6.638

  std::vector<std::vector<double>> ex5600_full_dets = {
    { 8.2, 14.63,  6.43, 11.42, 1.2725, 3.42311, 0.68006}, // 2
    {16.71, 23.05, 6.34, 19.88, 2.1557, 1.49325, 0.349}, // 3
    {24.33, 28.92, 4.59, 26.63, 2.0559, 0.79 * 4./3., 0.27} // 4
  };    


  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;
  std::vector<std::vector<TString>> labels;
  std::vector<TString> titles;

  data = {ex6136_6163_half_dets, ex6136_6163_half_dets, ex6136_6163_half_dets, ex6136_6163_half_dets,
	  ex6136_6163_half_dets, ex6136_6163_half_dets, ex6136_6163_half_dets, ex6136_6163_half_dets,
	  ex6136_6163_half_dets, ex6136_6163_half_dets, ex6136_6163_half_dets, ex6136_6163_half_dets,
	  ex6633_6643_half_dets, ex5600_full_dets, ex5600_full_dets};
  Ex_values = {6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.136, 6.633, 5.600, 5.600};
  labels = {
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"}, // 6136_6163 will probably be IAS of thr 5.34 (d) and 5.38 (s) states in O; there might be also addition of 6.108 state in F, I'm not sure what l is that, but also 6.108 is T=0 state so it's weaker
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"0 #font[42]{d}_{5/2}, 0^{+}", "1 #font[42]{s}_{1/2}, 3^{+}", "1 #font[42]{p}_{3/2}, 3^{-}"},
    {"1 #font[42]{f}_{7/2}, 2^{-}"}, // 6.633 state - parity and T unknown, 6.643 - 2- state but T=1 so stronger? I don't know which state could it possibly correspond to in oxygen - maybe the 6.35 state? It's also 2-, but very weak
    {"1 #font[42]{p}_{3/2}, 1^{-}"},
    {"1 #font[42]{s}_{3/2}, 2^{+}"}
  };
  titles = {"6.136 & 6.163 MeV", "6.136 & 6.163 MeV - 0.25", "6.136 & 6.163 MeV - 0.5", "6.136 & 6.163 MeV - 0.75",
	    "6.136 & 6.163 MeV - 1.", "6.136 & 6.163 MeV - 1.25", "6.136 & 6.163 MeV - 1.5", "6.136 & 6.163 MeV - 1.75",
	    "6.136 & 6.163 MeV - 2", "6.136 & 6.163 MeV - 2.25", "6.136 & 6.163 MeV - 2.5", "6.136 & 6.163 MeV",
	    "6.643 MeV", "5.6 MeV", "5.6 MeV"};

  std::vector<TGraph*> graphs;

  TFile *file;

  const auto& pot = potentials[0];

  TString filename = Form("../corrected_working/DWBA_17F_%s_he.root", pot.Data());
  TString outputDir = Form("plots_17Fw17Orecoils/minuit/ang_dist_newunc_sd/%s", pot.Data());
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
    {ex6136_6163_half_dets, {0, 1}},
    {ex6136_6163_half_dets, {7, 8}},
    {ex6136_6163_half_dets, {10, 11}},
    {ex6136_6163_half_dets, {13, 14}},
    {ex6136_6163_half_dets, {16, 17}},
    {ex6136_6163_half_dets, {19, 20}},
    {ex6136_6163_half_dets, {22, 23}},
    {ex6136_6163_half_dets, {25, 26}},
    {ex6136_6163_half_dets, {28, 29}},
    {ex6136_6163_half_dets, {31, 32}},
    {ex6136_6163_half_dets, {34, 35}},
    {ex6136_6163_half_dets, {37, 38}},
    {ex6633_6643_half_dets, {2}},
    {ex5600_full_dets, {40}},
    {ex5600_full_dets, {41}}
  };
  fitPairs =  {
    {{0, 1}},
    {{7, 8}},
    {{10, 11}},
    {{13, 14}},
    {{16, 17}},
    {{19, 20}},
    {{22, 23}},
    {{25, 26}},
    {{28, 29}},
    {{31, 32}},
    {{34, 35}},
    {{37, 38}},
    {{2}},
    {{40}},
    {{41}}
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

    experimentGraph->GetHistogram()->GetXaxis()->SetRangeUser(0, 45);
    experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.1, 100);
    if (mappingIndex == 2)
      experimentGraph->GetHistogram()->GetYaxis()->SetRangeUser(.4, 120);

    TAxis *axis = experimentGraph->GetXaxis();
    axis->SetLimits(0.,45.);
    
    //experimentGraph->SetTitle(Form("Ex = %.3f MeV;#theta_{CM} (deg);d#sigma/d#Omega (a. u.)", Ex_values[mappingIndex]));
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

    canvas->SetLogy();
    //canvas->Update();

    legend->Draw();
    experimentGraph->Draw("PE SAME");

    TString outputFilename = Form("%s/fit_Ex_%lu.png", outputDir.Data(), mappingIndex + 1);
    canvas->SaveAs(outputFilename);
    outputFilename = Form("%s/fit_Ex_%lu.pdf", outputDir.Data(), mappingIndex + 1);
    canvas->SaveAs(outputFilename);
  }
  
  resultFile.close();
}
