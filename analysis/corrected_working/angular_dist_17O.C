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

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};

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

void angular_dist_17O() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex1982 = {
    {7,		18.48,	11.48,	12.74,	2.53165228,	29.893, 1.886},
    {19.97,	25.09,	5.12,	22.53,	1.961815663,	19.140, 1.406},
    {26.16,	30.2,	4.04,	28.18,	1.907862135,	13.476, 1.237},
    {31.12,	34.6,	3.48,	32.86,	1.888206768,	9.113,  1.095},
    {35.37,	38.49,	3.12,	36.93,	1.874617232,	4.736,  0.938}
  }; // Theta min, theta max, delta theta, theta mean, sin(x)dx * 180/pi, integral corr, integral corr error

  /*std::vector<std::vector<double>> ex3552 = {
    {11.14,	19.83,	8.69,	15.485,	2.320109111,	26.5151},
    {21.25,	26.25,	5,	23.75,	2.013733447,	48.4343},
    {27.35,	31.37,	4.02,	29.36,	1.970987554,	35.6410},
    {32.24,	35.74,	3.5,	33.99,	1.9566687,	17.1285}
  };

  std::vector<std::vector<double>> ex3630 = {
    {10.07,	19.51,	9.44,	14.79,	2.409814986,	106.918    },
    {20.96,	26.03,	5.07,	23.495,	2.021252025,	36.4969    },
    {27.14,	31.2,	4.06,	29.17,	1.978854275,	11.09661333},
    {32.07,	35.59,	3.52,	33.83,	1.959691857,	11.15053333}
    };*/

  std::vector<std::vector<double>> ex3920 = {
    {7,		18.24,	11.24,	12.62,	2.455758891,	32.820, 2.176},
    {19.8,	25.17,	5.37,	22.485,	2.053711111,	20.230, 1.543},
    {26.32,	30.51,	4.19,	28.415,	1.993830287,	15.465, 1.564},
    {31.41,	35.02,	3.61,	33.215,	1.977493989,	11.746, 1.473}
  };

  std::vector<std::vector<double>> ex5255 = {
    {11.93,	20.37,	8.44,	16.15,	2.347611196,	24.781, 1.152},
    {21.86,	26.94,	5.08,	24.4,	2.098570501,	20.298, 1.696},
    {27.99,	32.09,	4.1,	30.04,	2.052478357,	24.758, 1.453}
  };

  std::vector<std::vector<double>> ex3530 = {
    {11.38,	19.93,	8.55,	15.65,	2.30645002,	89.22,  3.713 },
    {21.33,	26.32,	4.99,	23.82,	2.01528462,	55.403, 2.331},
    {27.41,	31.42,	4.01,	29.42,	1.972,		30.692, 2.216},
    {32.29,	35.79,	3.5,	34.04,	1.9552,		17.2,   1.818}
  };

   std::vector<std::vector<double>> ex0000 = {
     {8.00,	18.54,	10.54,	13.27,	2.4193, 4.097, 0.821},
     {19.97,	24.90,	4.92,	22.43,	1.8792,	2.961, 0.647},
     {25.94,	29.85,	3.91,	27.89,	1.8279,	1.733, 0.488},
     {30.71,	34.09,	3.38,	32.40,	1.8100,	0.534, 0.195}
   };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;

  if (checking_gs) {
    data = {ex5255, ex5255, ex5255, ex5255, ex5255, ex5255, ex5255};
    Ex_values = {5.255, 5.255, 5.255, 5.255, 5.255, 5.255, 5.255};
  } else {
    data = {ex1982, ex3530, ex3920, ex5255};
    Ex_values = {1.982, 3.553, 3.920, 5.255};
  }

  std::vector<TGraph*> graphs;

  //TFile *file = TFile::Open("../working/DWBA_17O_more_states.root");
  TFile *file = TFile::Open("DWBA_17O_pot_5255.root");
  
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

  if (checking_gs) {
    fitMappings = {
      {ex5255, {0, 7, 14, 21, 28}},
      {ex5255, {1, 8, 15, 22, 29}},
      {ex5255, {2, 9, 16, 23, 30}},
      {ex5255, {3, 10, 17, 24, 31}},
      {ex5255, {4, 11, 18, 25, 32}},
      {ex5255, {5, 12, 19, 26, 33}},
      {ex5255, {6, 13, 20, 27, 34}},
    };
    fitPairs =  {
      {{0}, {7}, {14}, {21}, {28}},
      {{1}, {8}, {15}, {22}, {29}},
      {{2}, {9}, {16}, {23}, {30}},
      {{3}, {10}, {17}, {24}, {31}},
      {{4}, {11}, {18}, {25}, {32}},
      {{5}, {12}, {19}, {26}, {33}},
      {{6}, {13}, {20}, {27}, {34}},
    };
  } else {
    fitMappings = {
      {ex1982, {0, 1, 2, 3, 4}},
      //{ex3552, {5, 6, 7, 8, 9}},
      //{ex3630, {10, 11, 12}},
      {ex3920, {13, 14, 15, 16, 17, 18}},
      {ex5255, {19, 20, 21, 22}},
      {ex3530, {5, 6, 7, 8, 9, 10, 11, 12}},
    };
    fitPairs =  {
      {{0, 2, 3}, {0, 2}, {0, 3}, {1, 4}},
      //{{5}, {6}, {5, 6}, {8, 9}},
      //{{10}, {11}, {12}, {10, 11}},
      {{13}, {17, 18}, {15, 16}, {14, 15}},
      {{19}, {20, 21}, {19, 22}},
      {{5}, {10}, {5, 6}, {8, 9}},
    };
  }

  std::cout << "After Mappings " << std::endl;

  std::vector<int> color = { 1, 629, 596, 418, 801, 905, 8, 9};
  int ncolor = 0;

  std::ofstream resultFile;

  if (checking_gs) {
    resultFile.open("plots_17O/ang_dist/checking_5255_DWBA_fit_results_17O.txt");
  } else {
    resultFile.open("plots_17O/ang_dist/DWBA_fit_results_17O.txt");
  }

  for (size_t mappingIndex = 0; mappingIndex < fitMappings.size(); ++mappingIndex) {
    std::cout << "Inside loop " << std::endl;
    const auto& dataSet = fitMappings[mappingIndex].first;
    const auto& fitIndices = fitMappings[mappingIndex].second;
    //auto allCombinations = generateCombinations(fitIndices);

    std::cout << "After combinations " << std::endl;

    double bestChi2 = 1e9;
    std::vector<int> bestCombination;
    TF1* bestFitFunction = nullptr;

    ncolor = 0;

    std::cout << "===========================================================================" << std::endl;

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
      int_unc.push_back(integral_unc / sin_x_dx);
      //int_unc.push_back(0.1);
      x_unc.push_back(0.71);
    }

    TGraphErrors* experimentGraph = new TGraphErrors(theta_means.size(), &theta_means[0], &int_corr_sins[0], &x_unc[0], &int_unc[0]);
    experimentGraph->SetTitle("");
    experimentGraph->SetMarkerStyle(20);
    experimentGraph->SetMarkerSize(1);
    experimentGraph->SetMarkerColor(kBlack);
    experimentGraph->SetLineColor(kBlack);
    experimentGraph->SetLineWidth(2);
    //experimentGraph->SetErrorSize(1.5);

    experimentGraph->GetXaxis()->SetRangeUser(0, 60);
    experimentGraph->GetYaxis()->SetRangeUser(.01, 100);

    TCanvas* canvas = new TCanvas(Form("fit_canvas_%lu", mappingIndex + 1), 
				  Form("Ex = %.3f MeV", Ex_values[mappingIndex]), 1600, 1200);

    experimentGraph->Draw("APE1");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(22);
    latex.DrawLatexNDC(0.5, 0.92, Form("Plot Ex = %.3f MeV", Ex_values[mappingIndex]));

    TLegend* legend = new TLegend(0.12, 0.15, 0.72, 0.35);
    legend->SetTextSize(0.03);
    legend->SetMargin(0.1);
    legend->SetBorderSize(0);
    legend->SetFillColor(0); 
    //legend->AddEntry(experimentGraph, "Data Points", "p");

    for (int fitIndex : fitIndices) {
      std::cout << "dataSet: " << mappingIndex << " fitIndex: " << fitIndex << std::endl;
      TGraph* fitGraph = graphsDWBA[fitIndex];
      fitGraph->SetLineColor(color[ncolor]);
      fitGraph->SetLineWidth(2);
      fitGraph->Draw("L");

      ncolor++;

      legend->AddEntry(fitGraph, Form("Func %d: %s", fitIndex, fitGraph->GetName()), "l");

    }
    bestChi2 = 1000;
    for (const auto& subset : fitPairs[mappingIndex]) {
      TF1* fitFunction = new TF1("fitFunction", combinedDWBA, 0, 60, subset.size());
      for (size_t i = 0; i < subset.size(); ++i) {
	fitFunction->SetParameter(i, 0.0);
	fitFunction->SetParLimits(i, 0.0, 100.0);
      }

      nr_of_functions = subset.size();
      for(int i=0; i<nr_of_functions; i++){
        functions_ids[i] = subset[i];
      }

      experimentGraph->Fit(fitFunction, "RL0");

      double chi2 = fitFunction->GetChisquare();
      int ndf = fitFunction->GetNDF();

      double x, y, xe, ye;
      double chi2ndf = 0;
      for(int i=0; i< experimentGraph->GetN(); i++){
        experimentGraph->GetPoint(i, x, y);
	xe = experimentGraph->GetErrorX(i);
	ye = experimentGraph->GetErrorY(i);
        cout<<"x: "<<x<<", y: "<<y<<", func: "<<fitFunction->Eval(x)<<", diff: "<<y - fitFunction->Eval(x)<<endl;
        chi2ndf += (y - fitFunction->Eval(x)) * (y - fitFunction->Eval(x)) / ye / ye;
      }
      //double chi2ndf = (ndf > 0) ? (chi2 / ndf) : chi2;
      //double chi2ndf = chi2 ;
      cout<<"Chi2 from calculation: "<<chi2ndf<<", chi2 from fit: "<<chi2<<endl;

      if (chi2ndf < bestChi2) {
	bestChi2 = chi2ndf;
	bestCombination = subset;
	bestFitFunction = fitFunction;
      }

      /*std::cout << std::endl << "Ex: " << Ex_values[mappingIndex] << ", fitted functions: " ;
      for (const auto& element : subset) {
        std::cout << element << " ";
      }
      
      std::cout << std::endl << "chi2ndf: " << chi2ndf << std::endl;*/

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
	resultFile << fitFunction->GetParameter(i) << " ";
      }
      resultFile << std::endl;

      resultFile << "Chi2 (from calc): " << chi2ndf << std::endl;
      resultFile << "Chi2 (from TF1): " << chi2 << ", NDF: " << ndf << std::endl;
      resultFile << "------------------------------------------------------------" << std::endl << std::endl;
    }

    nr_of_functions = bestCombination.size();
    
    for(int i=0; i<nr_of_functions; i++){
      functions_ids[i] = bestCombination[i];
    }


    bestFitFunction->SetLineColor(kRed);
    bestFitFunction->SetLineWidth(4);
    bestFitFunction->Draw("SAME");
    cout<<"Best combination: "<<endl;
    for(const auto& element : bestCombination)
      cout<<element<<endl;

    canvas->SetLogy();

    legend->Draw();
    if (checking_gs)
      canvas->SaveAs(Form("plots_17O/ang_dist/checking_gs/fit_5255_%lu.png", mappingIndex + 1));
    else
      canvas->SaveAs(Form("plots_17O/ang_dist/fit_Ex_%lu.png", mappingIndex + 1));
  }
  
  resultFile.close();
}
