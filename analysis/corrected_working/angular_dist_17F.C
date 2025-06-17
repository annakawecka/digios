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

void angular_dist_17F() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, integral_unc;

  std::vector<std::vector<double>> ex4963_positive_sigma_mean_fixedDist = {
    {8.0,	18.17,	10.17,	13.09,	2.3033, 13.5836 * 4./4., 0.},
    {19.75,	25.17,	5.41,	22.46,	2.0678, 7.72417 * 4./4., 0.},
    {26.32,	30.54,	4.22,	28.43,	2.0095, 5.00161 * 4./3., 0.},
    {31.45,	35.08,	3.63,	33.27,	1.9905, 4.83654 * 4./3., 0.}
  };
  // Theta min, theta max, delta theta, theta mean, sin(x)dx * 180/pi, integral corr, integral corr error
  // 0		1	   2		3	    4		       5

  std::vector<std::vector<double>> ex1041_positive_sigma_mean_fixedDist = {
    {  8.0,	18.58,	10.58,	13.29,	2.4321, 9.03912 * 4./4., 0.},
    {20.02,	24.96,	4.94,	22.49,	1.8900, 2.86952 * 4./4., 0.},
    {26.01,	29.93,	3.92,	27.97,	1.8379, 0.0673394 * 4./4., 0.},
    {30.79,	34.18,	3.39,	32.48,	1.8189, 2.82173e-10 * 4./4., 0.},
    {34.97,	38.01,	3.04,	36.49,	1.8085, 6.55208e-08 * 4./3., 0.},
    {38.69,	41.49,	2.80,	40.09,	1.8022, 4.66247e-13 * 4./3., 0.}
  };

  std::vector<std::vector<double>> ex1041_positive_sigma_low_fixedDist_likelihood = {
    {  8.0,	18.58,	10.58,	13.29,	2.4321, 25.0706 * 4./4., 0.},
    {20.02,	24.96,	4.94,	22.49,	1.8900, 17.9001 * 4./4., 0.},
    {26.01,	29.93,	3.92,	27.97,	1.8379, 0.00430497 * 4./4., 0.},
    {30.79,	34.18,	3.39,	32.48,	1.8189, 0.219477 * 4./4., 0.},
    {34.97,	38.01,	3.04,	36.49,	1.8085, 1.82631e-08 * 4./3., 0.},
    {38.69,	41.49,	2.80,	40.09,	1.8022, 7.92755e-08 * 4./3., 0.}
  };

  std::vector<std::vector<double>> ex3061_positive_sigma_mean_fixedDist_likelihood = {
    {  8.0,	18.32,	10.32,	13.29,	2.3496, 12.0029 * 4./4., 0.},
    {19.84,	25.02,	5.18, 	22.49,	1.9778, 10.0603 * 4./4., 0.},
    {26.10,	30.17,	4.07, 	27.97,	1.9210, 6.69614 * 4./4., 0.},
    {31.10,	34.60,	3.50, 	32.48,	1.8998, 2.70479 * 4./3., 0.},
    {35.38,	38.52,	3.15, 	36.49,	1.8905, 1.95997 * 4./3., 0.}
  };

  std::vector<std::vector<double>> ex4652_positive_sigma_mean_fixedDist_likelihood = {
    {10.05,	19.56,	9.51, 	14.81,	2.4293, 42.8655 * 4./4., 0.},
    {21.01,	26.11,	5.09, 	23.56,	2.0353, 23.033 * 4./4., 0.},
    {27.21,	31.29,	4.07, 	29.25,	1.9894, 8.79789 * 4./3., 0.},
    {32.17,	35.70,	3.53, 	33.93,	1.9711, 5.33186 * 4./3., 0.}
  };

  std::vector<std::vector<double>> ex4753_positive_sigma_mean_fixedDist_likelihood = {
    //{  8.0,	19.13,	11.13,	13.57,	2.6115, 1.11739e-06 * 4./4., 0.},
    {20.61,	25.81,	5.19,  	23.21,	2.0466, 6.80928 * 4./4., 0.},
    {26.93,	31.05,	4.12,  	28.99,	1.9961, 1.69489 * 4./3., 0.},
    {31.94,	35.50,	3.56,  	33.72,	1.9758, 1.87777 * 4./3., 0.}
  };

  std::vector<std::vector<std::vector<double>>> data;
  std::vector<double> Ex_values;

  if (checking_gs) {
    //data = {ex5255, ex5255, ex5255, ex5255, ex5255, ex5255, ex5255};
    //Ex_values = {5.255, 5.255, 5.255, 5.255, 5.255, 5.255, 5.255};
  } else {
    data = {ex4963_positive_sigma_mean_fixedDist, ex1041_positive_sigma_low_fixedDist_likelihood, ex3061_positive_sigma_mean_fixedDist_likelihood,
	    ex4652_positive_sigma_mean_fixedDist_likelihood, ex4753_positive_sigma_mean_fixedDist_likelihood};
    Ex_values = {4.963, 1.041, 3.061, 4.652, 4.753};
  }

  std::vector<TGraph*> graphs;

  //TFile *file = TFile::Open("../working/DWBA_17O_more_states.root");

  TFile *file;

  const auto& pot = potentials[24];

  TString filename = Form("DWBA_17F_%s.root", pot.Data());
  TString outputDir = Form("plots_17F/ang_dist/%s", pot.Data());
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

  if (checking_gs) {
    
  } else {
    fitMappings = {
      {ex4963_positive_sigma_mean_fixedDist, {4, 5, 6, 7, 8}},
      {ex1041_positive_sigma_mean_fixedDist, {0, 1, 2, 3}},
      {ex3061_positive_sigma_mean_fixedDist_likelihood, {9, 10, 11, 12}},
      {ex4652_positive_sigma_mean_fixedDist_likelihood, {13, 14, 15, 16}},
      {ex4753_positive_sigma_mean_fixedDist_likelihood, {17, 18, 19, 20}}
    };
    fitPairs =  {
      {{4}, {5}, {6}, {7}, {8}, {4, 5}, {4, 6}, {4, 8}, {5, 8}, {7, 8}},
      {{0}, {1}, {2}, {3}, {0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}},
      {{9}, {10}, {11}, {12}, {9, 10}, {9, 11}, {9, 12}, {10, 12}},
      {{13}, {14}, {15}, {16}, {13, 14}, {13, 15}, {13, 16}, {14, 16}},
      {{17}, {18}, {19}, {20}, {17, 18}, {17, 19}, {17, 20}, {19, 20}},
    };
  }

  std::cout << "After Mappings " << std::endl;

  std::vector<int> color = { 1, 629, 596, 418, 801, 905, 8, 9};
  int ncolor = 0;

  std::ofstream resultFile;

  if (checking_gs) {
    //resultFile.open("plots_17O/ang_dist/checking_5255_DWBA_fit_results_17O_newUnc.txt");
    resultFile.open("plots_17F/ang_dist/checking_gs_DWBA_fit_results_17O_newUnc.txt");
  } else {
    TString outputFilename = Form("%s/DWBA_fit_results_17F.txt", outputDir.Data());
    resultFile.open(outputFilename);
  }

  for (size_t mappingIndex = 0; mappingIndex < fitMappings.size(); ++mappingIndex) {
    const auto& dataSet = fitMappings[mappingIndex].first;
    const auto& fitIndices = fitMappings[mappingIndex].second;
    //auto allCombinations = generateCombinations(fitIndices);

    double bestChi2 = 1e9;
    std::vector<int> bestCombination;
    TF1* bestFitFunction = nullptr;

    ncolor = 0;

    std::cout << "===========================================================================" << std::endl;
    std::cout << "===========================================================================" << std::endl;
    std::cout << "===========================================================================" << std::endl;
    std::cout << "Energy " << Ex_values[mappingIndex] << std::endl << std::endl;

    /*resultFile << "===========================================================================" << endl;
    resultFile << "===========================================================================" << endl;
    resultFile << "===========================================================================" << endl;
    resultFile << "Energy " << Ex_values[mappingIndex] << endl << endl;*/

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
      //int_unc.push_back(integral_unc / sin_x_dx);
      //int_unc.push_back(0.1);
      x_unc.push_back(0.0);

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
    latex.DrawLatexNDC(0.5, 0.92, Form("Ex = %.3f MeV", Ex_values[mappingIndex]));

    TLegend* legend = new TLegend(0.12, 0.15, 0.72, 0.35);
    legend->SetTextSize(0.03);
    legend->SetMargin(0.1);
    legend->SetBorderSize(0);
    legend->SetFillColor(0); 
    //legend->AddEntry(experimentGraph, "Data Points", "p");

    for (int fitIndex : fitIndices) {
      std::cout << "dataSet: " << mappingIndex << " fitIndex: " << fitIndex  << "    " << graphsDWBA[fitIndex]->GetName() << std::endl;
      //resultFile << "dataSet: " << mappingIndex << " fitIndex: " << fitIndex  << "    " << graphsDWBA[fitIndex]->GetName() << endl;
      TGraph* fitGraph = graphsDWBA[fitIndex];
      fitGraph->SetLineColor(color[ncolor]);
      fitGraph->SetLineWidth(2);
      fitGraph->Draw("L");

      ncolor++;

      legend->AddEntry(fitGraph, Form("Func %d: %s", fitIndex, fitGraph->GetName()), "l");

    }
    bestChi2 = 1000;
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

      double chi2 = 1000000.;

      //experimentGraph->Fit(fitFunction, "RL0");
      TFitResultPtr fitResult = experimentGraph->Fit(fitFunction, "RLS0");
      if (fitResult->IsValid()) {
	chi2 = fitResult->Chi2();

	if (chi2 < bestChi2) {
	  bestChi2 = chi2;
	  bestCombination = subset;
	  bestFitFunction = (TF1*)fitFunction->Clone();
	}
      }
      else {
	std::cout << "Invalid fit " << std::endl;
      }

      //double chi2 = fitFunction->GetChisquare();
      //int ndf = fitFunction->GetNDF();

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

      if (chi2 < bestChi2) {
	bestChi2 = chi2;
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
      resultFile << "Chi2 (from TF1): " << chi2 << std::endl;
      resultFile << "------------------------------------------------------------" << std::endl << std::endl;
    }
    if (bestFitFunction) {
      bestFitFunction->SetLineColor(kRed);
      bestFitFunction->SetLineWidth(4);
      bestFitFunction->Draw("LSAME");
      legend->AddEntry(bestFitFunction, "Best Fit", "l");
    }

    nr_of_functions = bestCombination.size();
    
    for(int i=0; i<nr_of_functions; i++){
      functions_ids[i] = bestCombination[i];
    }


    /*bestFitFunction->SetLineColor(kRed);
    bestFitFunction->SetLineWidth(4);
    bestFitFunction->Draw("SAME");*/
    cout << "Best combination: " << endl;
    resultFile << "Best combination: " << endl;
    for(const auto& element : bestCombination){
      cout << element << "   " << graphsDWBA[element]->GetName() << endl;
      resultFile << element << "   " << graphsDWBA[element]->GetName() << endl;
    }

    canvas->SetLogy();
    canvas->Update();

    legend->Draw();
    if (checking_gs)
      canvas->SaveAs(Form("plots_17F/ang_dist/checking_gs/fit_gs_%lu_newUnc.png", mappingIndex + 1));
      //canvas->SaveAs(Form("plots_17O/ang_dist/checking_gs/fit_5255_%lu_newUnc.png", mappingIndex + 1));
    else {
      TString outputFilename = Form("%s/fit_Ex_%lu.png", outputDir.Data(), mappingIndex + 1);
      canvas->SaveAs(outputFilename);
    }
  }
  
  resultFile.close();
}
