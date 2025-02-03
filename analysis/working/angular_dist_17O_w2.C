#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};


double combinedDWBA(double *x, double *par) {
  double result = 0;
    //cout<<"Number of functions: "<<nr_of_functions<<endl;
  for (size_t i = 0; i < nr_of_functions; ++i) {
    double dwba_val = graphsDWBA[functions_ids[i]]->Eval(x[0]);
    //cout<<"dwba_val = "<<dwba_val<<endl;
    result += par[i] * dwba_val;
  //cout<<"Function: "<<functions_ids[i]<<endl;
  }
  //cout<<"result = "<<result<<endl;
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

void angular_dist_17O_w2() {

  double Tmin, Tmax, Dt, ThetaMean, sin_x_dx, integral_corr, int_corr_sin;

  std::vector<std::vector<double>> ex1982 = {
    {7,		18.48,	11.48,	12.74,	2.53165228,	45.5787,	115.3894198},
    {19.97,	25.09,	5.12,	22.53,	1.961815663,	29.2412,	57.36584416},
    {26.16,	30.2,	4.04,	28.18,	1.907862135,	20.9531,	39.97562609},
    {31.12,	34.6,	3.48,	32.86,	1.888206768,	14.06493333,	26.55750232},
    {35.37,	38.49,	3.12,	36.93,	1.874617232,	7.83856,	14.69429965}
  }; // Theta min, theta max, theta mean, delta theta, sin(x)dx * 180/pi, integral corr, int. corr * sin(x)dx * 180/pi

  std::vector<std::vector<double>> ex3552 = {
    {11.14,	19.83,	8.69,	15.485,	2.320109111,	26.5151,	61.5179251},
    {21.25,	26.25,	5,	23.75,	2.013733447,	48.4343,	97.5337699},
    {27.35,	31.37,	4.02,	29.36,	1.970987554,	35.64106667,	70.24809883},
    {32.24,	35.74,	3.5,	33.99,	1.9566687,	17.12853333,	33.51486505}
  };

  std::vector<std::vector<double>> ex3630 = {
    {10.07,	19.51,	9.44,	14.79,	2.409814986,	106.918,	257.6525987},
    {20.96,	26.03,	5.07,	23.495,	2.021252025,	36.4969,	73.76943302},
    {27.14,	31.2,	4.06,	29.17,	1.978854275,	11.09661333,	21.95858074},
    {32.07,	35.59,	3.52,	33.83,	1.959691857,	11.15053333,	21.85160937}
  };

  std::vector<std::vector<double>> ex3920 = {
    {7,		18.24,	11.24,	12.62,	2.455758891,	55.4257,	136.1121555},
    {19.8,	25.17,	5.37,	22.485,	2.053711111,	29.5664,	60.7208442},
    {26.32,	30.51,	4.19,	28.415,	1.993830287,	24.1908,	48.2323497},
    {31.41,	35.02,	3.61,	33.215,	1.977493989,	16.04573333,	31.73034122}
  };

  std::vector<std::vector<double>> ex5255 = {
    {11.93,	20.37,	8.44,	16.15,	2.347611196,	36.6722,	86.09206729},
    {21.86,	26.94,	5.08,	24.4,	2.098570501,	31.14733333,	65.36487493},
    {27.99,	32.09,	4.1,	30.04,	2.052478357,	36.7932,	75.51724669}
  };

  std::vector<std::vector<double>> ex3530 = {
    {11.38,	19.93,	8.55,	15.65,	2.30645002,	132.106,	307.0023357},
    {21.25,	26.25,	4.99,	23.82,	2.01528462,	83.5226,	168.3218109},
    {27.35,	31.37,	4.01,	29.42,	1.972,		46.512,		91.721664},
    {32.24,	35.74,	3.5,	34.04,	1.9552,		26.4529333,	51.72077525}
  };

  std::vector<std::vector<std::vector<double>>> data = {ex1982, ex3552, ex3630, ex3920, ex3530};
  std::vector<double> Ex_values = {1.982, 3.552, 3.630, 3.920, 5.255, 3.530};

  std::vector<TGraph*> graphs;

  TFile *file = TFile::Open("DWBA_17O_more_states.root");
  
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

  std::vector<std::pair<std::vector<std::vector<double>>, std::vector<int>>> fitMappings = {
    {ex1982, {0, 1, 2, 3, 4}},
    {ex3552, {5, 6, 7, 8, 9}},
    {ex3630, {10, 11, 12}},
    {ex3920, {13, 14, 15, 16, 17, 18}},
    {ex5255, {19, 20, 21, 22}},
    {ex3530, {5, 6, 7, 8, 9, 10, 11, 12}},
  };

  std::vector<std::vector<std::vector<int>>> fitPairs =  {
    {{0, 2, 3}, {0, 2}, {0, 3}, {1, 4}},
    {{5}, {6}, {5, 6}, {8, 9}},
    {{10}, {11}, {12}, {10, 11}},
    {{13}, {17, 18}, {15, 16}, {14, 15}},
    {{19}, {20, 21}, {19, 22}},
    {{5}, {10}, {5, 6}, {8, 9}},
  };

  std::vector<int> color = { 1, 629, 596, 418, 801, 905, 8, 9};
  int ncolor = 0;

  for (size_t mappingIndex = 0; mappingIndex < fitMappings.size(); ++mappingIndex) {
    const auto& dataSet = fitMappings[mappingIndex].first;
    const auto& fitIndices = fitMappings[mappingIndex].second;
    auto allCombinations = generateCombinations(fitIndices);

    double bestChi2 = 1e9;
    std::vector<int> bestCombination;
    TF1* bestFitFunction = nullptr;

    ncolor = 0;

    std::vector<double> theta_means;
    std::vector<double> int_corr_sins;

    for (const auto& det : dataSet) {
      Tmin = det[0];
      Tmax = det[1];
      Dt = det[2];
      ThetaMean = det[3];
      sin_x_dx = det[4];
      integral_corr = det[5];
      int_corr_sin = det[6];
      
      theta_means.push_back(det[3]);
      int_corr_sins.push_back(det[5] * det[4] / 10.0);
    }

    TGraph* experimentGraph = new TGraph(theta_means.size(), &theta_means[0], &int_corr_sins[0]);
    experimentGraph->SetTitle("");
    experimentGraph->SetMarkerStyle(20);
    experimentGraph->SetMarkerSize(1);
    experimentGraph->SetMarkerColor(kBlack);
    experimentGraph->SetLineWidth(0);

    experimentGraph->GetXaxis()->SetRangeUser(0, 60);
    experimentGraph->GetYaxis()->SetRangeUser(.01, 50);

    TCanvas* canvas = new TCanvas(Form("fit_canvas_%lu", mappingIndex + 1), 
				  Form("Ex = %.3f MeV", Ex_values[mappingIndex]), 1600, 1200);

    experimentGraph->Draw("AP");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(22);
    latex.DrawLatexNDC(0.5, 0.92, Form("Plot Ex = %.3f MeV", Ex_values[mappingIndex]));

    TLegend* legend = new TLegend(0.12, 0.15, 0.89, 0.35);
    legend->SetTextSize(0.03);
    legend->SetMargin(0.1);
    legend->SetBorderSize(0);
    legend->SetFillColor(0); 
    //legend->AddEntry(experimentGraph, "Data Points", "p");


    /*TF1* fitFunction = new TF1("fitFunction", combinedDWBA, 0, 60, fitIndices.size());
    for (size_t i = 0; i < fitIndices.size(); ++i) {
      fitFunction->SetParameter(i, 1.0);
    }

    experimentGraph->Fit(fitFunction, "R");

    std::vector<double> bestScales;
    for (size_t i = 0; i < fitIndices.size(); ++i) {
      bestScales.push_back(fitFunction->GetParameter(i));
      std::cout << "DWBA " << fitIndices[i] << " Scale Factor: " << bestScales[i] << std::endl;
      }*/

    for (int fitIndex : fitIndices) {
      std::cout << "dataSet: " << mappingIndex << " fitIndex: " << fitIndex << std::endl;
      TGraph* fitGraph = graphsDWBA[fitIndex];
      fitGraph->SetLineColor(color[ncolor]);
      fitGraph->SetLineWidth(2);
      fitGraph->Draw("L");

      ncolor++;

      /*double chi2 = 0.0;
	int nPoints = experimentGraph->GetN();
	for (int i = 0; i < nPoints; ++i) {
	double xData, yData;
	experimentGraph->GetPoint(i, xData, yData);

	double yFit = fitGraph->Eval(xData);

	double sigmaData = 0.1 * yData; 
	if (sigmaData == 0) sigmaData = 1.0;

	chi2 += pow((yData - yFit) / sigmaData, 2);
	}*/

      //std::cout << "Fit Index: " << fitIndex << ", Chi2: " << chi2 << ", NDF: " << nPoints << std::endl;

      //legend->AddEntry(fitGraph, Form("%s, #chi^{2}/ndf = %.3f", fitGraph->GetName(), chi2 / nPoints), "l");
      
      legend->AddEntry(fitGraph, Form("Func %d: %s", fitIndex, fitGraph->GetName()), "l");

    }
    bestChi2 = 1000;
    for (const auto& subset : fitPairs[mappingIndex]) {
      TF1* fitFunction = new TF1("fitFunction", combinedDWBA, 0, 60, subset.size());
      for (size_t i = 0; i < subset.size(); ++i) {
	fitFunction->SetParameter(i, 0.0);
    fitFunction->SetParLimits(i, 0.0, 5);
      }

      nr_of_functions = subset.size();
      for(int i=0; i<nr_of_functions; i++){
        functions_ids[i] = subset[i];
      }

      experimentGraph->Fit(fitFunction, "RL");

      double chi2 = fitFunction->GetChisquare();
      int ndf = fitFunction->GetNDF();

      double x, y;
      double chi2ndf = 0;
      for(int i=0; i< experimentGraph->GetN(); i++){
        experimentGraph->GetPoint(i, x, y);
        cout<<"x: "<<x<<", y: "<<y<<", func: "<<fitFunction->Eval(x)<<", diff: "<<y - fitFunction->Eval(x)<<endl;
        chi2ndf += (y - fitFunction->Eval(x))* (y - fitFunction->Eval(x));
      }
      //double chi2ndf = (ndf > 0) ? (chi2 / ndf) : chi2;
      //double chi2ndf = chi2 ;
      cout<<"Chi2 from calculation: "<<chi2ndf<<", chi2 from fit: "<<chi2<<endl;
    //fitFunction->SetLineColor(kBlue);
    //fitFunction->SetLineWidth(4);
    //fitFunction->Draw("SAME");


      if (chi2ndf < bestChi2) {
	bestChi2 = chi2ndf;
	bestCombination = subset;
	bestFitFunction = fitFunction;
      }

      std::cout << std::endl << "Ex: " << Ex_values[mappingIndex] << ", fitted functions: " ;
      for (const auto& element : subset) {
        std::cout << element << " ";
      }
      
      std::cout << std::endl << "chi2ndf: " << chi2ndf << std::endl;
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
    /*for (const auto& subset : allCombinations) {
      TF1* fitFunction = new TF1("fitFunction", combinedDWBA, 0, 60, subset.size());
      for (size_t i = 0; i < subset.size(); ++i) {
	fitFunction->SetParameter(i, 1.0);
      }

      experimentGraph->Fit(fitFunction, "R");

      double chi2 = fitFunction->GetChisquare();
      int ndf = fitFunction->GetNDF();
      double chi2ndf = (ndf > 0) ? (chi2 / ndf) : chi2;

      if (chi2ndf < bestChi2) {
	bestChi2 = chi2ndf;
	bestCombination = subset;
	bestFitFunction = fitFunction;
      }
    }

    bestFitFunction->SetLineColor(kRed);
    bestFitFunction->SetLineWidth(4);
    bestFitFunction->Draw("SAME");
    */

    canvas->SetLogy();

    legend->Draw();
    canvas->SaveAs(Form("ang_dist_17O/fit_Ex_%lu_w2.png", mappingIndex + 1));
  }
  
}
