#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

// Global vectors for DWBA graphs
std::vector<TGraph*> graphsDWBA_O;
std::vector<TGraph*> graphsDWBA_F;
int nr_of_functions = 0;
int functions_ids[10] = {0,0,0,0,0,0,0,0,0,0};
bool use_O_graphs = true;

// Combined DWBA function
double combinedDWBA(double *x, double *par) {
  double result = 0;
  std::vector<TGraph*>* graphsToUse = use_O_graphs ? &graphsDWBA_O : &graphsDWBA_F;
  for (size_t i = 0; i < nr_of_functions; ++i) {
    double dwba_val = (*graphsToUse)[functions_ids[i]]->Eval(x[0]);
    result += par[i] * dwba_val;
  }
  return result;
}

void combined_18O_18F_plot() {
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  // Data for 18O states
  std::vector<std::vector<double>> ex1982_half_dets = {
    {8.00 , 14.74, 6.74, 11.37, 1.3287, 28.4242, 1.42521},
    {14.74, 18.48, 3.73, 16.61, 1.0674, 17.2772, 1.11639},
    {19.97, 22.69, 2.71, 21.33, 0.9874, 16.2642, 1.07688},
    {22.69, 25.09, 2.40, 23.89, 0.9722, 12.7089, 0.99956},
    {26.16, 28.25, 2.09, 27.20, 0.9562, 10.2986, 0.86965},
    {28.25, 30.20, 1.95, 29.22, 0.9511, 10.5578, 0.88115},
    {31.12, 32.90, 1.78, 32.01, 0.9453, 6.15872 * 4./3., 0.559427 * 4./3.},
    {32.90, 34.60, 1.69, 33.75, 0.9413, 4.24311 * 4./3., 0.590164 * 4./3.},
    {35.37, 36.96, 1.59, 36.16, 0.9410, 3.42685 * 4./3., 0.524537 * 4./3.},
    {36.96, 38.49, 1.53, 37.73, 0.9365, 2.27111 * 4./3., 0.419665 * 4./3.}
  };

  std::vector<std::vector<double>> ex3552_3630_half_dets = {
    {10.56, 16.16, 5.60, 13.36, 1.2937, 63.0214 + 7.39285, 3.9},
    {16.16, 19.65, 3.49, 17.91, 1.0734, 66.6946, 2.21898},
    {21.09, 23.75, 2.67, 22.42, 1.0168, 48.8069, 1.86560},
    {23.75, 26.13, 2.38, 24.94, 1.0034, 37.5365, 1.63665},
    {27.23, 29.32, 2.09, 28.28, 0.9900, 20.0825 * 4./3., 1.0096 * 4./3.},
    {29.32, 31.27, 1.95, 30.30, 0.9855, 15.6256 * 4./3., 1.0757 * 4./3.},
    {32.15, 33.95, 1.80, 33.05, 0.9809, 11.3624 * 4./3., 0.9295 * 4./3.},
    {33.95, 35.66, 1.71, 34.80, 0.9765, 10.1256 * 4./3., 0.8921 * 4./3.}
  };

  std::vector<std::vector<double>> ex3920_half_dets = {
    { 8.00, 14.20, 6.20, 11.10, 1.1936, 26.7301, 1.43482},
    {14.20, 18.24, 4.03, 16.22, 1.1268, 22.9592, 1.34987},
    {19.80, 22.66, 2.86, 21.23, 1.0352, 14.5072, 1.03633},
    {22.66, 25.17, 2.51, 23.91, 1.0168, 15.1933, 1.05784},
    {26.32, 28.49, 2.17, 27.41, 1.0007, 8.45285 * 4./3., 0.67275 * 4./3.},
    {28.49, 30.51, 2.02, 29.50, 0.9951, 9.04236 * 4./3., 0.8328 * 4./3.},
    {31.41, 33.26, 1.85, 32.34, 0.9897, 7.50488 * 4./3., 0.7676 * 4./3.},
    {33.26, 35.02, 1.76, 34.14, 0.9879, 6.11097 * 4./3., 0.711327 * 4./3.}
  };

  // Data for 18F states
  std::vector<std::vector<double>> ex3061_half_dets = {
    {8.6  , 14.48, 5.88, 11.54, 1.1763,  5.83176, 0.658813},
    {14.48, 18.32, 3.84, 16.40, 1.0839,  6.27474, 0.702328},
    {19.84, 22.59, 2.75, 21.22, 0.9967,  4.93428, 0.60},
    {22.59, 25.02, 2.43, 23.81, 0.9804,  4.64378, 0.608},
    {26.10, 28.21, 2.11, 27.15, 0.9629,  3.41079, 0.51562},
    {28.21, 30.17, 1.96, 29.19, 0.9578,  3.23499, 0.497581},
    {31.10, 32.90, 1.80, 32.00, 0.9517,  1.07583 * 4./3., 0.309092 * 4./3.},
    {32.90, 34.61, 1.71, 33.75, 0.9478,  1.91486 * 4./3., 0.382529 * 4./3.},
    {35.38, 36.98, 1.61, 36.18, 0.9476,  0.76595 * 4./3., 0.272617 * 4./3.},
    {36.98, 38.52, 1.54, 37.75, 0.9427,  1.63772 * 4./3., 0.365782 * 4./3.}
  };

  std::vector<std::vector<double>> ex4652_4753_half_dets = {
    { 8.20, 15.72, 7.52, 11.96, 1.5583, (21.0052 + 0.436834), 1.6249},
    {15.72, 19.35, 3.63, 17.54, 1.0922, 20.5967, 1.27427},
    {20.82, 23.54, 2.72, 22.18, 1.0274, 15.1561, 1.05519},
    {23.54, 25.96, 2.42, 24.75, 1.0139, 13.0476, 1.00247},
    {27.08, 29.19, 2.12, 28.13, 0.9982, 6.08755 * 4./3., 0.673348 * 4./3.},
    {29.19, 31.17, 1.98, 30.18, 0.9937, 4.54548 * 4./3., 0.596529 * 4./3.},
    {32.05, 33.87, 1.82, 32.96, 0.9879, (3.86017) * 4./3., 0.562701 * 4./3.},
    {33.87, 35.60, 1.73, 34.73, 0.9849, (3.65818 + 0.517828) * 4./3., (sqrt( pow(1.27, 2) + pow(1.17135, 2) )) * 4./3.}
  };

  std::vector<std::vector<double>> ex4964_half_dets = {
    {8.2  , 14.07, 5.87, 11.13, 1.1336, 8.57289, 0.847075},
    {14.07, 18.17, 4.10, 16.12, 1.1386, 7.30343, 0.783341},
    {19.75, 22.64, 2.89, 21.20, 1.0432, 4.63121, 0.603035},
    {22.64, 25.17, 2.53, 23.90, 1.0240, 4.8697,  0.637558},
    {26.32, 28.51, 2.19, 27.42, 1.0075, 3.66391 * 4./3., 0.528969 * 4./3.},
    {28.51, 30.54, 2.03, 29.53, 1.0017, 2.2864 * 4./3.,  0.430266 * 4./3.},
    {31.45, 33.31, 1.86, 32.38, 0.9961, 3.12153 * 4./3., 0.515049 * 4./3.},
    {33.31, 35.08, 1.77, 34.19, 0.9942, 1.73812 * 4./3., 0.51603 * 4./3.}
  };

  TFile *fileO = TFile::Open("../corrected_working/DWBA_17O_AK.root");
  if (!fileO || fileO->IsZombie()) {
    std::cerr << "Error opening 18O ROOT file!" << std::endl;
    return;
  }

  TObjArray* objArrayO = (TObjArray*)fileO->Get("qList");
  if (objArrayO) {
    for (int i = 0; i < objArrayO->GetEntries(); ++i) {
      TObject* obj = objArrayO->At(i);
      if (obj && obj->InheritsFrom(TGraph::Class())) {
        TGraph* graph = (TGraph*)obj;
        graphsDWBA_O.push_back(graph);
      }
    }
  }

  TFile *fileF = TFile::Open("../corrected_working/DWBA_17F_AK.root");
  if (!fileF || fileF->IsZombie()) {
    std::cerr << "Error opening 18F ROOT file!" << std::endl;
    return;
  }

  TObjArray* objArrayF = (TObjArray*)fileF->Get("qList");
  if (objArrayF) {
    for (int i = 0; i < objArrayF->GetEntries(); ++i) {
      TObject* obj = objArrayF->At(i);
      if (obj && obj->InheritsFrom(TGraph::Class())) {
        TGraph* graph = (TGraph*)obj;
        graphsDWBA_F.push_back(graph);
      }
    }
  }

  TCanvas* canvas = new TCanvas("combined_canvas", "18O and 18F States", 2400, 3600);
  canvas->Divide(2, 3, 0.0, 0.0);

  std::vector<std::vector<std::vector<double>>> dataSets = {
    ex1982_half_dets, ex3061_half_dets,
    ex3552_3630_half_dets, ex4652_4753_half_dets,
    ex3920_half_dets, ex4964_half_dets
  };

  std::vector<std::string> titles = {
    "^{18}O: 1.982 MeV", "^{18}F: 3.062 MeV",
    "^{18}O: 3.552 & 3.630 MeV", "^{18}F: 4.652 & 4.753 MeV",
    "^{18}O: 3.920 MeV", "^{18}F: 4.964 MeV"
  };

  std::vector<std::vector<int>> fitIndices = {
    {0, 2},    // 18O 1.982 MeV
    {5, 6},    // 18F 3.061 MeV
    {36},      // 18O 3.552 & 3.630 MeV
    {29},      // 18F 4.652 & 4.753 MeV
    {15, 13},  // 18O 3.920 MeV
    {2, 3}     // 18F 4.964 MeV
  };

  std::vector<std::vector<TString>> labels = {
    {"1 s_{1/2} 2^{+}", "0 d_{5/2} 2^{+}"},        // 18O 1.982
    {"1 s_{1/2} 2^{+}", "0 d_{5/2} 2^{+}"},        // 18F 3.061
    {"0 d_{5/2} 4^{+}"},                           // 18O 3.552
    {"0 d_{5/2} 4^{+}"},                           // 18F 4.652
    {"1 s_{1/2} 2^{+}", "0 d_{5/2} 2^{+}"},        // 18O 3.920
    {"1 s_{1/2} 2^{+}", "0 d_{5/2} 2^{+}"}         // 18F 4.964
  };

  std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kCyan, kOrange};

  for (int iPad = 0; iPad < 6; ++iPad) {
    canvas->cd(iPad + 1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);

    bool isOxygen = (iPad % 2 == 0);
    use_O_graphs = isOxygen;

    if (isOxygen){
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.0);
    }
    else{
      gPad->SetLeftMargin(0.0);
      gPad->SetRightMargin(0.15);
    }

    const auto& dataSet = dataSets[iPad];
    std::vector<double> theta_means;
    std::vector<double> int_corr_sins;
    std::vector<double> int_unc;
    std::vector<double> x_unc;

    for (const auto& det : dataSet) {
      theta_means.push_back(det[3]);
      int_corr_sins.push_back(det[5] / det[4]);
      int_unc.push_back(det[6] / det[4]);
      x_unc.push_back(0.0);
    }

    TGraphErrors* experimentGraph = new TGraphErrors(theta_means.size(), 
						     &theta_means[0], 
						     &int_corr_sins[0], 
						     &x_unc[0], 
						     &int_unc[0]);
    
    experimentGraph->SetTitle("");
    experimentGraph->SetMarkerStyle(20);
    experimentGraph->SetMarkerSize(1.2);
    experimentGraph->SetMarkerColor(kBlack);
    experimentGraph->SetLineColor(kBlack);
    experimentGraph->SetLineWidth(2);

    experimentGraph->Draw("APE");
    
    experimentGraph->GetXaxis()->SetRangeUser(0, 45);
    experimentGraph->GetYaxis()->SetRangeUser(0.1, 200);
    
    if (iPad >= 4) {
      experimentGraph->GetXaxis()->SetTitle("#theta_{CM} (deg)");
      experimentGraph->GetXaxis()->SetTitleSize(0.05);
      experimentGraph->GetXaxis()->SetLabelSize(0.04);
    } else {
      experimentGraph->GetXaxis()->SetLabelSize(0);
    }
    
    if (iPad % 2 == 0) {
      experimentGraph->GetYaxis()->SetTitle("d#sigma/d#Omega (arb. units)");
      experimentGraph->GetYaxis()->SetTitleSize(0.05);
      experimentGraph->GetYaxis()->SetLabelSize(0.04);
    } else {
      experimentGraph->GetYaxis()->SetLabelSize(0);
    }

    const auto& indices = fitIndices[iPad];
    TF1* fitFunction = new TF1(Form("fit_%d", iPad), combinedDWBA, 0, 45, indices.size());
    
    for (size_t i = 0; i < indices.size(); ++i) {
      fitFunction->SetParameter(i, 1.0);
      fitFunction->SetParLimits(i, 0.0, 100.0);
    }

    nr_of_functions = indices.size();
    for(int i = 0; i < nr_of_functions; i++){
      functions_ids[i] = indices[i];
    }

    experimentGraph->Fit(fitFunction, "RQ0");
    
    fitFunction->SetLineColor(kRed);
    fitFunction->SetLineWidth(3);
    fitFunction->Draw("SAME");

    if (indices.size() > 1) {
      for (size_t comp = 0; comp < indices.size(); ++comp) {
        TF1* compFunc = new TF1(Form("comp_%d_%lu", iPad, comp), combinedDWBA, 0, 45, indices.size());
        
        for (size_t i = 0; i < indices.size(); ++i) {
          if (i == comp) compFunc->SetParameter(i, fitFunction->GetParameter(i));
          else compFunc->SetParameter(i, 0.0);
        }
        
        compFunc->SetLineColor(colors[comp]);
        compFunc->SetLineStyle(2);
        compFunc->SetLineWidth(2);
        compFunc->Draw("SAME");
      }
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.055);
    latex.SetTextFont(42);
    latex.DrawLatex(0.20, 0.85, titles[iPad].c_str());

    if (indices.size() > 1 && labels[iPad].size() > 0) {
      TLegend* legend = new TLegend(0.50, 0.65, 0.90, 0.82);
      legend->SetTextSize(0.04);
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      
      legend->AddEntry(fitFunction, "Total fit", "l");
      for (size_t i = 0; i < labels[iPad].size(); ++i) {
        TF1* dummy = new TF1("", "1", 0, 1);
        dummy->SetLineColor(colors[i]);
        dummy->SetLineStyle(2);
        dummy->SetLineWidth(2);
        legend->AddEntry(dummy, labels[iPad][i], "l");
      }
      legend->Draw();
    }
  }

  canvas->cd();
  TLatex mainTitle;
  mainTitle.SetNDC();
  mainTitle.SetTextSize(0.025);
  mainTitle.SetTextFont(42);
  mainTitle.DrawLatex(0.5, 0.98, "DWBA Fits for Selected States in ^{18}O and ^{18}F");

  canvas->SaveAs("combined_18O_18F_states.pdf");
  canvas->SaveAs("combined_18O_18F_states.png");
}
