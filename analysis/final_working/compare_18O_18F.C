#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <vector>
#include <iostream>

std::vector<TGraph*> graphsDWBA;
int nr_of_functions = 0;
int functions_ids[10];

double combinedDWBA(double *x, double *par)
{
  double result = 0;

  for (int i = 0; i < nr_of_functions; i++)
    {
      double val = graphsDWBA[functions_ids[i]]->Eval(x[0]);
      result += par[i] * val;
    }

  return result;
}

TGraphErrors* buildGraph(const std::vector<std::vector<double>>& data)
{
  std::vector<double> x,y,ex,ey;

  for(auto &d : data)
    {
      double theta = d[3];
      double sin_dx = d[4];
      double val = d[5];
      double unc = d[6];

      x.push_back(theta);
      y.push_back(val/sin_dx);
      ex.push_back(0);
      ey.push_back(unc/sin_dx);
    }

  auto g = new TGraphErrors(x.size(),&x[0],&y[0],&ex[0],&ey[0]);

  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.5);
  g->SetLineWidth(2);

  return g;
}

void drawPanel(
	       TPad* pad,
	       TGraphErrors* data,
	       std::vector<int> dwba_ids,
	       std::vector<TString> labels,
	       TString title,
	       bool drawYaxis,
	       bool drawLegend)
{
  pad->cd();
  pad->SetLogy();

  data->Draw("AP");

  data->GetXaxis()->SetLimits(0,60);
  data->GetHistogram()->GetYaxis()->SetRangeUser(0.01,200);

  if(!drawYaxis)
    data->GetYaxis()->SetLabelSize(0);

  data->GetXaxis()->SetTitle("#theta_{CM} (deg)");
  data->GetYaxis()->SetTitle("d#sigma/d#Omega (arb. units)");

  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.055);
  t.DrawLatex(0.15,0.88,title);

  TLegend *leg = nullptr;

  if(drawLegend)
    {
      leg = new TLegend(0.55,0.65,0.9,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
    }

  nr_of_functions = dwba_ids.size();

  for(int i=0;i<nr_of_functions;i++)
    functions_ids[i]=dwba_ids[i];

  TF1 *fit = new TF1("fit",combinedDWBA,0,60,nr_of_functions);

  for(int i=0;i<nr_of_functions;i++)
    fit->SetParameter(i,1.0);

  data->Fit(fit,"R0");

  fit->SetLineWidth(3);
  fit->SetLineColor(kRed);
  fit->Draw("same");

  if(drawLegend)
    leg->AddEntry(fit,"Total fit","l");

  for(int i=0;i<nr_of_functions;i++)
    {
      TF1 *comp = new TF1(Form("c%i",i),combinedDWBA,0,60,nr_of_functions);

      for(int j=0;j<nr_of_functions;j++)
	comp->SetParameter(j,(j==i)?fit->GetParameter(j):0);

      comp->SetLineStyle(2);
      comp->SetLineWidth(3);
      comp->SetLineColor(600+i*50);

      comp->Draw("same");

      if(drawLegend)
	leg->AddEntry(comp,labels[i],"l");
    }

  data->Draw("P same");

  if(drawLegend)
    leg->Draw();
}

void compare_18O_18F()
{
  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c","18O vs 18F",1600,2000);

  double w = 0.5;
  double h = 1.0/3.0;

  TPad *pads[6];

  for(int r=0;r<3;r++)
    {
      for(int col=0;col<2;col++)
        {
	  int i = r*2 + col;

	  double x1 = col*w;
	  double x2 = (col+1)*w;

	  double y2 = 1.0 - r*h;
	  double y1 = y2 - h;

	  pads[i] = new TPad(Form("pad%i",i),"",x1,y1,x2,y2);
	  pads[i]->SetMargin(
			     col==0 ? 0.18 : 0.02,
			     0.02,
			     r==2 ? 0.18 : 0.02,
			     0.02);

	  pads[i]->Draw();
        }
    }

  TFile *fileO = TFile::Open("DWBA_17O_AK.root");
  TFile *fileF = TFile::Open("DWBA_17F_AK.root");

  auto arrO = (TObjArray*)fileO->Get("qList");
  auto arrF = (TObjArray*)fileF->Get("qList");

  for(int i=0;i<arrO->GetEntries();i++)
    graphsDWBA.push_back((TGraph*)arrO->At(i));

  for(int i=0;i<arrF->GetEntries();i++)
    graphsDWBA.push_back((TGraph*)arrF->At(i));

  /* =======================
     DATA SELECTION
     ======================= */

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

  auto g_O_1982 = buildGraph(ex1982_half_dets);
  auto g_O_3552 = buildGraph(ex3552_3630_half_dets);
  auto g_O_3920 = buildGraph(ex3920_half_dets);

  auto g_F_3062 = buildGraph(ex3061_half_dets);
  auto g_F_4652 = buildGraph(ex4652_4753_half_dets);
  auto g_F_4964 = buildGraph(ex4964_half_dets);

  drawPanel(pads[0],g_O_1982,{0,2},
	    {"1s_{1/2}","0d_{5/2}"},
	    "^{18}O 1.982 MeV",
	    true,true);

  drawPanel(pads[1],g_F_3062,{5,6},
	    {"1s_{1/2}","0d_{5/2}"},
	    "^{18}F 3.062 MeV",
	    false,false);

  drawPanel(pads[2],g_O_3552,{36},
	    {"0d_{5/2}"},
	    "^{18}O 3.552 & 3.630 MeV",
	    true,false);

  drawPanel(pads[3],g_F_4652,{29},
	    {"0d_{5/2}"},
	    "^{18}F 4.652 & 4.753 MeV",
	    false,false);

  drawPanel(pads[4],g_O_3920,{15,13},
	    {"1s_{1/2}","0d_{5/2}"},
	    "^{18}O 3.920 MeV",
	    true,false);

  drawPanel(pads[5],g_F_4964,{2,3},
	    {"1s_{1/2}","0d_{5/2}"},
	    "^{18}F 4.964 MeV",
	    false,false);

  c->SaveAs("18O_18F_comparison.pdf");
}
