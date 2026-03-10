#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TKey.h>
#include <TCollection.h>
#include <TObject.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

void export_canvas_data(const char* filename, int npoints = 1000) {
  TFile *f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file " << filename << std::endl;
    return;
  }

  TCanvas *canvas = nullptr;
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom("TCanvas")) {
      canvas = (TCanvas*)obj;
      break;
    }
  }

  if (!canvas) {
    std::cerr << "No TCanvas found in file" << std::endl;
    return;
  }

  std::vector<TGraphErrors*> graphs;
  std::vector<TF1*> functions;

  TCollection *prims = canvas->GetListOfPrimitives();
  TIter primsIter(prims);
  TObject *prim;
  while ((prim = primsIter())) {
    if (prim->InheritsFrom("TGraphErrors")) {
      graphs.push_back((TGraphErrors*)prim);
    } else if (prim->InheritsFrom("TF1")) {
      functions.push_back((TF1*)prim);
    }
  }

  if (graphs.empty() && functions.empty()) {
    std::cerr << "No TGraphErrors or TF1 found on canvas" << std::endl;
    return;
  }

  double xmin = 1e30, xmax = -1e30;
  for (auto f1 : functions) {
    xmin = std::min(xmin, f1->GetXmin());
    xmax = std::max(xmax, f1->GetXmax());
  }

  std::string funcfile = std::string(filename);
  size_t pos = funcfile.find(".root");
  if (pos != std::string::npos) funcfile.replace(pos, 5, "_functions.txt");
  std::ofstream fout_func(funcfile);
  if (!fout_func.is_open()) {
    std::cerr << "Cannot open function output file " << funcfile << std::endl;
    return;
  }

  fout_func << "# x ";
  for (size_t i = 0; i < functions.size(); ++i) fout_func << "f" << i << " ";
  fout_func << "\n";

  double dx = (xmax - xmin) / (npoints - 1);
  for (int i = 0; i < npoints; ++i) {
    double x = xmin + i * dx;
    fout_func << x << " ";
    for (auto f1 : functions) fout_func << f1->Eval(x) << " ";
    fout_func << "\n";
  }
  fout_func.close();
  std::cout << "Functions exported to " << funcfile << " with " << npoints << " points." << std::endl;

  std::string graphfile = std::string(filename);
  if (pos != std::string::npos) graphfile.replace(pos, 5, "_graphs.txt");
  std::ofstream fout_graph(graphfile);
  if (!fout_graph.is_open()) {
    std::cerr << "Cannot open graph output file " << graphfile << std::endl;
    return;
  }

  fout_graph << "# ";
  for (size_t i = 0; i < graphs.size(); ++i) fout_graph << "graph" << i << "_x graph" << i << "_y graph" << i << "_err ";
  fout_graph << "\n";

  size_t max_points = 0;
  for (auto g : graphs) max_points = std::max(max_points, (size_t)g->GetN());

  for (size_t i = 0; i < max_points; ++i) {
    for (auto g : graphs) {
      if ((int)i < g->GetN()) {
	double x, y;
	g->GetPoint(i, x, y);
	double err = g->GetErrorY(i);
	fout_graph << x << " " << y << " " << err << " ";
      } else {
	fout_graph << "nan nan nan ";
      }
    }
    fout_graph << "\n";
  }
  fout_graph.close();
  std::cout << "Graphs exported to " << graphfile << std::endl;

  f->Close();
}
