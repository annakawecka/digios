#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

struct FitResult {
  string functionName;
  double chi2;
};

bool compareChi2(const FitResult& a, const FitResult& b) {
  return a.chi2 < b.chi2;
}

void extractAndSortFitData(const char* filename) {
  ifstream file(filename);
  if (!file.is_open()) {
    cerr << "Error opening file!" << endl;
    return;
  }

  string line;
  vector<FitResult> fitResults;
  FitResult currentResult;
  bool expectFunctionName = false;

  while (getline(file, line)) {
    if (line.find("Function Names:") != string::npos) {
      expectFunctionName = true;
      continue;
    }

    if (expectFunctionName) {
      if (line.length() >= 2) {
	currentResult.functionName = line.substr(line.length() - 2);
      } else {
	currentResult.functionName = line;
      }
      expectFunctionName = false;
    }

    if (line.find("Chi2 (from calc):") != string::npos) {
      stringstream ss(line);
      string word;
      ss >> word >> word;
      ss >> word;
      ss >> currentResult.chi2;

      fitResults.push_back(currentResult);
    }
  }

  file.close();

  sort(fitResults.begin(), fitResults.end(), compareChi2);

  cout << "Fit Results Sorted by Decreasing Chi2:" << endl;
  for (const auto& result : fitResults) {
    cout << "Function Name: " << result.functionName
	 << ", Chi2 (from calc): " << result.chi2 << endl;
  }

  cout << "\n% LaTeX table of Fit Results sorted by Chi2\n";
  cout << "\\begin{table}[htbp]\n";
  cout << "\\centering\n";
  cout << "\\begin{tabular}{|c|c|}\n";
  cout << "\\hline\n";
  cout << "Function Name & Chi2 \\\\ \n";
  cout << "\\hline\n";
  for (const auto& result : fitResults) {
    cout << result.functionName << " & " << result.chi2 << " \\\\ \n";
  }
  cout << "\\hline\n";
  cout << "\\end{tabular}\n";
  cout << "\\caption{Fit results sorted by decreasing $\\chi^2$}\n";
  cout << "\\label{tab:fit_results}\n";
  cout << "\\end{table}\n";
}

void sort_fit_data() {
  extractAndSortFitData("plots_17O/ang_dist/checking_5255_DWBA_fit_results_17O_newUnc.txt");
}
