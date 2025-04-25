#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <cmath>

void ExtractPeakInfo_17F() {
  std::vector<double> true_peaks = {
    0.0, 0.9372, 1.04155, 1.08054, 1.12136, 1.4, 1.6, 1.9, 1.70081,
    2.10061, 2.52335, 3.06184, 3.13387, 3.3582, 3.72419, 3.79149,
    3.83917, 4.1159, 4.2258, 4.36015, 4.3981, 4.652, 4.753,
    4.8483, 4.860, 4.9636, 5.2976
  };

  std::ifstream infile("plots_17F/fitting/fit_parameters_17F.txt");  // replace with your actual file path
  if (!infile.is_open()) {
    std::cerr << "Failed to open file\n";
    return;
  }

  std::ofstream outfile("plots_17F/fitting/peak_summary.txt");
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not create output file\n";
    return;
  }

  std::string line;
  std::string current_hist;
  std::vector<std::tuple<std::string, int, double, double, double, double, int>> peaks; // hist, peak#, pos, pos_err, int, int_err, matched_true_peak_index

  int peakNum = 0;
  double pos = 0, pos_err = 0, integral = 0, integral_err = 0;

  while (std::getline(infile, line)) {
    if (line.find("Histogram:") != std::string::npos) {
      current_hist = line.substr(line.find(":") + 2);
      peakNum = 0;
    }

    if (line.find("Peak") != std::string::npos && line.find("Position:") != std::string::npos) {
      size_t pos_colon = line.find(":");
      size_t paren_open = line.find("(");
      size_t paren_close = line.find(")");

      pos = std::stod(line.substr(pos_colon + 1, paren_open - pos_colon - 1));
      pos_err = std::stod(line.substr(paren_open + 1, paren_close - paren_open - 1));
    }

    if (line.find("Peak") != std::string::npos && line.find("Integral:") != std::string::npos) {
      size_t pos_colon = line.find(":");
      size_t paren_open = line.find("(");
      size_t paren_close = line.find(")");

      integral = std::stod(line.substr(pos_colon + 1, paren_open - pos_colon - 1));
      integral_err = std::stod(line.substr(paren_open + 1, paren_close - paren_open - 1));

      double min_diff = 1e9;
      int best_match = -1;
      for (size_t i = 0; i < true_peaks.size(); ++i) {
	double diff = std::abs(pos - true_peaks[i]);
	if (diff < min_diff) {
	  min_diff = diff;
	  best_match = i;
	}
      }

      peaks.emplace_back(current_hist, ++peakNum, pos, pos_err, integral, integral_err, best_match);
    }
  }

  for (auto* stream : outputs) {
    * stream << std::fixed << std::setprecision(5);
    * stream << std::setw(15) << "Histogram"
	     << std::setw(8) << "Peak#"
	     << std::setw(18) << "MatchedTruePeak"
	     << std::setw(12) << "Position"
	     << std::setw(12) << "PosErr"
	     << std::setw(12) << "Diff"
	     << std::setw(14) << "Integral"
	     << std::setw(12) << "IntErr" << "\n";

    * stream << std::string(91, '-') << "\n";

    for (const auto& peak : peaks) {
      * stream << std::setw(15) << std::get<0>(peak)
	       << std::setw(8)  << std::get<1>(peak)
	       << std::setw(18) << true_peaks[std::get<6>(peak)]
	       << std::setw(12) << std::get<2>(peak)
	       << std::setw(12) << std::get<3>(peak)
	       << std::setw(12) << std::get<2>(peak) - true_peaks[std::get<6>(peak)]
	       << std::setw(14) << std::get<4>(peak)
	       << std::setw(12) << std::get<5>(peak) << "\n";
    }
    *stream << "\n";
  }
}
