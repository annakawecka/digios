#include <TH2F.h>
#include <TFile.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TString.h>
#include <TObjArray.h>

void rdt_stats(){
  TChain * chain = new TChain("tree");
  chain->Add("gen_run022-052.root");

  chain->GetListOfFiles()->Print();

  TCanvas * cRDTStats = new TCanvas("RDTStats", "RDT Stats", 100, 100, 800, 800);

  TString varX, varY, expression[4];

  std::ofstream outFile("plots/rdt_entry_counts.txt");  // Open file to save the counts

  // Check if the file is opened successfully
  if (!outFile.is_open()) {
    std::cerr << "Error opening file to save counts!" << std::endl;
    return;
  }

  for (Int_t i = 0; i < 4; i++) {

      varX.Form("rdt[%d]",2 * i);
      varY.Form("rdt[%d]",2 * i + 1);

      TString histName = Form("h_%d", i);

      expression[i].Form("%s:%s>>%s(2096, 0, 10000, 2096, 0, 10000)", 
			 varY.Data(),
			 varX.Data(),
			 histName.Data());

      chain->Draw(expression[i], "", "colz");
     
      TH2F *hist = (TH2F*)gDirectory->Get(histName.Data());
      if (hist) {
	int entries = hist->GetEntries();
	printf("Number of entries in %s vs %s plot: %d\n", varY.Data(), varX.Data(), entries);

	// Write the count to the output file
	outFile << "Entries in " << varY.Data() << " vs " << varX.Data() << ": " << entries << "\n";
      }

      Int_t nEntriesX = chain->GetEntries(varX.Data());  // Get number of entries for this variable
      Int_t nEntriesY = chain->GetEntries(varY.Data());  // Get number of entries for this variable
        
      // Print the result to console
      std::cout << "Entries for " << varX.Data() << ": " << nEntriesX << std::endl;
      std::cout << "Entries for " << varY.Data() << ": " << nEntriesY << std::endl;

      // Save the result to the text file
      outFile << "Entries for " << varX.Data() << ": " << nEntriesX << std::endl;
      outFile << "Entries for " << varY.Data() << ": " << nEntriesY << std::endl;

      // Update canvas to display the plot
      cRDTStats->Modified();
      cRDTStats->Update();
   }

  
}
