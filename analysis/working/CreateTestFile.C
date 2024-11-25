#define CreateTestFile_cxx
#include "CreateTestFile.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void CreateTestFile::Loop()
{
//   In a ROOT session, you can do:
//      root> .L CreateTestFile.C
//      root> CreateTestFile t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //#################################################################### initialization
      eTemp    = TMath::QuietNaN();
      xTemp    = TMath::QuietNaN();
      zTemp    = TMath::QuietNaN();
   
      detIDTemp = -4;
      hitIDTemp = -4;
      multiHitTemp = 0;
   
      coinTimeTemp = TMath::QuietNaN();

      //#################################################################### Get Tree
      eventIDTemp = eventID;

      b_e->GetEntry(ientry,0);
      //b_Ring->GetEntry(entry,0);
      b_xf->GetEntry(ientry,0);
      b_xn->GetEntry(ientry,0);
      b_rdt->GetEntry(ientry,0);
      //b_EnergyTimestamp->GetEntry(entry,0);
      //b_RDTTimestamp->GetEntry(entry,0);
      b_coinTime->GetEntry(ientry,0);

      //=========== gate
   
      bool rejRDT1 = true; if( isRDTCutExist && cut[0] && cut[0]->IsInside( rdt[0], rdt[1] )) rejRDT1 = false;
      bool rejRDT2 = true; if( isRDTCutExist && cut[1] && cut[1]->IsInside( rdt[2], rdt[3] )) rejRDT2 = false;
      bool rejRDT3 = true; if( isRDTCutExist && cut[2] && cut[2]->IsInside( rdt[4], rdt[5] )) rejRDT3 = false;
      bool rejRDT4 = true; if( isRDTCutExist && cut[3] && cut[3]->IsInside( rdt[6], rdt[7] )) rejRDT4 = false;
   
      if( !isRDTCutExist ){
	rejRDT1 = false;
	rejRDT2 = false;
	rejRDT3 = false;
	rejRDT4 = false;
      }

      if( rejRDT1 && rejRDT2 && rejRDT3 && rejRDT4) continue;

      for(int idet = 0 ; idet < numDet; idet++){
      
	if( e[idet] < 0.010 ) continue;
	if( TMath::IsNaN(e[idet]) ) continue;
	if( TMath::IsNaN(xf[idet]) && TMath::IsNaN(xn[idet])  ) continue;

	if( e[idet] == 0 ) continue;
	if( xf[idet] == 0 && xn[idet] == 0 ) continue;

	eTemp = e[idet];
	xTemp = x[idet];
	zTemp = z[idet];

	if( !TMath::IsNaN(zTemp) ) {
	  multiHitTemp++;
	  detIDTemp = idet;
	}

	
      }

      if (detIDTemp < 0 || detIDTemp >= coinTimeFitParams.size()) {
	std::cerr << "Invalid detector ID: " << detIDTemp << std::endl;
	continue;
      }
      
      if (!TMath::IsNaN(coinTime)) {
	if ( coinTimeFitParams[detIDTemp].size() == 3 ) {
	  coinTimeTemp = coinTime - ( coinTimeFitParams[detIDTemp][0] + coinTimeFitParams[detIDTemp][1] * xTemp + coinTimeFitParams[detIDTemp][2] * TMath::Power(xTemp, 2) );
	} else if ( coinTimeFitParams[detIDTemp].size() == 4 ) {
	  coinTimeTemp = coinTime - ( coinTimeFitParams[detIDTemp][0] + coinTimeFitParams[detIDTemp][1] * xTemp + coinTimeFitParams[detIDTemp][2] * TMath::Power(xTemp, 2) + coinTimeFitParams[detIDTemp][3] * TMath::Power(xTemp, 3));
	} else if ( coinTimeFitParams[detIDTemp].size() == 5 ) {
	  coinTimeTemp = coinTime - ( coinTimeFitParams[detIDTemp][0] + coinTimeFitParams[detIDTemp][1] * xTemp + coinTimeFitParams[detIDTemp][2] * TMath::Power(xTemp, 2) + coinTimeFitParams[detIDTemp][3] * TMath::Power(xTemp, 3) + coinTimeFitParams[detIDTemp][4] * TMath::Power(xTemp, 4));
	}
	else {
	  std::cerr << "Unexpected number of fit parameters for detector " << detIDTemp << " (" << coinTimeFitParams[detIDTemp].size() << "). Skipping." << std::endl;
	  continue;
	}
      
      }

      if (coinTimeTemp < -20. || coinTimeTemp > 15.)
	continue;

      ExTemp = Ex;

      saveFile->cd(); //set focus on this file
      newTree->Fill();
   }

   saveFile->cd(); //set focus on this file
   newTree->Write(); 
   saveFile->Close();
}
