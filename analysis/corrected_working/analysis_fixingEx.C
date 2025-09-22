#include "utilities.h"

bool oxygen = false;

std::string isotope = oxygen ? "17O" : "17F";
std::string folderName = "plots_" + isotope + "/";

TString rdtCutFile = "rdtCuts_" + isotope + ".root";
TString rdtCutFileDiff = "rdtCuts_17F_with_17O_recoils.root";
TString saveFileHists = "rings_" + isotope + "_corrEx.root";
TString saveFileHistsHalfDets = "rings_" + isotope + "_half_dets_corrEx.root";
TString saveFileHistsSingle = "rings_" + isotope + "_single_dets_corrEx.root";
TString saveExHistsFile = "ex_hists_" + isotope + "_corrEx.root";
TString saveFileHists17Orecoils = "rings_17Fw17Orecoils_corrEx.root";
TString saveFileHistsHalfDets17Orecoils = "rings_17Fw17Orecoils_half_dets_corrEx.root";

TObjArray * cutList;
TObjArray * cutListDiff;
Bool_t isCutFileOpen, isCutFileOpenDiff;
int numCut, numCutDiff;

TCutG* cutG;

bool rdtgate = false;
bool diffrdtgate = false;
bool xgate = false;
bool cointimegate = false;
bool strict_cointimegate = false;
bool cointimegate_2turns = false;
bool cointimegate_3turns = false;

int n_coin = 0;
int n_coin_rdt = 0;

int recoil_n = 0;

//========correction parameters
DetGeo detGeo;

double eCorrections[24][2]; // e-correction (kinematics)

int numDet;
std::vector<double> pos;
std::vector<double> pos_half;
double Bfield;
double perpDist;
double length;
double firstPos;

bool isReaction;
double G, Z, H; // for excitation calcualtion
double mass, q;
double beta, gamma;
double alpha ;
double Et, massB;

bool plothist = true;
bool fitting = true;

int side1 = 0;
int side2 = 0;
int side3 = 0;
int side4 = 0;
int countIn[24];

std::ifstream file;

void analysis_fixingEx(){
  //================================= coinTime fit parameters

  std::vector<std::vector<double>> fitParams;
  std::string filePath = folderName + "/coinTime_x_fit_parameters.txt";
  readFitParameters(filePath.c_str(), fitParams);

  if (fitParams.empty()) {
    std::cerr << "No parameters loaded from file!" << std::endl;
    return;
  }

  //================================= rdt cut

  TFile * fCut = new TFile(rdtCutFile);
  isCutFileOpen = fCut->IsOpen();
  if(!isCutFileOpen) {
    printf( "Failed to open rdt-cutfile : %s\n" , rdtCutFile.Data());
    rdtCutFile = "";
  }
  numCut = 0 ;

  if( isCutFileOpen ){
    cutList = (TObjArray *) fCut->FindObjectAny("cutList");
    numCut = cutList->GetEntries();
    printf("=========== found %d cutG in %s \n", numCut, fCut->GetName());

    for(int i = 0; i < numCut ; i++){
      printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
	     cutList->At(i)->GetName(),
	     ((TCutG*)cutList->At(i))->GetVarX(),
	     ((TCutG*)cutList->At(i))->GetVarY(),
	     ((TCutG*)cutList->At(i))->GetN());
    }
  }

  //================================= different rdt cut (17F with 17O recoils)

  TFile * fCutDiff = new TFile(rdtCutFileDiff);
  isCutFileOpenDiff = fCutDiff->IsOpen();
  if(!isCutFileOpenDiff) {
    printf( "Failed to open rdt-cutfile : %s\n" , rdtCutFileDiff.Data());
    rdtCutFile = "";
  }
  numCutDiff = 0 ;

  if( isCutFileOpenDiff ){
    cutListDiff = (TObjArray *) fCutDiff->FindObjectAny("cutList");
    numCutDiff = cutListDiff->GetEntries();
    printf("=========== found %d cutG in %s \n", numCutDiff, fCutDiff->GetName());

    for(int i = 0; i < numCutDiff ; i++){
      printf("cut name : %s , VarX: %s, VarY: %s, numPoints: %d \n",
	     cutListDiff->At(i)->GetName(),
	     ((TCutG*)cutListDiff->At(i))->GetVarX(),
	     ((TCutG*)cutListDiff->At(i))->GetVarY(),
	     ((TCutG*)cutListDiff->At(i))->GetN());
    }
  }

  //========================================= detector Geometry
  printf("======================= loading parameters files .... \n");
  std::string detGeoFileName = "detectorGeo.txt";
  printf("loading detector geometery : %s.", detGeoFileName.c_str());

  TMacro * haha = new TMacro();
  if( haha->ReadFile(detGeoFileName.c_str()) > 0 ) {

    detGeo = LoadDetectorGeo(haha);

    PrintDetGeo(detGeo);

    Bfield = detGeo.Bfield;
    perpDist = detGeo.detPerpDist;
    firstPos = detGeo.firstPos;
    length = detGeo.detLength;
    pos = detGeo.detPos;

    printf("... done.\n");
  }else{
    printf("... fail\n");
    return;
  }

  numDet = detGeo.nDet * detGeo.mDet;

  for (auto it = pos.begin(); it != pos.end(); ++it) {
    pos_half.push_back(*it - 25.);
    pos_half.push_back(*it);
    std::cout << *it - 25. << std::endl << *it << std::endl;
  }

  //========================================= reaction parameters
  printf("loading reaction parameter.");
  file.open("reaction.dat");
  isReaction = false;
  if( file.is_open() ){
    std::string x;
    int i = 0;
    while( file >> x ){
      if( x.substr(0,2) == "//" )  continue;
      if( i == 0 ) mass = atof(x.c_str());
      if( i == 1 ) q    = atof(x.c_str());
      if( i == 2 ) beta = atof(x.c_str());
      if( i == 3 ) Et   = atof(x.c_str());
      if( i == 4 ) massB = atof(x.c_str());
      i = i + 1;
    }
    printf("................. done.\n");

    isReaction = true;
    alpha = 299.792458 * abs(Bfield) * q / TMath::TwoPi()/1000.;
    gamma = 1./TMath::Sqrt(1-beta*beta);
    G = alpha * gamma * beta * perpDist ;
    printf("============\n");
    printf("mass-b   : %f MeV/c2 \n", mass);
    printf("charge-b : %f \n", q);
    printf("E-total  : %f MeV \n", Et);
    printf("mass-B   : %f MeV/c2 \n", massB);
    printf("beta     : %f \n", beta);
    printf("B-field  : %f T \n", Bfield);
    printf("alpha    : %f MeV/mm \n", alpha);
    printf("perpDist : %f mm \n", perpDist);
    printf("G        : %f MeV \n", G);
  }else{
    printf("................. fail.\n");
    isReaction = false;
  }
  file.close();

  //========================================= e correction

  printf("loading e kinematic correction.");
  file.open("correction_e_KE_17O.dat");
  if( file.is_open() ){
    if (0){
      double a, b;
      int i = 0;
      while( file >> a >> b){
	if( i >= numDet) break;
	eCorrections[i][0] = a;  // 1/a1
	eCorrections[i][1] = b;  //  a0 , e' = e * a1 + a0
	//printf("\n%2d, e0: %9.4f, e1: %9.4f", i, eCorr[i][0], eCorr[i][1]);
	i = i + 1;
      }
      printf("....................... done.\n");
    }
  }else{
    printf("....................... fail.\n");
    for( int i = 0; i < numDet ; i++){
      eCorrections[i][0] = 1.;
      eCorrections[i][1] = 0.;
    }
    //return;
  }
  file.close();

  //================================= reading the tree

  TChain *chain = new TChain("tree");
  if (oxygen)
    chain->Add("trace_run055-066_corrected.root"); // 17O
  else
    chain->Add("trace_run022-052_corrected.root"); // 17F

  std::vector<TH1F*> correctedCoinTime;
  std::vector<TH1F*> correctedCoinTimeXgate;
  std::vector<TH1F*> correctedCoinTimeRDTCoin;
  std::vector<TH1F*> correctedCoinTimeXgateRDTCoin;
  std::vector<TH2F*> recoilDEE;

  TH1F* x_rdt_coinTime_gatedEx = new TH1F("x_rdt_coinTime_gatedEx", "Ex gated on x, recoils and coinTime", 200, -2, 12);
  TH1F* rdt_coinTime_gatedEx = new TH1F("rdt_coinTime_gatedEx", "Ex gated on recoils and coinTime", 200, -2, 12);
  TH1F* coinTime_gatedEx = new TH1F("coinTime_gatedEx", "Ex gated on coinTime", 200, -2, 12);
  TH1F* rdt_gatedEx = new TH1F("rdt_gatedEx", "Ex gated on recoils", 200, -2, 12);
  TH1F* Ex_nogates = new TH1F("Ex_nogates", "Ex, no gates", 200, -2, 12);
  TH2F* EZ_nogates = new TH2F("EZ_nogates", "e vs z, no gates", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated = new TH2F("EZ_gated", "e vs z, gated", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated_without0 = new TH2F("EZ_gated_without0", "e vs z, gated", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated_2turns = new TH2F("EZ_gated_2turns", "e vs z, gated, 2 turns", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated_3turns = new TH2F("EZ_gated_3turns", "e vs z, gated, 3 turns", 1000, -600, -200, 200, 0, 12);

  TH1F* x_diffrdt_coinTime_gatedEx = new TH1F("x_diffrdt_coinTime_gatedEx", "Ex gated on x, recoils and coinTime for another recoil cut", 200, -2, 12);
  
  std::vector<TH1F*> Ex_d; // Array to store histograms for Ex_d0, Ex_d1, ..., Ex_d5
  std::vector<TH1F*> Ex_d_half_dets;
  std::vector<TH1F*> Ex_d_strict_tc;
  std::vector<TH1F*> Ex_single; // Array to store histograms for Ex for all detectors individually
  std::vector<TH1F*> Ex_d_17Fw17Orecoils; // Array to store histograms for Ex_d0, Ex_d1, ..., Ex_d5
  std::vector<TH1F*> Ex_d_half_dets_17Fw17Orecoils;
  
  printf("Before initialising Ex_d histograms\n");

  for (int i = 0; i < 6; ++i) {
    TString histName;
    histName.Form("Ex_d%d", i);
    Ex_d.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_17Orecoil", i);
    Ex_d_17Fw17Orecoils.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_strict_tc", i);
    Ex_d_strict_tc.push_back(new TH1F(histName, histName, 200, -2, 12));
  }

  for (int i = 0; i < 12; ++i) {
    TString histName;
    histName.Form("Ex_d%d_half_dets", i);
    Ex_d_half_dets.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_half_dets_17Orecoil", i);
    Ex_d_half_dets_17Fw17Orecoils.push_back(new TH1F(histName, histName, 200, -2, 12));
  }

   printf("After initialising Ex_d histograms\n");

  for (int i = 0; i < numDet; ++i) {
    TString histName;
    histName.Form("corrected_coinTime_det%d", i);
    correctedCoinTime.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("corrected_coinTime_det%d_x_gate", i);
    correctedCoinTimeXgate.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("corrected_coinTime_det%d_rdt_coincidence", i);
    correctedCoinTimeRDTCoin.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("corrected_coinTime_det%d_x_gate_rdt_coincidence", i);
    correctedCoinTimeXgateRDTCoin.push_back(new TH1F(histName, histName, 400, -100, 200));
    histName.Form("Ex_single_det%d", i);
    Ex_single.push_back(new TH1F(histName, histName, 200, -2, 12));
  }

  for (int i = 0; i < 4; ++i) {
    TString histName;
    histName.Form("recoilDEE_%d_%d", i, i+1);
    recoilDEE.push_back(new TH2F(histName, histName, 2096, 0, 10000, 2096, 0, 10000));
  }

  // Set the branch addresses
  Float_t rdt[8], Ex, coinTime, x[24], e[24], z[24], eCorr, ExCorr, thetaCM, thetaLab;
  Int_t detID;
  Float_t coinTimeCorr;
  chain->SetBranchAddress("rdt", rdt);
  chain->SetBranchAddress("Ex", &Ex);
  chain->SetBranchAddress("coinTime", &coinTime);
  chain->SetBranchAddress("x", x);
  chain->SetBranchAddress("z", z);
  chain->SetBranchAddress("e", e);
  chain->SetBranchAddress("detID", &detID);

  for (size_t ii = 0; ii < pos_half.size(); ++ii) {
    double z_start = pos_half[ii] - 25.;
    double z_end = (ii + 1 < pos_half.size()) ? pos_half[ii] : -220.;

    std::cout << z_start << " :: " << z_end << std::endl;
  }

  // Loop over the events and apply the cut
  Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    rdtgate = false;
    diffrdtgate = false;
    xgate = false;
    cointimegate = false;
    cointimegate_2turns = false;
    cointimegate_3turns = false;
    strict_cointimegate = false;

    chain->GetEntry(i);

    if( isCutFileOpen ){
      for(int i = 0 ; i < numCut ; i++ ){
	cutG = (TCutG *)cutList->At(i) ;
	if(cutG->IsInside(rdt[2*i],rdt[2*i+1])) {
	  rdtgate = true;
	  recoil_n = 2*i;
	  break; /// only one is enough
	}
      }
    }

    if( isCutFileOpenDiff ){
      for(int i = 0 ; i < numCutDiff ; i++ ){
	cutG = (TCutG *)cutListDiff->At(i) ;
	if(cutG->IsInside(rdt[2*i],rdt[2*i+1])) {
	  diffrdtgate = true;
	  recoil_n = 2*i;
	  break; /// only one is enough
	}
      }
    }

    if ( x[detID] < 0.95 && x[detID] > -0.95)
      xgate = true;

    if ( fitParams[detID].size() == 3 ) {
      coinTimeCorr = coinTime - ( fitParams[detID][0] + fitParams[detID][1] * x[detID] + fitParams[detID][2] * TMath::Power(x[detID], 2) );
    } else if ( fitParams[detID].size() == 4 ) {
      coinTimeCorr = coinTime - ( fitParams[detID][0] + fitParams[detID][1] * x[detID] + fitParams[detID][2] * TMath::Power(x[detID], 2) + fitParams[detID][3] * TMath::Power(x[detID], 3));
    } else if ( fitParams[detID].size() == 5 ) {
      coinTimeCorr = coinTime - ( fitParams[detID][0] + fitParams[detID][1] * x[detID] + fitParams[detID][2] * TMath::Power(x[detID], 2) + fitParams[detID][3] * TMath::Power(x[detID], 3) + fitParams[detID][4] * TMath::Power(x[detID], 4));
    }
    else {
      std::cerr << "Unexpected number of fit parameters for detector " << detID << " (" << fitParams[detID].size() << "). Skipping." << std::endl;
      continue;
    }

    if (detID == 0)
      Ex = Ex + 1.22921;

    if (detID == 13) {
      Ex = Ex * 1.05429255;
    }
    if (detID == 15) {
      Ex = Ex * 0.975;
    }
    if (detID == 17) {
      Ex = Ex + 0.12693;
    }
    if (detID == 19) {
      Ex = Ex * 1.02528994;
    }
    if (detID == 22) {
      Ex = Ex - 0.06793;
    }

    EZ_nogates->Fill(z[detID], e[detID]);

    if (coinTimeCorr > -20 && coinTimeCorr < 15) {
      cointimegate = true;

      coinTime_gatedEx->Fill(Ex);
    }

    if (coinTimeCorr > -20 && coinTimeCorr < 2) {
      strict_cointimegate = true;
    }

    if (coinTimeCorr > 15 && coinTimeCorr < 40)
      cointimegate_2turns = true;

    if (coinTimeCorr > 40 && coinTimeCorr < 60)
      cointimegate_3turns = true;

    correctedCoinTime[detID]->Fill(coinTimeCorr);
    n_coin++;

    if (xgate)
      correctedCoinTimeXgate[detID]->Fill(coinTimeCorr);

    if (rdtgate) {
      correctedCoinTimeRDTCoin[detID]->Fill(coinTimeCorr);
      recoilDEE[recoil_n/2]->Fill(rdt[recoil_n], rdt[recoil_n+1]);
      n_coin_rdt++;

      rdt_gatedEx->Fill(Ex);
    }

    if (rdtgate && cointimegate)
      rdt_coinTime_gatedEx->Fill(Ex);

    if (xgate && rdtgate) {
      correctedCoinTimeXgateRDTCoin[detID]->Fill(coinTimeCorr);

      if (cointimegate) {
	
	x_rdt_coinTime_gatedEx->Fill(Ex);

	Ex_d[detID % 6]->Fill(Ex);

	Ex_single[detID]->Fill(Ex);

	EZ_gated->Fill(z[detID], e[detID]);

	if (!(detID == 0))
	  EZ_gated_without0->Fill(z[detID], e[detID]);

	if (detID <= 5)
	  side1++;
	else if (detID >= 6 && detID <= 11)
	  side2++;
	else if (detID >= 12 && detID <= 17)
	  side3++;
	else if (detID >= 18 && detID <= 23)
	  side4++;

	countIn[detID]++;

	for (size_t ii = 0; ii < pos_half.size(); ii++) {
	  double z_start = pos_half[ii] - 25.;
	  double z_end = (ii + 1 < pos_half.size()) ? pos_half[ii] : -220.;

	  // std::cout << z_start << " :: " << z_end << std::endl;

	  if (z[detID] >= z_start && z[detID] < z_end) {
	    Ex_d_half_dets[ii]->Fill(Ex);
	    break;
	  }
	}
      }

    }


    if (xgate && diffrdtgate && cointimegate && !oxygen) {
      x_diffrdt_coinTime_gatedEx->Fill(Ex);

      Ex_d_17Fw17Orecoils[detID % 6]->Fill(Ex);

      for (size_t ii = 0; ii < pos_half.size(); ii++) {
	double z_start = pos_half[ii] - 25.;
	double z_end = (ii + 1 < pos_half.size()) ? pos_half[ii] : -220.;

	if (z[detID] >= z_start && z[detID] < z_end) {
	  Ex_d_half_dets_17Fw17Orecoils[ii]->Fill(Ex);
	  break;
	}
      }
    }

    Ex_nogates->Fill(Ex);


  } // end of loop over entries

  printf("Counts in each detector side:\nside1:\t%d\tside2:\t%d\tside3:\t%d\tside4:\t%d\n", side1, side2, side3, side4);
  printf("Scaled to number of working dets:\nside1:\t%d\tside2:\t%d\tside3:\t%d\tside4:\t%d\n", side1, side2*6/4, side3, side4);
  printf("Individual sides:\n0:\t%d\t1:\t%d\t2:\t%d\t3:\t%d\t4:\t%d\t5:\t%d\n6:\t%d\t7:\t%d\t8:\t%d\t9:\t%d\t10:\t%d\t11:\t%d\n12:\t%d\t13:\t%d\t14:\t%d\t15:\t%d\t16:\t%d\t17:\t%d\n18:\t%d\t19:\t%d\t20:\t%d\t21:\t%d\t22:\t%d\t23:\t%d\n", countIn[0], countIn[1], countIn[2], countIn[3], countIn[4], countIn[5], countIn[6], countIn[7], countIn[8], countIn[9], countIn[10], countIn[11], countIn[12], countIn[13], countIn[14], countIn[15], countIn[16], countIn[17], countIn[18], countIn[19], countIn[20], countIn[21], countIn[22], countIn[23]);

  //================================= saving histograms to a root file

  TFile* outputFile = new TFile(saveFileHists, "RECREATE");

  for (int ii = 0; ii < 6; ++ii)
    Ex_d[ii]->Write();

  outputFile->Close();

  TFile* outputFileHalfDets = new TFile(saveFileHistsHalfDets, "RECREATE");

  for (int ii = 0; ii < 12; ++ii)
    Ex_d_half_dets[ii]->Write();

  outputFileHalfDets->Close();

   TFile* outputFileSingle = new TFile(saveFileHistsSingle, "RECREATE");

  for (int ii = 0; ii < 24; ++ii)
    Ex_single[ii]->Write();

  outputFileSingle->Close();

  if (!oxygen) {
    TFile* outputFile17Orecoils = new TFile(saveFileHists17Orecoils, "RECREATE");

    for (int ii = 0; ii < 6; ++ii)
      Ex_d_17Fw17Orecoils[ii]->Write();

    outputFile17Orecoils->Close();

    TFile* outputFileHalfDets17Orecoils = new TFile(saveFileHistsHalfDets17Orecoils, "RECREATE");

    for (int ii = 0; ii < 12; ++ii)
      Ex_d_half_dets_17Fw17Orecoils[ii]->Write();

    outputFileHalfDets17Orecoils->Close();
  }
}
