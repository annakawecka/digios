#include "utilities.h"

bool oxygen = true;

std::string isotope = oxygen ? "17O" : "17F";
std::string folderName = "plots_" + isotope + "_up_line/";

TString rdtCutFile = "rdtCuts_" + isotope + "_up_line.root";
TString rdtCutFileDiff = "rdtCuts_17F_with_17O_recoils.root";
TString saveFileHists = "rings_" + isotope + "_up_line.root";
TString saveFileHistsHalfDets = "rings_" + isotope + "_half_dets_up_line.root";
TString saveFileHistsStrictTC = "rings_" + isotope + "_200bins_strict_tc_up_line.root";
TString saveFileHistsSingle = "rings_" + isotope + "_single_dets_up_line.root";
TString saveExHistsFile = "ex_hists_" + isotope + "_up_line.root";
TString saveFileHists17Orecoils = "rings_17Fw17Orecoils_up_line.root";
TString saveFileHistsHalfDets17Orecoils = "rings_17Fw17Orecoils_half_dets_up_line.root";

TString saveFileHistsMoreBins = "rings_" + isotope + "_morebins_up_line.root";
TString saveFileHistsHalfDetsMoreBins = "rings_" + isotope + "_half_dets_morebins_up_line.root";
TString saveFileHists17OrecoilsMoreBins = "rings_17Fw17Orecoils_morebins_up_line.root";
TString saveFileHistsHalfDets17OrecoilsMoreBins = "rings_17Fw17Orecoils_half_dets_morebins_up_line.root";

TObjArray * cutList;
TObjArray * cutListDiff;
Bool_t isCutFileOpen, isCutFileOpenDiff;
int numCut, numCutDiff;

TCutG* cutG;

bool rdtgate = false;
bool rdtgate0 = false;
bool rdtgate1 = false;
bool rdtgate2 = false;
bool rdtgate3 = false;
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

void SetThesisStyle() {
  gStyle->SetOptStat(0);                 // Hide stat box
  gStyle->SetOptFit(0);                  // Hide fit box unless wanted

  gStyle->SetPalette(kRainBow);          // Nicer color palette (choose one)

  // Canvas
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);

  // Pads
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Frame
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameBorderMode(0);

  // Fonts
  gStyle->SetTextFont(42);               // Helvetica (thesis-friendly)
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");

  // Margins
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.05);

  // Marker & line styles
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.25);

  gStyle->SetImageScaling(3.0);
}

void StyleHistogram(TH1 *h, const char *mainTitle,
                    const char *xTitle, const char *yTitle, int width = 2)
{
    if (!h) return;

    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    h->SetTitle("");

    h->SetLineWidth(width);

    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.037);
    h->GetYaxis()->SetLabelSize(0.037);
    h->GetXaxis()->SetTickLength(0.02);
    h->GetYaxis()->SetTickLength(0.02);

    h->GetXaxis()->SetTitleOffset(1.25);
    h->GetYaxis()->SetTitleOffset(1.35);
    h->GetXaxis()->SetLabelOffset(0.02);
    h->GetYaxis()->SetLabelOffset(0.02);

    //h->GetXaxis()->SetNdivisions(505);
    //h->GetYaxis()->SetNdivisions(505);

    gPad->SetTopMargin(0.12);

    TLatex *latex = new TLatex();
    latex->SetTextFont(42);
    latex->SetTextSize(0.05);
    latex->SetTextAlign(22);
    latex->DrawLatexNDC(0.5, 0.95, mainTitle);
}


void analysis(){

  SetThesisStyle();
  
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
  TH2F* EZ_nogates_without0 = new TH2F("EZ_nogates_without0", "e vs z, no gates", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated = new TH2F("EZ_gated", "e vs z, gated", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated_without0 = new TH2F("EZ_gated_without0", "e vs z, gated", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated_2turns = new TH2F("EZ_gated_2turns", "e vs z, gated, 2 turns", 1000, -600, -200, 200, 0, 12);
  TH2F* EZ_gated_3turns = new TH2F("EZ_gated_3turns", "e vs z, gated, 3 turns", 1000, -600, -200, 200, 0, 12);

  TH1F* x_diffrdt_coinTime_gatedEx = new TH1F("x_diffrdt_coinTime_gatedEx", "Ex gated on x, recoils and coinTime for another recoil cut", 200, -2, 12);
  
  std::vector<TH1F*> Ex_d; // Array to store histograms for Ex_d0, Ex_d1, ..., Ex_d5
  std::vector<TH1F*> Ex_d_morebins;
  std::vector<TH1F*> Ex_d_half_dets;
  std::vector<TH1F*> Ex_d_half_dets_morebins;
  std::vector<TH1F*> Ex_d_strict_tc;
  std::vector<TH1F*> Ex_single; // Array to store histograms for Ex for all detectors individually
  std::vector<TH1F*> Ex_single_rdt0;
  std::vector<TH1F*> Ex_single_rdt1;
  std::vector<TH1F*> Ex_single_rdt2;
  std::vector<TH1F*> Ex_single_rdt3;
  std::vector<TH1F*> Ex_d_17Fw17Orecoils; // Array to store histograms for Ex_d0, Ex_d1, ..., Ex_d5
  std::vector<TH1F*> Ex_d_17Fw17Orecoils_morebins;
  std::vector<TH1F*> Ex_d_half_dets_17Fw17Orecoils;
  std::vector<TH1F*> Ex_d_half_dets_17Fw17Orecoils_morebins;
  
  printf("Before initialising Ex_d histograms\n");

  for (int i = 0; i < 6; ++i) {
    TString histName;
    histName.Form("Ex_d%d", i);
    Ex_d.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_morebins", i);
    Ex_d_morebins.push_back(new TH1F(histName, histName, 600, -2, 12));
    histName.Form("Ex_d%d_17Orecoil", i);
    Ex_d_17Fw17Orecoils.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_17Orecoil_morebins", i);
    Ex_d_17Fw17Orecoils_morebins.push_back(new TH1F(histName, histName, 600, -2, 12));
    histName.Form("Ex_d%d_strict_tc", i);
    Ex_d_strict_tc.push_back(new TH1F(histName, histName, 200, -2, 12));
  }

  for (int i = 0; i < 12; ++i) {
    TString histName;
    histName.Form("Ex_d%d_half_dets", i);
    Ex_d_half_dets.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_half_dets_17Orecoil", i);
    Ex_d_half_dets_17Fw17Orecoils.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_d%d_half_dets_morebins", i);
    Ex_d_half_dets_morebins.push_back(new TH1F(histName, histName, 600, -2, 12));
    histName.Form("Ex_d%d_half_dets_17Orecoil_morebins", i);
    Ex_d_half_dets_17Fw17Orecoils_morebins.push_back(new TH1F(histName, histName, 600, -2, 12));
  }

   printf("After initialising Ex_d histograms\n");

  for (int i = 0; i < numDet; ++i) {
    TString histName;
    histName.Form("corrected_coinTime_det%d", i);
    correctedCoinTime.push_back(new TH1F(histName, histName, 300, -100, 200));
    histName.Form("corrected_coinTime_det%d_x_gate", i);
    correctedCoinTimeXgate.push_back(new TH1F(histName, histName, 300, -100, 200));
    histName.Form("corrected_coinTime_det%d_rdt_coincidence", i);
    correctedCoinTimeRDTCoin.push_back(new TH1F(histName, histName, 300, -100, 200));
    histName.Form("corrected_coinTime_det%d_x_gate_rdt_coincidence", i);
    correctedCoinTimeXgateRDTCoin.push_back(new TH1F(histName, histName, 300, -100, 200));
    histName.Form("Ex_single_det%d", i);
    Ex_single.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_single_det%d_rdt0", i);
    Ex_single_rdt0.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_single_det%d_rdt1", i);
    Ex_single_rdt1.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_single_det%d_rdt2", i);
    Ex_single_rdt2.push_back(new TH1F(histName, histName, 200, -2, 12));
    histName.Form("Ex_single_det%d_rdt3", i);
    Ex_single_rdt3.push_back(new TH1F(histName, histName, 200, -2, 12));
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
    rdtgate0 = false;
    rdtgate1 = false;
    rdtgate2 = false;
    rdtgate3 = false;
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

    if( isCutFileOpen ){
      for(int i = 0 ; i < numCut ; i++ ){
	cutG = (TCutG *)cutList->At(i) ;
	if(cutG->IsInside(rdt[2*i],rdt[2*i+1])) {
	  if (i == 0)
	    rdtgate0 = true;
	  if (i == 1)
	    rdtgate1 = true;
	  if (i == 2)
	    rdtgate2 = true;
	  if (i == 3)
	    rdtgate3 = true;
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

    EZ_nogates->Fill(z[detID], e[detID]);

    if (!(detID == 0))
      EZ_nogates_without0->Fill(z[detID], e[detID]);

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

      //if (cointimegate && !(detID == 0)) {
      if (cointimegate) {
	x_rdt_coinTime_gatedEx->Fill(Ex);

	Ex_d[detID % 6]->Fill(Ex);
	Ex_d_morebins[detID % 6]->Fill(Ex);

	Ex_single[detID]->Fill(Ex);

	if (rdtgate0)
	  Ex_single_rdt0[detID]->Fill(Ex);
	if (rdtgate1)
	  Ex_single_rdt1[detID]->Fill(Ex);
	if (rdtgate2)
	  Ex_single_rdt2[detID]->Fill(Ex);
	if (rdtgate3)
	  Ex_single_rdt3[detID]->Fill(Ex);

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
	    Ex_d_half_dets_morebins[ii]->Fill(Ex);
	    break;
	  }
	}
      }

      if (strict_cointimegate && !(detID == 0))
	Ex_d_strict_tc[detID % 6]->Fill(Ex);

      if (cointimegate_2turns && !(detID == 0))
	EZ_gated_2turns->Fill(z[detID], e[detID]);

      if (cointimegate_3turns && !(detID == 0))
	EZ_gated_3turns->Fill(z[detID], e[detID]);
    }

    if (xgate && diffrdtgate && cointimegate) {
      x_diffrdt_coinTime_gatedEx->Fill(Ex);
    }

    if (xgate && diffrdtgate && cointimegate && !oxygen) {
      x_diffrdt_coinTime_gatedEx->Fill(Ex);

      Ex_d_17Fw17Orecoils[detID % 6]->Fill(Ex);
      Ex_d_17Fw17Orecoils_morebins[detID % 6]->Fill(Ex);

      for (size_t ii = 0; ii < pos_half.size(); ii++) {
	double z_start = pos_half[ii] - 25.;
	double z_end = (ii + 1 < pos_half.size()) ? pos_half[ii] : -220.;

	if (z[detID] >= z_start && z[detID] < z_end) {
	  Ex_d_half_dets_17Fw17Orecoils[ii]->Fill(Ex);
	  Ex_d_half_dets_17Fw17Orecoils_morebins[ii]->Fill(Ex);
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

  TFile* outputFileStrictTC = new TFile(saveFileHistsStrictTC, "RECREATE");

  for (int ii = 0; ii < 6; ++ii)
    Ex_d_strict_tc[ii]->Write();

  outputFileStrictTC->Close();

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

    TFile* outputFile17OrecoilsMoreBins = new TFile(saveFileHists17OrecoilsMoreBins, "RECREATE");

    for (int ii = 0; ii < 6; ++ii)
      Ex_d_17Fw17Orecoils_morebins[ii]->Write();

    outputFile17OrecoilsMoreBins->Close();

    TFile* outputFileHalfDets17OrecoilsMoreBins = new TFile(saveFileHistsHalfDets17OrecoilsMoreBins, "RECREATE");

    for (int ii = 0; ii < 12; ++ii)
      Ex_d_half_dets_17Fw17Orecoils_morebins[ii]->Write();

    outputFileHalfDets17OrecoilsMoreBins->Close();
  }

  TFile* outputFileMoreBins = new TFile(saveFileHistsMoreBins, "RECREATE");

  for (int ii = 0; ii < 6; ++ii)
    Ex_d_morebins[ii]->Write();

  outputFileMoreBins->Close();

  TFile* outputFileHalfDetsMoreBins = new TFile(saveFileHistsHalfDetsMoreBins, "RECREATE");

  for (int ii = 0; ii < 12; ++ii)
    Ex_d_half_dets_morebins[ii]->Write();

  outputFileHalfDetsMoreBins->Close();

  //================================= fitting 17O
  if (fitting) {

    printf("Start of fitting\n");

    TCanvas *cExdet = new TCanvas("cExdet", "Ex for different detector rings", 1000, 600);
    cExdet->Divide(3, 2);

    for (int i = 0; i < 6; ++i) {
      cExdet->cd(i + 1);
      Ex_d[i]->Draw();
      StyleHistogram(Ex_d[i],
		     Form("E_{x}, ring %d", i),
		     Form("E_{x} (MeV)"),
		     Form("Counts / 70 keV"),
		     1);
    }

    cExdet->SaveAs((folderName + "/Ex_for_rings.png").c_str());
    cExdet->SaveAs((folderName + "/Ex_for_rings.pdf").c_str());

    //for (int detectorId = 0; detectorId <= 5; detectorId++) {
    //  fitSpectra(Ex_d[detectorId], detectorId);
    //} // the 'real' fitting part is done in a separate root script
    
  }

  //================================= plotting and saving histograms

  if (plothist)  {

    TCanvas *cAllDetectors = new TCanvas("cAllDetectors", "Corrected CoinTime for All Detectors", 1200, 800);
    cAllDetectors->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectors->cd(i + 1);
      correctedCoinTime[i]->Draw();
      StyleHistogram(correctedCoinTime[i],
		     Form("Detector %d", i),
		     Form("coincidence time (ns)"),
		     Form("Counts / 1 ns"),
		     1);
    }

    cAllDetectors->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors.png").c_str());
    cAllDetectors->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors.pdf").c_str());

    TCanvas *cAllDetectorsXgate = new TCanvas("cAllDetectorsXgate", "Corrected coincidence time for All Detectors with x gate", 1200, 800);
    cAllDetectorsXgate->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectorsXgate->cd(i + 1);
      correctedCoinTimeXgate[i]->Draw();
      StyleHistogram(correctedCoinTimeXgate[i],
		     Form("Detector %d", i),
		     Form("coincidence time (ns)"),
		     Form("Counts / 1 ns"),
		     1);
    }

    cAllDetectorsXgate->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors_Xgate.png").c_str());
    cAllDetectorsXgate->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors_Xgate.pdf").c_str());

    TCanvas *cAllDetectorsRDT = new TCanvas("cAllDetectorsRDT", "Corrected coincidence time for All Detectors (RDT coin)", 1200, 800);
    cAllDetectorsRDT->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectorsRDT->cd(i + 1);
      correctedCoinTimeRDTCoin[i]->Draw();
      StyleHistogram(correctedCoinTimeRDTCoin[i],
		     Form("Detector %d", i),
		     Form("coincidence time (ns)"),
		     Form("Counts / 1 ns"),
		     1);
    }

    cAllDetectorsRDT->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors_RDTCoin.png").c_str());
    cAllDetectorsRDT->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors_RDTCoin.pdf").c_str());

    TCanvas *cAllDetectorsXgateRDT = new TCanvas("cAllDetectorsXgateRDT", "Corrected coincidence time for All Detectors (RDT coin && x gate)", 1200, 800);
    cAllDetectorsXgateRDT->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectorsXgateRDT->cd(i + 1);
      correctedCoinTimeXgateRDTCoin[i]->Draw();
      StyleHistogram(correctedCoinTimeXgateRDTCoin[i],
		     Form("Detector %d", i),
		     Form("coinTime (ns)"),
		     Form("Counts / 1 ns"),
		     1);
    }

    cAllDetectorsRDT->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors_XgateRDTCoin.png").c_str());
    cAllDetectorsRDT->SaveAs((folderName + "/Corrected_CoinTime_AllDetectors_XgateRDTCoin.pdf").c_str());

    std::cout << "coin: " << n_coin << " coin_rdt: " << n_coin_rdt << std::endl;

    TCanvas *cRDT = new TCanvas("cRDT", "RDT dE-E gated", 1200, 300);
    cRDT->Divide(4);

    for (int i = 0; i < 4; ++i) {
      cRDT->cd(i + 1);
      recoilDEE[i]->Draw();
    }

    cRDT->SaveAs((folderName + "/RDT_DEE_gated.png").c_str());
  

    TH1F *combinedHist = (TH1F*)correctedCoinTimeXgateRDTCoin[0]->Clone("combinedHist");
    combinedHist->Reset();
    for (int i = 0; i < numDet; ++i) {
      combinedHist->Add(correctedCoinTimeXgateRDTCoin[i]);
    }

    TCanvas *cCombined = new TCanvas("cCombined", "Combined Corrected CoinTime with recoil && x gate", 800, 600);
    combinedHist->Draw();
    StyleHistogram(combinedHist,
		   Form("Corrected coincidence time"),
		   Form("time difference (ns)"),
		   Form("Counts / 1 ns"));
    cCombined->SaveAs((folderName + "/Combined_Corrected_CoinTime_XgateRDTCoin.png").c_str());
    cCombined->SaveAs((folderName + "/Combined_Corrected_CoinTime_XgateRDTCoin.pdf").c_str());
    cCombined->SaveAs((folderName + "/Combined_Corrected_CoinTime_XgateRDTCoin.root").c_str());
    

    TCanvas *cGatedEx = new TCanvas("cGatedEx", "Ex gated on x, recoils and coinTime", 800, 600);
    x_rdt_coinTime_gatedEx->Draw();
    StyleHistogram(x_rdt_coinTime_gatedEx,
		   Form("E_{x} spectrum, gated on position, RDT and coincidence time"),
		   Form("E_{x} (MeV)"),
		   Form("Counts / 70 keV"));
    if (oxygen == false){
      gPad->Update();
      double x_line = 5.6071; // 5.6 MeV
      double y_min = gPad->GetUymin();
      double y_max = gPad->GetUymax();

      TLine *vline = new TLine(x_line, y_min, x_line, y_max);
      vline->SetLineColor(kBlack);
      vline->SetLineStyle(2); // kreskowana
      vline->SetLineWidth(4);
      vline->Draw("same");

      TLatex *label = new TLatex(
				 x_line + 0.15,              // small x-offset to avoid overlap
				 y_max * 0.85,               // near top of spectrum
				 "S_{p}(^{18}F)"
				 );

      label->SetTextFont(42);
      label->SetTextSize(0.04);
      label->SetTextColor(kBlack);
      label->SetTextAlign(12);        // left, vertically centered
      label->Draw();
    }
    cGatedEx->SaveAs((folderName + "/Ex_x_recoil_coinTime_gated.png").c_str());
    cGatedEx->SaveAs((folderName + "/Ex_x_recoil_coinTime_gated.pdf").c_str());
    x_rdt_coinTime_gatedEx->SaveAs((folderName + "/Ex_x_recoil_coinTime_gated.root").c_str());

    TCanvas *cGatedExDiffRdt = new TCanvas("cGatedExDiffRdt", "Ex gated on x, recoils and coinTime, diff rdt", 800, 600);
    x_diffrdt_coinTime_gatedEx->Draw();
    StyleHistogram(x_diffrdt_coinTime_gatedEx,
		   Form("E_{x} spectrum, gated on position, RDT and coincidence time"),
		   Form("E_{x} (MeV)"),
		   Form("Counts / 70 keV"));
    if (oxygen == false){
      gPad->Update();
      double x_line = 5.6071; // 5.6 MeV
      double y_min = gPad->GetUymin();
      double y_max = gPad->GetUymax();

      TLine *vline = new TLine(x_line, y_min, x_line, y_max);
      vline->SetLineColor(kBlack);
      vline->SetLineStyle(2); // kreskowana
      vline->SetLineWidth(4);
      vline->Draw("same");

       TLatex *label = new TLatex(
				 x_line - 1.6,              // small x-offset to avoid overlap
				 y_max * 0.85,               // near top of spectrum
				 "S_{p}(^{18}F)"
				 );

      label->SetTextFont(42);
      label->SetTextSize(0.04);
      label->SetTextColor(kBlack);
      label->SetTextAlign(12);        // left, vertically centered
      label->Draw();
    }
    cGatedExDiffRdt->SaveAs((folderName + "/Ex_x_diffrecoil_coinTime_gated.png").c_str());
    cGatedExDiffRdt->SaveAs((folderName + "/Ex_x_diffrecoil_coinTime_gated.pdf").c_str());
    x_diffrdt_coinTime_gatedEx->SaveAs((folderName + "/Ex_x_diffrecoil_coinTime_gated.root").c_str());

    TCanvas *cRDTCoinTimeGatedEx = new TCanvas("cRDTCoinTimeGatedEx", "Ex gated on recoils and coinTime", 800, 600);
    rdt_coinTime_gatedEx->Draw();
    StyleHistogram(rdt_coinTime_gatedEx,
		   Form("E_{x} spectrum, gated on RDT and coincidence time"),
		   Form("E_{x} (MeV)"),
		   Form("Counts / 70 keV"));
    if (oxygen == false){
      gPad->Update();
      double x_line = 5.6071; // 5.6 MeV
      double y_min = gPad->GetUymin();
      double y_max = gPad->GetUymax();

      TLine *vline = new TLine(x_line, y_min, x_line, y_max);
      vline->SetLineColor(kBlack);
      vline->SetLineStyle(2); // kreskowana
      vline->SetLineWidth(4);
      vline->Draw("same");

       TLatex *label = new TLatex(
				 x_line + 0.15,              // small x-offset to avoid overlap
				 y_max * 0.85,               // near top of spectrum
				 "S_{p}(^{18}F)"
				 );

      label->SetTextFont(42);
      label->SetTextSize(0.04);
      label->SetTextColor(kBlack);
      label->SetTextAlign(12);        // left, vertically centered
      label->Draw();
    }
    cRDTCoinTimeGatedEx->SaveAs((folderName + "/Ex_recoil_coinTime_gated.png").c_str());
    cRDTCoinTimeGatedEx->SaveAs((folderName + "/Ex_recoil_coinTime_gated.pdf").c_str());

    TCanvas *cCoinTimeGatedEx = new TCanvas("cCoinTimeGatedEx", "Ex gated on coinTime", 800, 600);
    coinTime_gatedEx->Draw();
    StyleHistogram(coinTime_gatedEx,
		   Form("E_{x} spectrum, gated on coincidence time"),
		   Form("E_{x} (MeV)"),
		   Form("Counts / 70 keV"));
    if (oxygen == false){
      gPad->Update();
      double x_line = 5.6071; // 5.6 MeV
      double y_min = gPad->GetUymin();
      double y_max = gPad->GetUymax();

      TLine *vline = new TLine(x_line, y_min, x_line, y_max);
      vline->SetLineColor(kBlack);
      vline->SetLineStyle(2); // kreskowana
      vline->SetLineWidth(4);
      vline->Draw("same");

       TLatex *label = new TLatex(
				 x_line + 0.15,              // small x-offset to avoid overlap
				 y_max * 0.85,               // near top of spectrum
				 "S_{p}(^{18}F)"
				 );

      label->SetTextFont(42);
      label->SetTextSize(0.04);
      label->SetTextColor(kBlack);
      label->SetTextAlign(12);        // left, vertically centered
      label->Draw();
    }
    cCoinTimeGatedEx->SaveAs((folderName + "/Ex_coinTime_gated.png").c_str());
    cCoinTimeGatedEx->SaveAs((folderName + "/Ex_coinTime_gated.pdf").c_str());

    TCanvas *cRDTGatedEx = new TCanvas("cRDTGatedEx", "Ex gated on rdt", 800, 600);
    rdt_gatedEx->Draw();
    StyleHistogram(rdt_gatedEx,
		   Form("E_{x} spectrum, gated on RDT"),
		   Form("E_{x} (MeV)"),
		   Form("Counts / 70 keV"));
    cRDTGatedEx->SaveAs((folderName + "/Ex_rdt_gated.png").c_str());
    cRDTGatedEx->SaveAs((folderName + "/Ex_rdt_gated.pdf").c_str());

    TCanvas *cExnogates = new TCanvas("cExnogates", "Ex, no gates", 800, 600);
    Ex_nogates->Draw();
    StyleHistogram(Ex_nogates,
		   Form("E_{x} spectrum, no gates"),
		   Form("E_{x} (MeV)"),
		   Form("Counts / 70 keV"));
    cExnogates->SaveAs((folderName + "/Ex_nogates.png").c_str());
    cExnogates->SaveAs((folderName + "/Ex_nogates.pdf").c_str());

    TCanvas *cEZnogates = new TCanvas("cEZnogates", "e vs z, no gates", 800, 600);
    EZ_nogates->Draw();
    StyleHistogram(EZ_nogates,
		   Form("E vs z, no gates"),
		   Form("z (mm)"),
		   Form("E (MeV)"));
    cEZnogates->SaveAs((folderName + "/EZ_nogates.png").c_str());
    cEZnogates->SaveAs((folderName + "/EZ_nogates.pdf").c_str());

    TCanvas *cEZnogatesNo0 = new TCanvas("cEZnogateswithout0", "e vs z, no gates", 800, 600);
    EZ_nogates_without0->Draw();
    StyleHistogram(EZ_nogates_without0,
		   Form("E vs z, no gates"),
		   Form("z (mm)"),
		   Form("E (MeV)"));
    cEZnogatesNo0->SaveAs((folderName + "/EZ_nogates_without0.pdf").c_str());
    cEZnogatesNo0->SaveAs((folderName + "/EZ_nogates_without0.png").c_str());

    gStyle->SetOptStat(000);

    TCanvas *cEZgated = new TCanvas("cEZgated", "e vs z, gated", 800, 600);
    EZ_gated->Draw();
    StyleHistogram(EZ_gated,
		   Form("E vs z, gated"),
		   Form("z (mm)"),
		   Form("E (MeV)"));
    cEZgated->SaveAs((folderName + "/EZ_gated.png").c_str());
    cEZgated->SaveAs((folderName + "/EZ_gated.pdf").c_str());

    TCanvas *cEZgatedNo0 = new TCanvas("cEZgatedNo0", "e vs z, gated", 800, 600);
    EZ_gated_without0->Draw();
    StyleHistogram(EZ_gated_without0,
		   Form("E vs z, gated"),
		   Form("z (mm)"),
		   Form("E (MeV)"));
    cEZgatedNo0->SaveAs((folderName + "/EZ_gated_without_det0.png").c_str());
    cEZgatedNo0->SaveAs((folderName + "/EZ_gated_without_det0.pdf").c_str());

    TCanvas *cEZgated2turns = new TCanvas("cEZgated2turns", "e vs z, gated, 2 turns", 800, 600);
    EZ_gated_2turns->Draw();
    StyleHistogram(EZ_gated_2turns,
		   Form("E vs z, gated, two turns"),
		   Form("z (mm)"),
		   Form("E (MeV)"));
    cEZgated2turns->SaveAs((folderName + "/EZ_gated_2turns.png").c_str());
    cEZgated2turns->SaveAs((folderName + "/EZ_gated_2turns.pdf").c_str());

    TCanvas *cEZgated3turns = new TCanvas("cEZgated3turns", "e vs z, gated, 3 turns", 800, 600);
    EZ_gated_3turns->Draw();
    StyleHistogram(EZ_gated_3turns,
		   Form("E vs z, gated, three turns"),
		   Form("z (mm)"),
		   Form("E (MeV)"));
    cEZgated3turns->SaveAs((folderName + "/EZ_gated_3turns.png").c_str());
    cEZgated3turns->SaveAs((folderName + "/EZ_gated_3turns.pdf").c_str());

    TCanvas *cExSingleDets = new TCanvas("cExSingleDets", "Ex for each detector", 1200, 800);
    cExSingleDets->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cExSingleDets->cd(i + 1);
      Ex_single[i]->Draw();
      StyleHistogram(Ex_single[i],
		     Form("E_{x} spectrum, gated"),
		     Form("E_{x} (MeV)"),
		     Form("Counts / 70 keV"),
		     1);
    }

    cExSingleDets->SaveAs((folderName + "/Ex_each_detector.png").c_str());
    cExSingleDets->SaveAs((folderName + "/Ex_each_detector.pdf").c_str());

    TCanvas *cExSingleDetsRDT0 = new TCanvas("cExSingleDetsRDT0", "Ex for each detector, rdt0", 1200, 800);
    cExSingleDetsRDT0->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cExSingleDetsRDT0->cd(i + 1);
      Ex_single_rdt0[i]->Draw();
      StyleHistogram(Ex_single_rdt0[i],
		     Form("E_{x} spectrum, gated on RDT0"),
		     Form("E_{x} (MeV)"),
		     Form("Counts / 70 keV"),
		     1);
    }

    cExSingleDetsRDT0->SaveAs((folderName + "/Ex_each_detector_rdt0.png").c_str());
    cExSingleDetsRDT0->SaveAs((folderName + "/Ex_each_detector_rdt0.pdf").c_str());

    TCanvas *cExSingleDetsRDT1 = new TCanvas("cExSingleDetsRDT1", "Ex for each detector, rdt1", 1200, 800);
    cExSingleDetsRDT1->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cExSingleDetsRDT1->cd(i + 1);
      Ex_single_rdt1[i]->Draw();
      StyleHistogram(Ex_single_rdt1[i],
		     Form("E_{x} spectrum, gated on RDT1"),
		     Form("E_{x} (MeV)"),
		     Form("Counts / 70 keV"),
		     1);
    }

    cExSingleDetsRDT1->SaveAs((folderName + "/Ex_each_detector_rdt1.png").c_str());
    cExSingleDetsRDT1->SaveAs((folderName + "/Ex_each_detector_rdt1.pdf").c_str());

    TCanvas *cExSingleDetsRDT2 = new TCanvas("cExSingleDetsRDT2", "Ex for each detector, rdt2", 1200, 800);
    cExSingleDetsRDT2->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cExSingleDetsRDT2->cd(i + 1);
      Ex_single_rdt2[i]->Draw();
      StyleHistogram(Ex_single_rdt2[i],
		     Form("E_{x} spectrum, gated on RDT2"),
		     Form("E_{x} (MeV)"),
		     Form("Counts / 70 keV"),
		     1);
    }

    cExSingleDetsRDT2->SaveAs((folderName + "/Ex_each_detector_rdt2.png").c_str());
    cExSingleDetsRDT2->SaveAs((folderName + "/Ex_each_detector_rdt2.pdf").c_str());

    TCanvas *cExSingleDetsRDT3 = new TCanvas("cExSingleDetsRDT3", "Ex for each detector, rdt3", 1200, 800);
    cExSingleDetsRDT3->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cExSingleDetsRDT3->cd(i + 1);
      Ex_single_rdt3[i]->Draw();
      StyleHistogram(Ex_single_rdt3[i],
		     Form("E_{x} spectrum, gated on RDT3"),
		     Form("E_{x} (MeV)"),
		     Form("Counts / 70 keV"),
		     1);
    }

    cExSingleDetsRDT3->SaveAs((folderName + "/Ex_each_detector_rdt3.png").c_str());
    cExSingleDetsRDT3->SaveAs((folderName + "/Ex_each_detector_rdt3.pdf").c_str());
    
  }

}
