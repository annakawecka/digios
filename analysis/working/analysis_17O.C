#include "utilities.h"

#define N 5

TString rdtCutFile = "rdtCuts_17O.root";

TObjArray * cutList;
Bool_t isCutFileOpen;
int numCut;

TCutG* cutG;

bool rdtgate = false;
bool xgate = false;
bool cointimegate = false;

int n_coin = 0;
int n_coin_rdt = 0;

int recoil_n = 0;

//========correction parameters
DetGeo detGeo;

double eCorrections[24][2]; // e-correction (kinematics)

int numDet;
std::vector<double> pos;
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
bool fitting = false;

std::ifstream file;

const std::vector<double> peakPositions = {0.0, 1.982, 3.552, 3.630, 3.920, 5.255};

void analysis_17O(){
  //================================= coinTime fit parameters

  std::vector<std::vector<double>> fitParams;
  readFitParameters("coinTimeCorr_17O/coinTime_x_fit_parameters.txt", fitParams);

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
  chain->Add("trace_run055-066.root");

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

  TH1F* Ex_d0 = new TH1F("Ex_d0", "Ex, det 0, 6, 12, 18", 200, -2, 12);
  TH1F* Ex_d1 = new TH1F("Ex_d1", "Ex, det 1, 7, 13, 19", 200, -2, 12);
  TH1F* Ex_d2 = new TH1F("Ex_d2", "Ex, det 2, 8, 14, 20", 200, -2, 12);
  TH1F* Ex_d3 = new TH1F("Ex_d3", "Ex, det 3, 9, 15, 21", 200, -2, 12);
  TH1F* Ex_d4 = new TH1F("Ex_d4", "Ex, det 4, 10, 16, 22", 200, -2, 12);
  TH1F* Ex_d5 = new TH1F("Ex_d5", "Ex, det 5, 11, 17, 23", 200, -2, 12);

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

  // Loop over the events and apply the cut
  Long64_t nEntries = chain->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    rdtgate = false;
    xgate = false;
    cointimegate = false;

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

    // correcting e

    eCorr = e[detID]/eCorrections[detID][0] + eCorrections[detID][1];

    // calculating corrected Ex and ThetaCM

    double y = eCorr + mass;
    Z = alpha * gamma * beta * z[detID];
    H = TMath::Sqrt(TMath::Power(gamma * beta,2) * (y*y - mass * mass) ) ;

    if( TMath::Abs(Z) < H ) {
      // copied from Cali_e_trace
      ///using Newton's method to solve 0 ==  H * sin(phi) - G * tan(phi) - Z = f(phi)
      double tolerrence = 0.001;
      double phi = 0; ///initial phi = 0 -> ensure the solution has f'(phi) > 0
      double nPhi = 0; /// new phi

      int iter = 0;
      do{
	phi = nPhi;
	nPhi = phi - (H * TMath::Sin(phi) - G * TMath::Tan(phi) - Z) / (H * TMath::Cos(phi) - G /TMath::Power( TMath::Cos(phi), 2));
	iter ++;
	if( iter > 10 || TMath::Abs(nPhi) > TMath::PiOver2()) break;
      }while( TMath::Abs(phi - nPhi ) > tolerrence);
      phi = nPhi;

      /// check f'(phi) > 0
      double Df = H * TMath::Cos(phi) - G / TMath::Power( TMath::Cos(phi),2);
      if( Df > 0 && TMath::Abs(phi) < TMath::PiOver2()  ){
	double K = H * TMath::Sin(phi);
	double xx = TMath::ACos( mass / ( y * gamma - K));
	double k = mass * TMath::Tan(xx); /// momentum of particel b or B in CM frame
	double EB = TMath::Sqrt(mass*mass + Et*Et - 2*Et*TMath::Sqrt(k*k + mass * mass));
	ExCorr = EB - massB;

	double hahaha1 = gamma* TMath::Sqrt(mass * mass + k * k) - y;
	double hahaha2 = gamma* beta * k;
	thetaCM = TMath::ACos(hahaha1/hahaha2) * TMath::RadToDeg();

	double pt = k * TMath::Sin(thetaCM * TMath::DegToRad());
	double pp = gamma*beta*TMath::Sqrt(mass*mass + k*k) - gamma * k * TMath::Cos(thetaCM * TMath::DegToRad());

	thetaLab = TMath::ATan(pt/pp) * TMath::RadToDeg();
      }
    }

    //end of Ex and ThetaCM calculation

    if (coinTimeCorr > -20 && coinTimeCorr < 15) {
      cointimegate = true;

      coinTime_gatedEx->Fill(ExCorr);
    }

    correctedCoinTime[detID]->Fill(coinTimeCorr);
    n_coin++;

    if (xgate)
      correctedCoinTimeXgate[detID]->Fill(coinTimeCorr);

    if (rdtgate) {
      correctedCoinTimeRDTCoin[detID]->Fill(coinTimeCorr);
      recoilDEE[recoil_n/2]->Fill(rdt[recoil_n], rdt[recoil_n+1]);
      n_coin_rdt++;

      rdt_gatedEx->Fill(ExCorr);
    }

    if (rdtgate && cointimegate)
      rdt_coinTime_gatedEx->Fill(ExCorr);

    if (xgate && rdtgate) {
      correctedCoinTimeXgateRDTCoin[detID]->Fill(coinTimeCorr);

      if (cointimegate) {
	x_rdt_coinTime_gatedEx->Fill(ExCorr);

	if (detID % 6 == 0)
	  Ex_d0->Fill(ExCorr);
	if (detID % 6 == 1)
	  Ex_d1->Fill(ExCorr);
	if (detID % 6 == 2)
	  Ex_d2->Fill(ExCorr);
	if (detID % 6 == 3)
	  Ex_d3->Fill(ExCorr);
	if (detID % 6 == 4)
	  Ex_d4->Fill(ExCorr);
	if (detID % 6 == 5)
	  Ex_d5->Fill(ExCorr);
      }
    }

    Ex_nogates->Fill(ExCorr);

  } // end of loop over entries

  //================================= fitting 17O
  if (fitting) {
    TCanvas *cExdet = new TCanvas("cExdet", "Ex for different detector rings", 1000, 800);
    cExdet->Divide(3, 2);

    cExdet->cd(1);
    Ex_d0->Draw();

    cExdet->cd(2);
    Ex_d1->Draw();

    cExdet->cd(3);
    Ex_d2->Draw();

    cExdet->cd(4);
    Ex_d3->Draw();

    cExdet->cd(5);
    Ex_d4->Draw();

    cExdet->cd(6);
    Ex_d5->Draw();

    TF1 *f1 = new TF1("Fit",FitNPeaks,-1,7,2*N+3);

    f1->SetParameter(0,0.1);
    f1->SetParLimits(0,0.05,0.15);
    f1->SetParameter(1, 0);
    f1->SetParameter(2, 0); 

    f1->SetParName(0,"sigma");
    f1->SetParName(1,"p0");
    f1->SetParName(2,"p1");
    for(int i=0;i<N;i++){
      f1->SetParName(3+2*i,Form("Area%i",i+1));        
      f1->SetParName(4+2*i,Form("Cent%i",i+1));
    }

    // Peak 1
    f1->SetParameter(3,10);
    f1->SetParLimits(3,0,1000);
    f1->SetParameter(4,0);
    //f1->SetParLimits(4,-0.5,0.5);
    // Peak 2
    f1->SetParameter(5,10);
    f1->SetParLimits(5,0,1000);
    f1->SetParameter(6,2.0);
    //f1->SetParLimits(6,1.8,2.5);
    //f1->FixParameter(6,0.9376);
    // Peak 3
    f1->SetParameter(7,10);
    f1->SetParLimits(7,0,1000);
    f1->SetParameter(8,3.7);
    //f1->SetParLimits(8,2.9,3.9);
    //f1->FixParameter(8,1.0416);
    // Peak 4
    //f1->SetParameter(9,10);
    //f1->FixParameter(9,0);
    //f1->SetParLimits(9,0,1000);
    //f1->SetParameter(10,3.8);
    //f1->SetParLimits(10,3.1,3.9);
    //f1->FixParameter(10,3.9);
    // Peak 5
    f1->SetParameter(9,10);
    f1->SetParLimits(9,0,1000);
    f1->SetParameter(10,4.1);
    //f1->SetParLimits(10,3.9,4.3);
    // Peak 6
    f1->SetParameter(11,10);
    f1->SetParLimits(11,0,1000);
    f1->SetParameter(12,5.5);
    //f1->SetParLimits(12,5.3,6.0);
    //f1->FixParameter(12,1.1214);

  
    TCanvas *cExd5 = new TCanvas("cExd5", "Ex, det 5...", 800, 600);
    f1->SetNpx(1000);

    Ex_d5->Fit(f1, "RB0");

    Ex_d5->SetBinErrorOption(TH1::kPoisson);
  
    Ex_d5->Draw("e1");
    Ex_d5->SetMarkerStyle(20);
    Ex_d5->SetMarkerColor(kRed);
    Ex_d5->SetMarkerSize(1);
    f1->SetLineWidth(4);
    f1->Draw("same");

    printf("Detectors 5, 11, 17, 23 \n");

    TF1 *indpeak[N];
    for(int i=0;i<N;i++){
      indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
      indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
      if(i==0)
	indpeak[i]->SetParameter(1,f1->GetParameter(4));
      else
	indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
      indpeak[i]->SetParameter(2,f1->GetParameter(0));
      indpeak[i]->SetLineStyle(7);
      indpeak[i]->Draw("same");
      if (i==0)
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
      else
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
    }

    TCanvas *cExd4 = new TCanvas("cExd4", "Ex, det 4...", 800, 600);
    f1->SetNpx(1000);

    Ex_d4->Fit(f1, "RB0");

    Ex_d4->SetBinErrorOption(TH1::kPoisson);

    Ex_d4->Draw("e1");
    Ex_d4->SetMarkerStyle(20);
    Ex_d4->SetMarkerColor(kRed);
    Ex_d4->SetMarkerSize(1);
    f1->SetLineWidth(4);
    f1->Draw("same");

    printf("Detectors 4, 10, 16, 22 \n");

    for(int i=0;i<N;i++){
      indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
      indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
      if(i==0)
	indpeak[i]->SetParameter(1,f1->GetParameter(4));
      else
	indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
      indpeak[i]->SetParameter(2,f1->GetParameter(0));
      indpeak[i]->SetLineStyle(7);
      indpeak[i]->Draw("same");
      if (i==0)
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
      else
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
    }

    TCanvas *cExd3 = new TCanvas("cExd3", "Ex, det 3...", 800, 600);
    f1->SetNpx(1000);

    Ex_d3->Fit(f1, "RB0");

    Ex_d3->SetBinErrorOption(TH1::kPoisson);

    Ex_d3->Draw("e1");
    Ex_d3->SetMarkerStyle(20);
    Ex_d3->SetMarkerColor(kRed);
    Ex_d3->SetMarkerSize(1);
    f1->SetLineWidth(4);
    f1->Draw("same");

    printf("Detectors 3, 9, 15, 21 \n");

    for(int i=0;i<N;i++){
      indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
      indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
      if(i==0)
	indpeak[i]->SetParameter(1,f1->GetParameter(4));
      else
	indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
      indpeak[i]->SetParameter(2,f1->GetParameter(0));
      indpeak[i]->SetLineStyle(7);
      indpeak[i]->Draw("same");
      if (i==0)
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
      else
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
    }

    TCanvas *cExd2 = new TCanvas("cExd2", "Ex, det 2...", 800, 600);
    f1->SetNpx(1000);

    Ex_d2->Fit(f1, "RB0");

    Ex_d2->SetBinErrorOption(TH1::kPoisson);

    Ex_d2->Draw("e1");
    Ex_d2->SetMarkerStyle(20);
    Ex_d2->SetMarkerColor(kRed);
    Ex_d2->SetMarkerSize(1);
    f1->SetLineWidth(4);
    f1->Draw("same");

    printf("Detectors 2, 8, 14, 20 \n");

    for(int i=0;i<N;i++){
      indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
      indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
      if(i==0)
	indpeak[i]->SetParameter(1,f1->GetParameter(4));
      else
	indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
      indpeak[i]->SetParameter(2,f1->GetParameter(0));
      indpeak[i]->SetLineStyle(7);
      indpeak[i]->Draw("same");
      if (i==0)
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
      else
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
    }

    TCanvas *cExd1 = new TCanvas("cExd1", "Ex, det 1...", 800, 600);
    f1->SetNpx(1000);

    Ex_d1->Fit(f1, "RB0");

    Ex_d1->SetBinErrorOption(TH1::kPoisson);

    Ex_d1->Draw("e1");
    Ex_d1->SetMarkerStyle(20);
    Ex_d1->SetMarkerColor(kRed);
    Ex_d1->SetMarkerSize(1);
    f1->SetLineWidth(4);
    f1->Draw("same");

    printf("Detectors 1, 7, 13, 19 \n");

    for(int i=0;i<N;i++){
      indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
      indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
      if(i==0)
	indpeak[i]->SetParameter(1,f1->GetParameter(4));
      else
	indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
      indpeak[i]->SetParameter(2,f1->GetParameter(0));
      indpeak[i]->SetLineStyle(7);
      indpeak[i]->Draw("same");
      if (i==0)
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
      else
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
    }

    TCanvas *cExd0 = new TCanvas("cExd0", "Ex, det 0...", 800, 600);
    f1->SetNpx(1000);

    Ex_d0->Fit(f1, "RB0");

    Ex_d0->SetBinErrorOption(TH1::kPoisson);

    Ex_d0->Draw("e1");
    Ex_d0->SetMarkerStyle(20);
    Ex_d0->SetMarkerColor(kRed);
    Ex_d0->SetMarkerSize(1);
    f1->SetLineWidth(4);
    //f1->Draw("same");

    printf("Detectors 0, 6, 12, 18 \n");

    for(int i=0;i<N;i++){
      indpeak[i] = new TF1(Form("Peak_%i",i+1),Peak,f1->GetParameter(4+i*2)-0.5,f1->GetParameter(4+i*2)+0.5,3);
      indpeak[i]->SetParameter(0,f1->GetParameter(3+i*2));
      if(i==0)
	indpeak[i]->SetParameter(1,f1->GetParameter(4));
      else
	indpeak[i]->SetParameter(1,f1->GetParameter(4+i*2)+f1->GetParameter(4));
      indpeak[i]->SetParameter(2,f1->GetParameter(0));
      indpeak[i]->SetLineStyle(7);
      indpeak[i]->Draw("same");
      if (i==0)
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4));
      else
	printf("Peak %d: area: %f, centr: %f \n", i, f1->GetParameter(3+i*2), f1->GetParameter(4+i*2)+f1->GetParameter(4));
    }
  }

  //================================= plotting and saving histograms

  if (plothist)  {

    TCanvas *cAllDetectors = new TCanvas("cAllDetectors", "Corrected CoinTime for All Detectors", 1200, 800);
    cAllDetectors->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectors->cd(i + 1);
      correctedCoinTime[i]->Draw();
    }

    cAllDetectors->SaveAs("17O_analysis_ExCorr/Corrected_CoinTime_AllDetectors.png");

    TCanvas *cAllDetectorsXgate = new TCanvas("cAllDetectorsXgate", "Corrected CoinTime for All Detectors with x gate", 1200, 800);
    cAllDetectorsXgate->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectorsXgate->cd(i + 1);
      correctedCoinTimeXgate[i]->Draw();
    }

    cAllDetectorsXgate->SaveAs("17O_analysis_ExCorr/Corrected_CoinTime_AllDetectors_Xgate.png");

    TCanvas *cAllDetectorsRDT = new TCanvas("cAllDetectorsRDT", "Corrected CoinTime for All Detectors (RDT coin)", 1200, 800);
    cAllDetectorsRDT->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectorsRDT->cd(i + 1);
      correctedCoinTimeRDTCoin[i]->Draw();
    }

    cAllDetectorsRDT->SaveAs("17O_analysis_ExCorr/Corrected_CoinTime_AllDetectors_RDTCoin.png");

    TCanvas *cAllDetectorsXgateRDT = new TCanvas("cAllDetectorsXgateRDT", "Corrected CoinTime for All Detectors (RDT coin && x gate)", 1200, 800);
    cAllDetectorsXgateRDT->Divide(6, 4);

    for (int i = 0; i < numDet; ++i) {
      cAllDetectorsXgateRDT->cd(i + 1);
      correctedCoinTimeXgateRDTCoin[i]->Draw();
    }

    cAllDetectorsRDT->SaveAs("17O_analysis_ExCorr/Corrected_CoinTime_AllDetectors_XgateRDTCoin.png");

    std::cout << "coin: " << n_coin << " coin_rdt: " << n_coin_rdt << std::endl;

    TCanvas *cRDT = new TCanvas("cRDT", "RDT dE-E gated", 1200, 300);
    cRDT->Divide(4);

    for (int i = 0; i < 4; ++i) {
      cRDT->cd(i + 1);
      recoilDEE[i]->Draw();
    }

    cRDT->SaveAs("17O_analysis_ExCorr/RDT_DEE_gated.png");

    TH1F *combinedHist = (TH1F*)correctedCoinTimeXgateRDTCoin[0]->Clone("combinedHist");
    combinedHist->Reset();
    for (int i = 0; i < numDet; ++i) {
      combinedHist->Add(correctedCoinTimeXgateRDTCoin[i]);
    }

    TCanvas *cCombined = new TCanvas("cCombined", "Combined Corrected CoinTime with recoil && x gate", 800, 600);
    combinedHist->Draw();
    cCombined->SaveAs("17O_analysis_ExCorr/Combined_Corrected_CoinTime_XgateRDTCoin.png");

    TCanvas *cGatedEx = new TCanvas("cGatedEx", "Ex gated on x, recoils and coinTime", 800, 600);
    x_rdt_coinTime_gatedEx->Draw();
    cGatedEx->SaveAs("17O_analysis_ExCorr/Ex_x_recoil_coinTime_gated.png");

    TCanvas *cRDTCoinTimeGatedEx = new TCanvas("cRDTCoinTimeGatedEx", "Ex gated on recoils and coinTime", 800, 600);
    rdt_coinTime_gatedEx->Draw();
    cRDTCoinTimeGatedEx->SaveAs("17O_analysis_ExCorr/Ex_recoil_coinTime_gated.png");

    TCanvas *cCoinTimeGatedEx = new TCanvas("cCoinTimeGatedEx", "Ex gated on coinTime", 800, 600);
    coinTime_gatedEx->Draw();
    cCoinTimeGatedEx->SaveAs("17O_analysis_ExCorr/Ex_coinTime_gated.png");

    TCanvas *cRDTGatedEx = new TCanvas("cRDTGatedEx", "Ex gated on rdt", 800, 600);
    rdt_gatedEx->Draw();
    cRDTGatedEx->SaveAs("17O_analysis_ExCorr/Ex_rdt_gated.png");

    TCanvas *cExnogates = new TCanvas("cExnogates", "Ex, no gates", 800, 600);
    Ex_nogates->Draw();
    cExnogates->SaveAs("17O_analysis_ExCorr/Ex_nogates.png");

  }

}
