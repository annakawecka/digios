//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 25 17:40:16 2024 by ROOT version 6.33.01
// from TTree tree/tree
// found on file: trace_run055-066.root
//////////////////////////////////////////////////////////

#ifndef CreateTestFile_h
#define CreateTestFile_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <TBenchmark.h>
#include <TF1.h>
#include <string>
#include <fstream>
#include <TObjArray.h>
#include <TCutG.h>

#include "../Armory/AnalysisLibrary.h"

void readFitParameters(const TString &fileName, std::vector<std::vector<double>> &params) {
  std::ifstream inFile(fileName);
  if (!inFile.is_open()) {
    std::cerr << "Failed to open fit parameters file!" << std::endl;
    return;
  }

  std::string line;
  while (std::getline(inFile, line)) {
    std::istringstream iss(line);
    std::string temp;
    int detectorIndex;

    iss >> temp >> detectorIndex >> temp;

    std::vector<double> fitParams;
    double param;
    while (iss >> param) {
      fitParams.push_back(param);
    }

    if (detectorIndex >= params.size()) {
      params.resize(detectorIndex + 1);
    }
    params[detectorIndex] = fitParams;
  }
  inFile.close();
}


// Header file for the classes stored in the TTree if any.

class CreateTestFile {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventID;
   Int_t           run;
   Float_t         e[24];
   Float_t         xf[24];
   Float_t         xn[24];
   Float_t         ring[24];
   Float_t         x[24];
   Float_t         z[24];
   Int_t           detID;
   Int_t           hitID[24];
   Int_t           multiHit;
   Float_t         Ex;
   Float_t         thetaCM;
   Float_t         thetaLab;
   ULong64_t       e_t[24];
   Float_t         rdt[8];
   ULong64_t       rdt_t[8];
   Int_t           rdtID;
   Int_t           rdtdEMultiHit;
   Float_t         elum[32];
   ULong64_t       elum_t[32];
   Int_t           coin_t;
   Float_t         tcoin_t;
   Float_t         coinTimeUC;
   Float_t         coinTime;
   Float_t         te[24];
   Float_t         te_r[24];
   Float_t         te_t[24];
   Float_t         trdt[8];
   Float_t         trdt_t[8];
   Float_t         trdt_r[8];

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_run;   //!
   TBranch        *b_e;   //!
   TBranch        *b_xf;   //!
   TBranch        *b_xn;   //!
   TBranch        *b_x;   //!
   TBranch        *b_z;   //!
   TBranch        *b_detID;   //!
   TBranch        *b_hitID;   //!
   TBranch        *b_multiHit;   //!
   TBranch        *b_Ex;   //!
   TBranch        *b_thetaCM;   //!
   TBranch        *b_thetaLab;   //!
   TBranch        *b_e_t;   //!
   TBranch        *b_rdt;   //!
   TBranch        *b_rdt_t;   //!
   TBranch        *b_rdtID;   //!
   TBranch        *b_rdtdEMultiHit;   //!
   TBranch        *b_elum;   //!
   TBranch        *b_elum_t;   //!
   TBranch        *b_coin_t;   //!
   TBranch        *b_tcoin_t;   //!
   TBranch        *b_coinTimeUn;   //!
   TBranch        *b_coinTime;   //!
   TBranch        *b_te;   //!
   TBranch        *b_te_r;   //!
   TBranch        *b_te_t;   //!
   TBranch        *b_trdt;   //!
   TBranch        *b_trdt_t;   //!
   TBranch        *b_trdt_r;   //!

   CreateTestFile(TTree *tree=0);
   virtual ~CreateTestFile();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);

  //=============================== 
   TFile * saveFile;
   TTree * newTree;
   TString saveFileName;
   int totnumEntry; // of original root
  
   //tree  
   Int_t eventIDTemp;
   double eTemp;
   double xTemp; // unadjusted position, range (-1,1)
   double zTemp;
  double ExTemp;
   int detIDTemp;
   int hitIDTemp; // is e, xf, xn are all fired.
   int multiHitTemp; // multipicity of z
   int rdtdEMultiHitTemp; // multipicity of recoil-dE
   
   Float_t coinTimeTemp;
   
   TObjArray * fList; //!

  Int_t not_rej_rdt;
  Int_t rej_rdt;
  Int_t not_rej_coin;
  Int_t rej_coin;

  //========correction parameters
  int numDet;
  int iDet; // number of detector at different position
  int jDet; // number of detector at same position
  vector<double> pos;
  double Bfield;
  double a;
  double length;
  double firstPos;

  DetGeo detGeo;
   
  double xnCorr[30]; // xn correction for xn = xf
  double xfxneCorr[30][2]; //xf, xn correction for e = xf + xn
  double xScale[30]; //scale x to full range
  double cTCorr[30][9]; // coinTime correction
  std::vector<std::vector<double>> coinTimeFitParams;
  bool isCoinTimeLoaded;
  TF1 ** f7 ; //!
   
  //============ RDT cut
  bool isRDTCutExist;
  TCutG ** cut = NULL;
};

#endif

#ifdef CreateTestFile_cxx
CreateTestFile::CreateTestFile(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("trace_run055-066.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("trace_run055-066.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

CreateTestFile::~CreateTestFile()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CreateTestFile::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CreateTestFile::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CreateTestFile::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("xf", xf, &b_xf);
   fChain->SetBranchAddress("xn", xn, &b_xn);
   //fChain->SetBranchAddress("ring", ring, &b_xn);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("detID", &detID, &b_detID);
   fChain->SetBranchAddress("hitID", hitID, &b_hitID);
   fChain->SetBranchAddress("multiHit", &multiHit, &b_multiHit);
   fChain->SetBranchAddress("Ex", &Ex, &b_Ex);
   fChain->SetBranchAddress("thetaCM", &thetaCM, &b_thetaCM);
   fChain->SetBranchAddress("thetaLab", &thetaLab, &b_thetaLab);
   fChain->SetBranchAddress("e_t", e_t, &b_e_t);
   fChain->SetBranchAddress("rdt", rdt, &b_rdt);
   fChain->SetBranchAddress("rdt_t", rdt_t, &b_rdt_t);
   fChain->SetBranchAddress("rdtID", &rdtID, &b_rdtID);
   fChain->SetBranchAddress("rdtdEMultiHit", &rdtdEMultiHit, &b_rdtdEMultiHit);
   fChain->SetBranchAddress("elum", elum, &b_elum);
   fChain->SetBranchAddress("elum_t", elum_t, &b_elum_t);
   fChain->SetBranchAddress("coin_t", &coin_t, &b_coin_t);
   fChain->SetBranchAddress("tcoin_t", &tcoin_t, &b_tcoin_t);
   fChain->SetBranchAddress("coinTimeUC", &coinTimeUC, &b_coinTimeUn);
   fChain->SetBranchAddress("coinTime", &coinTime, &b_coinTime);
   fChain->SetBranchAddress("te", te, &b_te);
   fChain->SetBranchAddress("te_r", te_r, &b_te_r);
   fChain->SetBranchAddress("te_t", te_t, &b_te_t);
   fChain->SetBranchAddress("trdt", trdt, &b_trdt);
   fChain->SetBranchAddress("trdt_t", trdt_t, &b_trdt_t);
   fChain->SetBranchAddress("trdt_r", trdt_r, &b_trdt_r);

   //===================================================== loading parameter
   //========================================= detector Geometry
   string detGeoFileName = "detectorGeo.txt";
   printf("----- loading detector geometery : %s.", detGeoFileName.c_str());

   TMacro * haha = new TMacro();
   if( haha->ReadFile(detGeoFileName.c_str()) > 0 ) {

      detGeo = LoadDetectorGeo(haha);

      PrintDetGeo(detGeo);

      iDet = detGeo.nDet;
      jDet = detGeo.mDet;

      Bfield = detGeo.Bfield;
      a = detGeo.detPerpDist;
      length = detGeo.detLength;
      firstPos = detGeo.firstPos;
	pos = detGeo.detPos;
      
      printf("... done.\n");
   }else{
      printf("... fail\n");
      //Terminate();
   }

   numDet = iDet * jDet;

   //====================================== make tree
   saveFileName = "temp.root"; 
   saveFile = new TFile( saveFileName,"recreate");
   newTree =  new TTree("tree","tree");
   
   eventIDTemp = -1;
   
   newTree->Branch("eventID",&eventIDTemp,"eventID/I"); 
   
   newTree->Branch("e" ,  &eTemp, "eTemp/D");
   newTree->Branch("x" ,  &xTemp, "xTemp/D");
   newTree->Branch("z" ,  &zTemp, "zTemp/D");
   newTree->Branch("Ex" ,  &ExTemp, "ExTemp/D");
   newTree->Branch("detID", &detIDTemp, "detIDTemp/I");
   newTree->Branch("hitID", &hitIDTemp, "hitID/I");
   newTree->Branch("multiHit", &multiHitTemp, "multiHit/I");
   newTree->Branch("coinTime", &coinTimeTemp, "coinTime/F");

   //========================================= coinTime correction
   readFitParameters("coinTimeCorr_17O/coinTime_x_fit_parameters.txt", coinTimeFitParams);

   //====================================== load RDT cut
   TFile * fileCut = new TFile("rdtCuts_17O.root");   
   TObjArray * cutList = NULL;
   isRDTCutExist = false;
   if( fileCut->IsOpen() ){
      TObjArray * cutList = (TObjArray*) fileCut->FindObjectAny("cutList");
      
      if( cutList != NULL){
         isRDTCutExist = true;
         const int numCut = cutList->GetEntries();
         cut = new TCutG * [numCut];
         printf("=========== found %d cuts in %s \n", numCut, fileCut->GetName());
         for( int i = 0 ; i < numCut; i++){
            cut[i] = (TCutG* ) cutList->At(i);
            printf("cut name: %s , VarX: %s, VarY: %s\n", cut[i]->GetName(), cut[i]->GetVarX(), cut[i]->GetVarY()); 
         }
      }
   }
   
   Notify();
}

bool CreateTestFile::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void CreateTestFile::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CreateTestFile::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CreateTestFile_cxx
