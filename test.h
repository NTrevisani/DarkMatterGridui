////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////            DUMMY            ////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#ifndef MUONFR_NEWJETETDEF_H
#define MUONFR_NEWJETETDEF_H 1

#include "PAFAnalysis.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include "TCounterUI.h"
#include <TLorentzVector.h>
#include "Riostream.h"  
#include "PUWeight.h"

const Double_t ZMASS = 91.1876;

const UInt_t numberMetCuts = 5;
const UInt_t numberDYMVACuts = 5;

Double_t MetCut[] = {20, 25, 30, 45, 1000};

Double_t DYMVACut_0j[numberDYMVACuts] = {-0.9, -0.86, -0.6, 0.88, 1000};

Double_t DYMVACut_1j[numberDYMVACuts] = {-0.9, -0.86, -0.6, 0.84, 1000};

Bool_t runAtOviedo = false;
Bool_t runAtIfca   = !runAtOviedo;

Float_t SelectedChannel;

class test: public PAFAnalysis{

 public:
   test(TTree *tree=0);
   virtual ~test() {}

 protected:
   virtual void              Initialise();
   virtual void              InsideLoop();
   virtual void              SetDataMembersAtTermination();
   virtual void              Summary(); 

   TH1F* h_n_PV; 

   TTree *tree;

   // WW level histograms                                                         
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1WWLevel[4];
   TH1F* hPtLepton2WWLevel[4];
   TH1F* hPtDiLeptonWWLevel[4];
   TH1F* hMinvWWLevel[4];
   TH1F* hMtWWLevel[4];
   TH1F* hMt1WWLevel[4];
   TH1F* hMt2WWLevel[4];
   TH1F* hpfMetWWLevel[4];
   TH1F* hpminMetWWLevel[4];
   TH1F* hDeltaRLeptonsWWLevel[4];
   TH1F* hDeltaPhiLeptonsWWLevel[4];
   TH1F* hDPhiPtllJetWWLevel[4];
   TH1F* hPtWWWWLevel[4];
   TH1F* hMcWWLevel[4];
   TH1F* hTrkMetWWLevel[4];
   TH1F* hHtWWLevel[4];

   // TwoLeptons level histograms                                            
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1TwoLeptonsLevel;
   TH1F* hPtLepton2TwoLeptonsLevel ;
   TH1F* hPtDiLeptonTwoLeptonsLevel;
   TH1F* hMinvTwoLeptonsLevel;
   TH1F* hMtTwoLeptonsLevel;
   TH1F* hMt1TwoLeptonsLevel;
   TH1F* hMt2TwoLeptonsLevel;
   TH1F* hpfMetTwoLeptonsLevel;
   TH1F* hpminMetTwoLeptonsLevel;
   TH1F* hDeltaRLeptonsTwoLeptonsLevel;
   TH1F* hDeltaPhiLeptonsTwoLeptonsLevel;
   TH1F* hDPhiPtllJetTwoLeptonsLevel;
   TH1F* hDPhillMetTwoLeptonsLevel;
   TH1F* hMcTwoLeptonsLevel;
   TH1F* hPtWWTwoLeptonsLevel;
   TH1F* hTrkMetTwoLeptonsLevel;
   TH1F* hHtTwoLeptonsLevel;

 public:
 
   // My Declarations:
   // Define global variables

   PUWeight* fPUWeight;

   bool G_Debug_DefineAnalysisVariables;
 
   //* Histograms 
  
   // * Input parameters
   TString Signal; // Type of Signal
   int NEvents; // Total number of events in the sample before skim
   double Luminosity; // Total luminosity
   double XSection; // Process cross section
   bool IsDATA; // True if is Data, False in case MC
   int WhichRun; // 1 in case of RunI samples. 2 In case of RunII samples.
   int MultiplicateXS;//multiplicate Dark Matter XS by this factor
   TString TheSample; //path to the input files
   TString flavorChannel; //selected decay channel 
   int jetChannel; //number of jets in the event
   ClassDef(test,0);
};
#endif

