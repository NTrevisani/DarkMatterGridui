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

   TTree *mvaTree;
   TTree *tree;

   Float_t Mt1;
   Float_t Mt2;
   Float_t Mc;
   Float_t ptWW;
   Float_t Ht;
   Float_t mpmet;
   Float_t dphijetjet;
   Float_t detajetjet;

   // Counting histograms 
   //---------------------------------------------------------------------------- 

   TH1F* hWExtraLepton;   
   TH1F* hWeffExtraLepton;
   TH1F* hWMetCut;
   TH1F* hWeffMetCut;
   TH1F* hWLowMinv;
   TH1F* hWeffLowMinv;
   TH1F* hWZVeto;
   TH1F* hWeffZVeto;
   TH1F* hWpMetCut;
   TH1F* hWeffpMetCut;
   TH1F* hWDeltaPhiJet;
   TH1F* hWeffDeltaPhiJet;
   TH1F* hWPtll;
   TH1F* hWeffPtll;
   TH1F* hWnJets;
   TH1F* hWeffnJets;
   TH1F* hWnBtaggedJets;
   TH1F* hWeffnBtaggedJets;
   TH1F* hWTopTagging;
   TH1F* hWeffTopTagging;
   TH1F* hWSoftMuVeto;
   TH1F* hWeffSoftMuVeto;


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
   TH1F* hdphijetjetTwoLeptonsLevel;
   TH1F* hdetajetjetTwoLeptonsLevel;

   // WW level histograms                                            
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1WWLevel;
   TH1F* hPtLepton2WWLevel ;
   TH1F* hPtDiLeptonWWLevel;
   TH1F* hMinvWWLevel;
   TH1F* hMtWWLevel;
   TH1F* hMt1WWLevel;
   TH1F* hMt2WWLevel;
   TH1F* hpfMetWWLevel;
   TH1F* hpminMetWWLevel;
   TH1F* hDeltaRLeptonsWWLevel;
   TH1F* hDeltaPhiLeptonsWWLevel;
   TH1F* hDPhiPtllJetWWLevel;
   TH1F* hDPhillMetWWLevel;
   TH1F* hMcWWLevel;
   TH1F* hPtWWWWLevel;
   TH1F* hTrkMetWWLevel;
   TH1F* hHtWWLevel;
   TH1F* hdphijetjetWWLevel;
   TH1F* hdetajetjetWWLevel;

   // monoH level histograms                                                         
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1monoHLevel[4];
   TH1F* hPtLepton2monoHLevel[4];
   TH1F* hPtDiLeptonmonoHLevel[4];
   TH1F* hMinvmonoHLevel[4];
   TH1F* hMtmonoHLevel[4];
   TH1F* hMt1monoHLevel[4];
   TH1F* hMt2monoHLevel[4];
   TH1F* hpfMetmonoHLevel[4];
   TH1F* hpminMetmonoHLevel[4];
   TH1F* hDeltaRLeptonsmonoHLevel[4];
   TH1F* hDeltaPhiLeptonsmonoHLevel[4];
   TH1F* hDPhiPtllJetmonoHLevel[4];
   TH1F* hDPhillMetmonoHLevel[4];
   TH1F* hPtWWmonoHLevel[4];
   TH1F* hMcmonoHLevel[4];
   TH1F* hTrkMetmonoHLevel[4];
   TH1F* hHtmonoHLevel[4];
   TH1F* hdphijetjetmonoHLevel[4];
   TH1F* hdetajetjetmonoHLevel[4];

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
   TString TheSample; //path to the input files
   TString flavorChannel; //selected decay channel 
   int jetChannel; //number of jets in the event
   ClassDef(test,0);
};
#endif

