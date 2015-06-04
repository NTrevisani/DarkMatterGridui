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


   // Counting histograms                                                                  //----------------------------------------------------------------------------       

   TH1F* hWTrigger;
   TH1F* hWMetCut;
   TH1F* hWLowMinv;
   TH1F* hWZVeto;
   TH1F* hWpMetCut;
   TH1F* hWJetVeto;
   TH1F* hWnJetsBeforeBtag;
   TH1F* hWeffnJetsBeforeBtag;
   TH1F* hWnJets;
   TH1F* hWeffnJets;
   TH1F* hWnBtaggedJets;
   TH1F* hWeffnBtaggedJets;
   TH1F* hWnJetsBveto;
   TH1F* hWeffnJetsBveto;
   TH1F* hWnJetsBvetoAfterHt;
   TH1F* hWeffnJetsBvetoAfterHt;
   TH1F* hWmth;
   TH1F* hWeffmth;
   TH1F* hWDeltaPhiJet;
   TH1F* hWSoftMuVeto;
   TH1F* hWExtraLepton;
   TH1F* hWPtll;
   TH1F* hWTopTagging;
   TH1F* hWdphill;
   TH1F* hWtrkMet;
   TH1F* hWefftrkMet;
   TH1F* hWdrll;
   TH1F* hWeffdrll;

   TH1F* hWeffTrigger;
   TH1F* hWeffMetCut;
   TH1F* hWeffLowMinv;
   TH1F* hWeffZVeto;
   TH1F* hWeffpMetCut;
   TH1F* hWeffJetVeto;
   TH1F* hWeffDeltaPhiJet;
   TH1F* hWeffSoftMuVeto;
   TH1F* hWeffExtraLepton;
   TH1F* hWeffPtll;
   TH1F* hWeffTopTagging;
   TH1F* hWeffdphill;

   // WW level histograms                                                         
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1WWLevel[4];
   TH1F* hPtLepton2WWLevel[4];
   TH1F* hPtDiLeptonWWLevel[4];
   TH1F* hMinvWWLevel[4];
   TH1F* hMtWWLevel[4];
   TH1F* hNJets30WWLevel[4];
   TH1F* hpfMetWWLevel[4];
   TH1F* hppfMetWWLevel[4];
   TH1F* hchMetWWLevel[4];
   TH1F* hpchMetWWLevel[4];
   TH1F* hpminMetWWLevel[4];
   TH1F* hDeltaRLeptonsWWLevel[4];
   TH1F* hDeltaPhiLeptonsWWLevel[4];
   TH1F* hDPhiPtllJetWWLevel[4];
   TH1F* hDataEvents[4];
   TH1F* hBackgroundEvents[4];
   TH1F* hEff[4];
   TH1F* hPtWW[4];
   TH1F* hMcWWLevel[4];
   TH1F* hTrkMetWWLevel[4];

   TH1F* hPtLepton1WWLevelNoHt[4];
   TH1F* hPtLepton2WWLevelNoHt[4];
   TH1F* hPtDiLeptonWWLevelNoHt[4];
   TH1F* hMinvWWLevelNoHt[4];
   TH1F* hMtWWLevelNoHt[4];
   TH1F* hNJets30WWLevelNoHt[4];
   TH1F* hpfMetWWLevelNoHt[4];
   TH1F* hppfMetWWLevelNoHt[4];
   TH1F* hchMetWWLevelNoHt[4];
   TH1F* hpchMetWWLevelNoHt[4];
   TH1F* hpminMetWWLevelNoHt[4];
   TH1F* hDeltaRLeptonsWWLevelNoHt[4];
   TH1F* hDeltaPhiLeptonsWWLevelNoHt[4];
   TH1F* hDPhiPtllJetWWLevelNoHt[4];
   TH1F* hDataEventsNoHt[4];
   TH1F* hBackgroundEventsNoHt[4];
   TH1F* hEffNoHt[4];
   TH1F* hPtWWNoHt[4];
   TH1F* hMcWWLevelNoHt[4];

   TH1F* hPtLepton1WWLevelHtPlus[4];
   TH1F* hPtLepton2WWLevelHtPlus[4];
   TH1F* hPtDiLeptonWWLevelHtPlus[4];
   TH1F* hMinvWWLevelHtPlus[4];
   TH1F* hMtWWLevelHtPlus[4];
   TH1F* hNJets30WWLevelHtPlus[4];
   TH1F* hpfMetWWLevelHtPlus[4];
   TH1F* hppfMetWWLevelHtPlus[4];
   TH1F* hchMetWWLevelHtPlus[4];
   TH1F* hpchMetWWLevelHtPlus[4];
   TH1F* hpminMetWWLevelHtPlus[4];
   TH1F* hDeltaRLeptonsWWLevelHtPlus[4];
   TH1F* hDeltaPhiLeptonsWWLevelHtPlus[4];
   TH1F* hDPhiPtllJetWWLevelHtPlus[4];
   TH1F* hDataEventsHtPlus[4];
   TH1F* hBackgroundEventsHtPlus[4];
   TH1F* hEffHtPlus[4];
   TH1F* hPtWWHtPlus[4];
   TH1F* hMcWWLevelHtPlus[4];

   TH1F* hHt[4];
   TH1F* hHtAfter[4];
   TH1F* hPtWWNoBtagging[4];

   TH1F* h_WWLevel_TightFailEvents;
   TH1F* h_WWLevel_TightTightEvents;
   TH1F* h_WWLevel_TightLooseEvents;


   // TwoLeptons level histograms                                            
   //---------------------------------------------------------------------------- 
   
   TH1F* hPtLepton1TwoLeptonsLevel;
   TH1F* hPtLepton2TwoLeptonsLevel ;
   TH1F* hPtDiLeptonTwoLeptonsLevel;
   TH1F* hMinvTwoLeptonsLevel;
   TH1F* hMtTwoLeptonsLevel;
   TH1F* hMt1TwoLeptonsLevel;
   TH1F* hMt2TwoLeptonsLevel;
   TH1F* hNJets30TwoLeptonsLevel;
   TH1F* hpfMetTwoLeptonsLevel;
   TH1F* hppfMetTwoLeptonsLevel;
   TH1F* hchMetTwoLeptonsLevel;
   TH1F* hpchMetTwoLeptonsLevel;
   TH1F* hpminMetTwoLeptonsLevel;
   TH1F* hDeltaRLeptonsTwoLeptonsLevel;
   TH1F* hDeltaPhiLeptonsTwoLeptonsLevel;
   TH1F* hDPhiPtllJetTwoLeptonsLevel;
   TH1F* hDPhillMetTwoLeptonsLevel;
   TH1F* hMcTwoLeptonsLevel;
   TH1F* hPtWWTwoLeptonsLevel;
   TH1F* hTrkMetTwoLeptonsLevel;

   TH1F* h_TwoLeptons_TightFailEvents;
   TH1F* h_TwoLeptons_TightTightEvents;
   TH1F* h_TwoLeptons_TightLooseEvents;

   //Isolation Plots
   //-----------------------------------------------------------------------------

   TH1F* hIsoMu;
   TH1F* hIsoEl;

   TH1F* hchHadronMu;
   TH1F* hchHadronEl;
   
   TH1F* hneuHadronMu;
   TH1F* hneuHadronEl;
   
   TH1F* hphotonMu;
   TH1F* hphotonEl;
   
   TH1F* hPUMu;
   TH1F* hPUEl;


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
int WhichRun; // 1 in case of RunI samples. 2 In case of RunII samples.;
 TString TheSample; //path to the input files
 TString flavorChannel; //selected decay channel 
 int jetChannel; //number of jets in the event
   ClassDef(test,0);
};
#endif

