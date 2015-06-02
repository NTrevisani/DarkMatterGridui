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

   Float_t hPtLepton1;
   Float_t hPtLepton2;
   Float_t hPtDiLepton;
   Float_t hMinv;
   Float_t hMth;
   Float_t hMt1;
   Float_t hMt2;
   Float_t hNJets30;
   Float_t hpfMet;
   Float_t hmpMet;
   Float_t hDeltaRLeptons;
   Float_t hDeltaPhiLeptons;
   Float_t hDPhiPtllMet;
   Float_t hDPhiPtllJet;
   Float_t hPtWW;
   Float_t hMc;
   Float_t hTrkMet;
   Float_t hHt;
   Float_t htotalW;

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
