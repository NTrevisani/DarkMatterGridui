///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////                                                                                             /////////////
/////////////                                 mono Higgs ANALYSIS SELECTOR                                /////////////
/////////////                                                                                             /////////////
/////////////                                    Nicolò Trevisani (IFCA)                                  /////////////
/////////////                                          Jun 2016                                           /////////////
/////////////                                                                                             /////////////
/////////////                              -> Adjust to a 120 width window <-                             /////////////
/////////////                                                                                             /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "/gpfs/csic_projects/cms/sw/PAF/include/PAFChainItemSelector.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <vector>
#include "Riostream.h"  

const Double_t ZMASS = 91.1876;

class monoHiggsSelector: public PAFChainItemSelector{
  
 public:
  virtual ~monoHiggsSelector() {}
  
  virtual void                Initialise();
  virtual void                InsideLoop();
  virtual void                Summary();
  
 protected:
  // See description in implementation

  //PUWeight* fPUWeight;

  float Testing(int k);
  bool  IsTightLepton(int k, TString _MuonID_);
  float MuonIsolation(int k);
  float ElectronIsolation(int k);
  bool  IsIsolatedLepton(int k);

  bool G_Debug_DefineAnalysisVariables;

  // My Declarations:
  // Define data members

  int SelectedChannel;

  // VARIABLES FOR EACH EVENT (to be initialized after every event)
 
  std::vector<float> std_vector_lepton_pt;      
  std::vector<float> std_vector_jet_pt;         
  std::vector<float> std_vector_lepton_eta;      
  std::vector<float> std_vector_jet_eta;         
  std::vector<float> std_vector_lepton_phi;      
  std::vector<float> std_vector_jet_phi;         
  std::vector<float> std_vector_lepton_muSIP3D; 
  std::vector<float> std_vector_lepton_elSIP3D; 
  std::vector<float> std_vector_lepton_flavour;
  std::vector<float> std_vector_lepton_isTightMuon;
  std::vector<float> std_vector_lepton_chargedHadronIso;
  std::vector<float> std_vector_lepton_photonIso;
  std::vector<float> std_vector_lepton_neutralHadronIso;
  std::vector<float> std_vector_lepton_sumPUPt;
  std::vector<float> std_vector_electron_effectiveArea;
  std::vector<float> std_vector_lepton_BestTrackdxy;
  std::vector<float> std_vector_lepton_BestTrackdz;
  std::vector<float> std_vector_lepton_isMediumMuon;
  std::vector<float> std_vector_lepton_eleIdMedium;
  std::vector<float> std_vector_jet_csvv2ivf;
  std::vector<float> std_vector_jet_softMuPt;

  float jetphi2;
  float jetpt2;
  float jeteta1;
  float jeteta2;
  float phi1;
  float phi2;
  float jetRho;
  float puW;
  float effW;
  float triggW;
  float trigger;
  float dataset;
  float baseW;
  float dphilljetjet;
  float dphilmet1;
  float dphilmet2;
  float pfType1Met;
  float pt1;
  float pt2;
  float ptll;
  float mll;
  float mth;
  float drll;
  float dphill;
  float dphilljet;
  float trkMet;
  float ch1;
  float ch2;
  float channel;
  float njet;
  float nbjet;
  float nbjettche;
  float jetpt1;
  float jetphi1;
  float pfType1Metphi;
  float dphillmet;

  int   dphiv;
  int   nvtx;
  int   nextra;
  int   bveto_ip;
  int   bveto_mu;

  // VARIABLES FOR ALL EVENTS (to be initialized only once)

  TTree *mvaTree;
 
  // Histograms 
  TH1F *h_n_PV;
  
  // Counting histograms
  //----------------------------------------------------------------------------       
  
  TH1F* hsoftMuPt;
  TH1F* hjetPt;
  TH1F* hWTrigger;
  TH1F* hWeffTrigger;   
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
  
  // Control Region histograms
  //---------------------------------------------------------------------------- 
  
  TH1F* hpfMetCR;
  TH1F* htrkMetCR;
  TH1F* hmpMetCR;
  
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
  
  //Additional Variables
  
  Float_t fullpmet;
  Float_t trkpmet;
  Float_t mpmet;
  Float_t Ht;
  Float_t dphijet1met;
  Float_t ratioMet;
  Float_t ptWW;
  Float_t metvar;
  Float_t njetGen;
  Float_t Mt1;
  Float_t Mt2;
  Float_t Mc;
  Float_t dphijetjet;
  Float_t detajetjet;
  
  // My Declarations:OA
  // Define global variables
  
  // Input parameters
  TString                     _Signal;               // Type of Signal
  int                         _NEvents;              // Total number of events in the sample before skim
  float                       _Luminosity;           // Total luminosity
  float                       _XSection;             // Process cross section
  bool                        _IsDATA;               // True if is Data, False in case MC
  int                         _WhichRun;             // 1 in case of RunI samples. 2 In case of RunII samples.
  bool                        _Debug;                // True for verbose while debugging
  bool                        _Report;               // Count events and print final report
  TString                     _SameSign;             // Choose the type of events looking at the lepton's charge
  TString                     _FlavorChannel;        // Choose the type of events looking at the lepton's charge
  TString                     _outPath;              // Output folder
  TString                     _MuonID;               // ID requirement to consider a muon

  // Weights
  float                       _factN;                // Normalization factor



 public:  
 monoHiggsSelector() : 
  
  // VARIABLES FOR ALL EVENTS (to be initialized only once)

  h_n_PV(),
  
  // Counting histograms
  //----------------------------------------------------------------------------       
  
    hWExtraLepton(),   
    hWeffExtraLepton(),
    hWMetCut(),
    hWeffMetCut(),
    hWLowMinv(),
    hWeffLowMinv(),
    hWZVeto(),
    hWeffZVeto(),
    hWpMetCut(),
    hWeffpMetCut(),
    hWDeltaPhiJet(),
    hWeffDeltaPhiJet(),
    hWPtll(),
    hWeffPtll(),
    hWnJets(),
    hWeffnJets(),
    hWnBtaggedJets(),
    hWeffnBtaggedJets(),
    hWTopTagging(),
    hWeffTopTagging(),
    hWSoftMuVeto(),
    hWeffSoftMuVeto(),


   // TwoLeptons level histograms 
   //---------------------------------------------------------------------------- 
   
    hPtLepton1TwoLeptonsLevel(),
    hPtLepton2TwoLeptonsLevel (),
    hPtDiLeptonTwoLeptonsLevel(),
    hMinvTwoLeptonsLevel(),
    hMtTwoLeptonsLevel(),
    hMt1TwoLeptonsLevel(),
    hMt2TwoLeptonsLevel(),
    hpfMetTwoLeptonsLevel(),
    hpminMetTwoLeptonsLevel(),
    hDeltaRLeptonsTwoLeptonsLevel(),
    hDeltaPhiLeptonsTwoLeptonsLevel(),
    hDPhiPtllJetTwoLeptonsLevel(),
    hDPhillMetTwoLeptonsLevel(),
    hMcTwoLeptonsLevel(),
    hPtWWTwoLeptonsLevel(),
    hTrkMetTwoLeptonsLevel(),
    hHtTwoLeptonsLevel(),
    hdphijetjetTwoLeptonsLevel(),
    hdetajetjetTwoLeptonsLevel(),

   // WW level histograms                                            
   //---------------------------------------------------------------------------- 
   
    hPtLepton1WWLevel(),
    hPtLepton2WWLevel (),
    hPtDiLeptonWWLevel(),
    hMinvWWLevel(),
    hMtWWLevel(),
    hMt1WWLevel(),
    hMt2WWLevel(),
    hpfMetWWLevel(),
    hpminMetWWLevel(),
    hDeltaRLeptonsWWLevel(),
    hDeltaPhiLeptonsWWLevel(),
    hDPhiPtllJetWWLevel(),
    hDPhillMetWWLevel(),
    hMcWWLevel(),
    hPtWWWWLevel(),
    hTrkMetWWLevel(),
    hHtWWLevel(),
    hdphijetjetWWLevel(),
    hdetajetjetWWLevel(),

   // monoH level histograms
   //---------------------------------------------------------------------------- 

    hpfMetCR(),
    htrkMetCR(),
    hmpMetCR(),

   // monoH level histograms                                                         
   //---------------------------------------------------------------------------- 
   
    hPtLepton1monoHLevel(),
    hPtLepton2monoHLevel(),
    hPtDiLeptonmonoHLevel(),
    hMinvmonoHLevel(),
    hMtmonoHLevel(),
    hMt1monoHLevel(),
    hMt2monoHLevel(),
    hpfMetmonoHLevel(),
    hpminMetmonoHLevel(),
    hDeltaRLeptonsmonoHLevel(),
    hDeltaPhiLeptonsmonoHLevel(),
    hDPhiPtllJetmonoHLevel(),
    hDPhillMetmonoHLevel(),
    hPtWWmonoHLevel(),
    hMcmonoHLevel(),
    hTrkMetmonoHLevel(),
    hHtmonoHLevel(),
    hdphijetjetmonoHLevel(),
    hdetajetjetmonoHLevel(),

    _Signal(),
    _NEvents(),
    _Luminosity(),
    _XSection(),
    _IsDATA(),
    _WhichRun(),
    _Debug(),
    _Report(),
    _MuonID(),
    _factN()
      { }
  
  ClassDef(monoHiggsSelector,0);
};
