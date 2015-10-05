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

#include "monoHiggsSelector.h"

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include <vector>
#include "TROOT.h"
#include <iostream>
#include "TSystem.h"

#include "TDatabasePDG.h"

ClassImp(monoHiggsSelector)

float dist(float phi1, float eta1, float phi2, float eta2){
  float dphi = fabs(phi1 - phi2);
  float deta = fabs(eta1 - eta2);
  if (dphi > TMath::Pi())
    dphi = 2*TMath::Pi() - dphi;
  return sqrt(dphi*dphi + deta*deta);
}


// Initialise input parameters and data members for all events
void monoHiggsSelector::Initialise() {

  _Signal          = GetParam<TString>("Signal");
  _IsDATA          = GetParam<bool>("IsDATA");
  _NEvents         = GetParam<int>("NEvents");
  _Luminosity      = GetParam<float>("Luminosity");
  _XSection        = GetParam<float>("XSection");
  _WhichRun        = GetParam<int>("WhichRun"); 
  _Debug           = GetParam<bool>("Debug");
  _Report          = GetParam<bool>("Report");
  _SameSign        = GetParam<TString>("SameSign"); 
  _FlavorChannel   = GetParam<TString>("FlavorChannel");
  _outPath         = GetParam<TString>("OutPath");
  _MuonID          = GetParam<TString>("MuonID");

  gSystem -> Exec("mkdir " + _outPath);

  //Define weights
  _factN = 1.;
  //  if (!_IsDATA && _XSection > 0) _factN = _XSection * _Luminosity / _NEvents;

  //For counting
  if (_Report) {

  }

  SelectedChannel = -999;
  
  if      (_FlavorChannel == "MuMu") SelectedChannel =  0;
  else if (_FlavorChannel == "EE"  ) SelectedChannel =  1;
  else if (_FlavorChannel == "EMu" ) SelectedChannel =  2;
  else if (_FlavorChannel == "MuE" ) SelectedChannel =  3;
  else if (_FlavorChannel == "OF"  ) SelectedChannel =  4;
  else if (_FlavorChannel == "SF"  ) SelectedChannel =  5;
  
  else if (_FlavorChannel == "All" ) SelectedChannel = -1;
 

  //------------------------------------------------------------------------------
  // Create tree and branches
  //------------------------------------------------------------------------------
  mvaTree = CreateTree("nt","nt");

  //Variables to Be Studied
  mvaTree->Branch("pt1",&pt1);       
  mvaTree->Branch("pt2",&pt2);       
  mvaTree->Branch("ptll",&ptll);      
  mvaTree->Branch("mll",&mll);       
  mvaTree->Branch("mth",&mth);       
  mvaTree->Branch("Mt1",&Mt1);       
  mvaTree->Branch("Mt2",&Mt2);       
  mvaTree->Branch("pfType1Met",&pfType1Met);
  mvaTree->Branch("mpmet",&mpmet);     
  mvaTree->Branch("drll",&drll);      
  mvaTree->Branch("dphill",&dphill);    
  mvaTree->Branch("dphilljet",&dphilljet); 
  mvaTree->Branch("dphillmet",&dphillmet); 
  mvaTree->Branch("Mc",&Mc);        
  mvaTree->Branch("ptWW",&ptWW);      
  mvaTree->Branch("trkMet",&trkMet);    
  mvaTree->Branch("Ht",&Ht);        
  mvaTree->Branch("baseW",&baseW);        

  //------------------------------------------------------------------------------
  // Create histos
  //------------------------------------------------------------------------------
  
  h_n_PV = CreateH1F("h_n_PV","h_n_PV",50,0.,10.);
  
  // Counting histograms    
  //---------------------------------------------------------------------------- 
  
  hWExtraLepton     = CreateH1F("hWExtraLepton",     "", 10, 0, 10);    
  hWeffExtraLepton  = CreateH1F("hWeffExtraLepton",  "", 10, 0, 10);
  hWMetCut          = CreateH1F("hWMetCut",          "", 10, 0, 10);
  hWeffMetCut       = CreateH1F("hWeffMetCut",       "", 10, 0, 10);
  hWLowMinv         = CreateH1F("hWLowMinv",         "", 10, 0, 10);
  hWeffLowMinv      = CreateH1F("hWeffLowMinv",      "", 10, 0, 10);
  hWZVeto           = CreateH1F("hWZVeto",           "", 10, 0, 10);
  hWeffZVeto        = CreateH1F("hWeffZVeto",        "", 10, 0, 10);
  hWpMetCut         = CreateH1F("hWpMetCut",         "", 10, 0, 10);
  hWeffpMetCut      = CreateH1F("hWeffpMetCut",      "", 10, 0, 10);
  hWDeltaPhiJet     = CreateH1F("hWDeltaPhiJet",     "", 10, 0, 10);
  hWeffDeltaPhiJet  = CreateH1F("hWeffDeltaPhiJet",  "", 10, 0, 10);
  hWPtll            = CreateH1F("hWPtll",            "", 10, 0, 10);
  hWeffPtll         = CreateH1F("hWeffPtll",         "", 10, 0, 10);
  hWnJets           = CreateH1F("hWnJets",           "", 10, 0, 10);
  hWeffnJets        = CreateH1F("hWeffnJets",        "", 10, 0, 10);
  hWnBtaggedJets    = CreateH1F("hWnBtaggedJets",    "", 10, 0, 10);
  hWeffnBtaggedJets = CreateH1F("hWeffnBtaggedJets", "", 10, 0, 10);
  hWTopTagging      = CreateH1F("hWTopTagging",      "", 10, 0, 10);
  hWeffTopTagging   = CreateH1F("hWeffTopTagging",   "", 10, 0, 10);
  hWSoftMuVeto      = CreateH1F("hWSoftMuVeto",      "", 10, 0, 10);
  hWeffSoftMuVeto   = CreateH1F("hWeffSoftMuVeto",   "", 10, 0, 10);
  
  // TwoLeptons level histograms   
  //----------------------------------------------------------------------------
  
  hPtLepton1TwoLeptonsLevel       = CreateH1F("hPtLepton1TwoLeptonsLevel",       "", 3000, 0, 3000);
  hPtLepton2TwoLeptonsLevel       = CreateH1F("hPtLepton2TwoLeptonsLevel",       "", 3000, 0, 3000);
  hPtDiLeptonTwoLeptonsLevel      = CreateH1F("hPtDiLeptonTwoLeptonsLevel",      "", 3000, 0, 3000);
  hMinvTwoLeptonsLevel            = CreateH1F("hMinvTwoLeptonsLevel",            "", 3000, 0, 3000);
  hMtTwoLeptonsLevel              = CreateH1F("hMtTwoLeptonsLevel",              "", 3000, 0, 3000);
  hMt1TwoLeptonsLevel             = CreateH1F("hMt1TwoLeptonsLevel",             "", 3000, 0, 3000);
  hMt2TwoLeptonsLevel             = CreateH1F("hMt2TwoLeptonsLevel",             "", 3000, 0, 3000);
  hpfMetTwoLeptonsLevel           = CreateH1F("hpfMetTwoLeptonsLevel",           "", 3000, 0, 3000);
  hpminMetTwoLeptonsLevel         = CreateH1F("hpminMetTwoLeptonsLevel",         "", 3000, 0, 3000);
  hDeltaRLeptonsTwoLeptonsLevel   = CreateH1F("hDeltaRLeptonsTwoLeptonsLevel",   "",   50, 0,    5);
  hDeltaPhiLeptonsTwoLeptonsLevel = CreateH1F("hDeltaPhiLeptonsTwoLeptonsLevel", "",   32, 0,  3.2);
  hDPhiPtllJetTwoLeptonsLevel     = CreateH1F("hDPhiPtllJetTwoLeptonsLevel",     "",   32, 0,  3.2);
  hDPhillMetTwoLeptonsLevel       = CreateH1F("hDPhillMetTwoLeptonsLevel",       "",   32, 0,  3.2);
  hMcTwoLeptonsLevel              = CreateH1F("hMcTwoLeptonsLevel",              "", 3000, 0, 3000);  
  hPtWWTwoLeptonsLevel            = CreateH1F("hPtWWTwoLeptonsLevel",            "", 3000, 0, 3000);
  hTrkMetTwoLeptonsLevel          = CreateH1F("hTrkMetTwoLeptonsLevel",          "", 3000, 0, 3000);
  hHtTwoLeptonsLevel              = CreateH1F("hHtTwoLeptonsLevel",              "", 3000, 0, 3000);
  hdphijetjetTwoLeptonsLevel      = CreateH1F("hdphijetjetTwoLeptonsLevel",      "",   32, 0,  3.2);
  hdetajetjetTwoLeptonsLevel      = CreateH1F("hdetajetjetTwoLeptonsLevel",      "",   50, 0,  5.0);


  // WW level histograms     
  //----------------------------------------------------------------------------  
  
  hPtLepton1WWLevel       = CreateH1F("hPtLepton1WWLevel",       "", 3000, 0, 3000);
  hPtLepton2WWLevel       = CreateH1F("hPtLepton2WWLevel",       "", 3000, 0, 3000);
  hPtDiLeptonWWLevel      = CreateH1F("hPtDiLeptonWWLevel",      "", 3000, 0, 3000);
  hMinvWWLevel            = CreateH1F("hMinvWWLevel",            "", 3000, 0, 3000);
  hMtWWLevel              = CreateH1F("hMtWWLevel",              "", 3000, 0, 3000);
  hMt1WWLevel             = CreateH1F("hMt1WWLevel",             "", 3000, 0, 3000);
  hMt2WWLevel             = CreateH1F("hMt2WWLevel",             "", 3000, 0, 3000);
  hpfMetWWLevel           = CreateH1F("hpfMetWWLevel",           "", 3000, 0, 3000);
  hpminMetWWLevel         = CreateH1F("hpminMetWWLevel",         "", 3000, 0, 3000);
  hDeltaRLeptonsWWLevel   = CreateH1F("hDeltaRLeptonsWWLevel",   "",   50, 0,    5);
  hDeltaPhiLeptonsWWLevel = CreateH1F("hDeltaPhiLeptonsWWLevel", "",   32, 0,  3.2);
  hDPhiPtllJetWWLevel     = CreateH1F("hDPhiPtllJetWWLevel",     "",   32, 0,  3.2);
  hDPhillMetWWLevel       = CreateH1F("hDPhillMetWWLevel",       "",   32, 0,  3.2);
  hMcWWLevel              = CreateH1F("hMcWWLevel",              "", 3000, 0, 3000);  
  hPtWWWWLevel            = CreateH1F("hPtWWWWLevel",            "", 3000, 0, 3000);
  hTrkMetWWLevel          = CreateH1F("hTrkMetWWLevel",          "", 3000, 0, 3000);
  hHtWWLevel              = CreateH1F("hHtWWLevel",              "", 3000, 0, 3000);
  hdphijetjetWWLevel      = CreateH1F("hdphijetjetWWLevel",      "",   32, 0,  3.2);
  hdetajetjetWWLevel      = CreateH1F("hdetajetjetWWLevel",      "",   50, 0,  5.0);
    
  //Control Region Plots
  //----------------------------------------------------------------------------                

  hpfMetCR  = CreateH1F("hpfMetCR",           "", 3000, 0, 3000);
  htrkMetCR = CreateH1F("htrkMetCR",          "", 3000, 0, 3000);
  hmpMetCR  = CreateH1F("hmpMetCR",           "", 3000, 0, 3000);


  // monoH level histograms     
  //----------------------------------------------------------------------------  

  for (Int_t nC=0; nC<4; nC++) {

    hPtLepton1monoHLevel[nC]       = CreateH1F(Form("hPtLepton1monoHLevel%.1i", nC),       "", 3000, 0, 3000);
    hPtLepton2monoHLevel[nC]       = CreateH1F(Form("hPtLepton2monoHLevel%.1i", nC),       "", 3000, 0, 3000);
    hPtDiLeptonmonoHLevel[nC]      = CreateH1F(Form("hPtDiLeptonmonoHLevel%.1i", nC),      "", 3000, 0, 3000);
    hMinvmonoHLevel[nC]            = CreateH1F(Form("hMinvmonoHLevel%.1i", nC),            "", 3000, 0, 3000);
    hMtmonoHLevel[nC]              = CreateH1F(Form("hMtmonoHLevel%.1i", nC),              "", 3000, 0, 3000);
    hMt1monoHLevel[nC]             = CreateH1F(Form("hMt1monoHLevel%.1i", nC),             "", 3000, 0, 3000);
    hMt2monoHLevel[nC]             = CreateH1F(Form("hMt2monoHLevel%.1i", nC),             "", 3000, 0, 3000);
    hpfMetmonoHLevel[nC]           = CreateH1F(Form("hpfMetmonoHLevel%.1i", nC),           "", 3000, 0, 3000);
    hpminMetmonoHLevel[nC]         = CreateH1F(Form("hpminMetmonoHLevel%.1i", nC),         "", 3000, 0, 3000);
    hDeltaRLeptonsmonoHLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsmonoHLevel%.1i", nC),   "",   50, 0,    5);
    hDeltaPhiLeptonsmonoHLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsmonoHLevel%.1i", nC), "",   32, 0,  3.2);
    hDPhiPtllJetmonoHLevel[nC]     = CreateH1F(Form("hDPhiPtllJetmonoHLevel%.1i", nC),     "",   32, 0,  3.2);
    hDPhillMetmonoHLevel[nC]       = CreateH1F(Form("hDPhillMetmonoHLevel%.1i", nC),       "",   32, 0,  3.2);
    hPtWWmonoHLevel[nC]            = CreateH1F(Form("hPtWWmonoHLevel%.1i", nC),            "", 3000, 0, 3000);
    hMcmonoHLevel[nC]              = CreateH1F(Form("hMcmonoHLevel%.1i", nC),              "", 3000, 0, 3000);
    hTrkMetmonoHLevel[nC]          = CreateH1F(Form("hTrkMetmonoHLevel%.1i", nC),          "", 3000, 0, 3000);
    hHtmonoHLevel[nC]              = CreateH1F(Form("hHtmonoHLevel%.1i", nC),              "", 3000, 0, 3000);
    hdphijetjetmonoHLevel[nC]      = CreateH1F(Form("hdphijetjetmonoHLevel%.1i", nC),      "",   32, 0,  3.2);
    hdetajetjetmonoHLevel[nC]      = CreateH1F(Form("hdetajetjetmonoHLevel%.1i", nC),      "",   50, 0,  5.0);
  }
}

void monoHiggsSelector::InsideLoop() {

  //Assigning values to vectors
  //----------------------------------------------------------------------------
  Assign("std_vector_lepton_pt",std_vector_lepton_pt);
  Assign("std_vector_jet_pt",std_vector_jet_pt);    
  Assign("std_vector_lepton_muSIP3D",std_vector_lepton_muSIP3D);
  Assign("std_vector_lepton_elSIP3D",std_vector_lepton_elSIP3D);
  Assign("std_vector_lepton_flavour",std_vector_lepton_flavour);     
  Assign("std_vector_lepton_isTightMuon",std_vector_lepton_isTightMuon);
  Assign("std_vector_lepton_eleIdMedium",std_vector_lepton_eleIdMedium);
  Assign("std_vector_lepton_chargedHadronIso",std_vector_lepton_chargedHadronIso);
  Assign("std_vector_lepton_photonIso",std_vector_lepton_photonIso);       
  Assign("std_vector_lepton_neutralHadronIso",std_vector_lepton_neutralHadronIso);
  Assign("std_vector_lepton_sumPUPt",std_vector_lepton_sumPUPt); 
  Assign("std_vector_lepton_chargedHadronIso",std_vector_lepton_chargedHadronIso);
  Assign("std_vector_lepton_photonIso",std_vector_lepton_photonIso);  
  Assign("std_vector_lepton_neutralHadronIso",std_vector_lepton_neutralHadronIso);
  Assign("std_vector_electron_effectiveArea",std_vector_electron_effectiveArea); 
  Assign("std_vector_lepton_BestTrackdxy",std_vector_lepton_BestTrackdxy);
  Assign("std_vector_lepton_BestTrackdz",std_vector_lepton_BestTrackdz);  
  Assign("std_vector_lepton_phi",std_vector_lepton_phi);
  Assign("std_vector_lepton_eta",std_vector_lepton_eta);
  Assign("std_vector_jet_eta",std_vector_jet_eta);
  Assign("std_vector_jet_phi",std_vector_jet_phi);  
  Assign("std_vector_lepton_isMediumMuon",std_vector_lepton_isMediumMuon);
  Assign("std_vector_jet_csvv2ivf",std_vector_jet_csvv2ivf);
  Assign("std_vector_jet_softMuPt",std_vector_jet_softMuPt);

  //Assigning values to float variables
  //----------------------------------------------------------------------------
  
  Assign("phi1",phi1);
  Assign("phi2",phi2);
  Assign("puW",puW);  
  Assign("effW",effW);
  Assign("triggW",triggW);
  Assign("trigger",trigger);
  Assign("dataset",dataset);
  Assign("baseW",baseW);
  Assign("dphilljetjet",dphilljetjet);
  Assign("dphilmet1",dphilmet1);
  Assign("dphilmet2",dphilmet2);
  Assign("pfType1Met",pfType1Met);  
  Assign("pt1",pt1);
  Assign("pt2",pt2);
  Assign("ptll",ptll);
  Assign("mll",mll);
  Assign("mth",mth);
  Assign("drll",drll);
  Assign("dphill",dphill);
  Assign("dphilljet",dphilljet);
  Assign("trkMet",trkMet);
  Assign("ch1",ch1);
  Assign("ch2",ch2);
  Assign("jetRho",jetRho);
  Assign("channel",channel);
  Assign("njet",njet);
  Assign("nbjet",nbjet);
  Assign("nbjettche",nbjettche);
  Assign("jetpt1",jetpt1);
  Assign("jetphi1",jetphi1);
  Assign("pfType1Metphi",pfType1Metphi);
  Assign("dphillmet",dphillmet);

  //Assigning values to integer variables
  //----------------------------------------------------------------------------
  Assign("nvtx",nvtx);
  Assign("bveto_ip",bveto_ip);

  //Creating the variables we need
  //----------------------------------------------------------------------------
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * _Luminosity;
  
  if (_Signal.Contains("Data"))
    totalW = 1.0;

  h_n_PV -> Fill(1,efficiencyW);  
  
  dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
  int jetbin = njet;
  
  //Building mpmet
  Float_t dphimin = (min(dphilmet1,dphilmet2));
  
  if (dphimin < TMath::Pi() / 2)
    fullpmet = pfType1Met * sin(dphimin);
  else 
    fullpmet = pfType1Met;
  
  if (dphimin < TMath::Pi() / 2)
    trkpmet = trkMet * sin(dphimin);
  else
    trkpmet = trkMet;
  
  mpmet = min(trkpmet,fullpmet);
  
  metvar = (njet <= 1) ? mpmet : pfType1Met;
  
  //building Ht
  Ht = std_vector_lepton_pt.at(0) + std_vector_lepton_pt.at(1) + pfType1Met;
  
  if(njet > 10) njet = 10;
  for (int i = 0; i < njet; ++i)
    if(std_vector_jet_pt.at(i) > 0)
      Ht += std_vector_jet_pt.at(i);

  //building ptWW
  TLorentzVector L1,L2;
  TLorentzVector MET;
  
  if (pt1>0 && pt2>0) {  
    L1.SetPtEtaPhiM(pt1, 0, phi1, 0.);
    L2.SetPtEtaPhiM(pt2, 0, phi2, 0.);
    MET.SetPtEtaPhiM(pfType1Met, 0, pfType1Metphi, 0.);
    ptWW = (L1+L2+MET).Pt();
  }

  //defining nextra
  nextra = 0;
  for (unsigned int i = 0; i < std_vector_lepton_pt.size(); ++i)
    if(std_vector_lepton_pt.at(i) > 10)
      ++nextra;
  nextra = nextra - 2;
       
  //building Mt of the two Ws
  if( pfType1Met > 0){
    if( pt1 > 0 )
      Mt1 = sqrt(2*pt1*pfType1Met*(1 - cos(dphilmet1)));  
    if( pt2 > 0 )
      Mt2 = sqrt(2*pt2*pfType1Met*(1 - cos(dphilmet2)));  
  }
  
  //building Mc
  if (ptll > 0 && mll > 0 && pfType1Met > 0)
    Mc = sqrt( pow(sqrt(ptll*ptll + mll*mll) + pfType1Met,2) - pow(ptll + pfType1Met,2) );
  
  //building dphijetjet
  if (jetpt1 > 0 && jetpt2 > 0)
    if (jetphi1 > -3.5 && jetphi1 < 3.5)
      if (jetphi2 > -3.5 && jetphi2 < 3.5){
	dphijetjet = fabs(jetphi2 - jetphi1);
	if (dphijetjet > TMath::Pi())
	  dphijetjet = 2*TMath::Pi() - dphijetjet;
      }

  //building detajetjet
  if (jetpt1 > 0 && jetpt2 > 0)
    detajetjet = fabs(jeteta1 - jeteta2);

  //Building dphijet1met
  dphijet1met = 0.;
  if (jetphi1 > 0 && pfType1Metphi > 0){
    dphijet1met = fabs(jetphi1 - pfType1Metphi);
    if (dphijet1met > TMath::Pi()) dphijet1met = 2*TMath::Pi() - dphijet1met;
  }
  
  //Building RatioMet   
  ratioMet = 0.;
  if (pfType1Met > 0 && trkMet > 0)
    ratioMet = pfType1Met / sqrt(pfType1Met + trkMet);
  
   //Building b-veto csvv2ivf Loose (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfLoose = false;

   for (int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.423)
       bvetocsvv2ivfLoose = true;

   //Building b-veto csvv2ivf Medium (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfMedium = false;

   for (int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.814)
       bvetocsvv2ivfMedium = true;

   //Building b-veto csvv2ivf Tight (false: no b-jet -> good!! true: b-jet detected -> reject event!!)
   bool bvetocsvv2ivfTight = false;

   for (int i = 0; i < njet; ++i)
     if (std_vector_jet_csvv2ivf.at(i) > 0.941)
       bvetocsvv2ivfTight = true;

   //Building a Soft Muon B-Veto (Quite Rudely)
   bool bvetoMu = false;
   
   for (unsigned int i = 0; i < std_vector_jet_softMuPt.size(); ++i)
     if (std_vector_jet_pt.at(i) > 10 && std_vector_jet_pt.at(i) < 30)
       if (std_vector_jet_softMuPt.at(i) > 3){
	 //hsoftMuPt->Fill(std_vector_jet_softMuPt.at(i), totalW);
	 //hjetPt   ->Fill(std_vector_jet_pt.at(i),       totalW);
	 bvetoMu = true;
	 break;
       }

   //Building Number of Gen Jet
   njetGen = 0;
   /*
   for (int i = 0; i < std_vector_jetGen_pt.size(); ++i)
   if (std_vector_jetGen_pt.at(i) > 30)
   ++njetGen;       
   */

   // The selection begins here
   //--------------------------------------------------------------------------
   if (trigger == 1)
     if (std_vector_lepton_pt.at(0) > 20)
       if (std_vector_lepton_pt.at(1) > 20) 
	 if ((_SameSign == "SS" && ch1*ch2 > 0) || (_SameSign == "OS" && ch1*ch2 < 0))
	   if ( (SelectedChannel == -1)                                     || 
		(channel == SelectedChannel)                                || 
		(SelectedChannel == 4 && (channel == 2 || channel == 3) )   || 
		(SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
		){
	     /*	    
		    if (IsTightLepton(0,_MuonID) && !IsTightLepton(1,_MuonID))
		    hLooseIso -> Fill(ElectronIsolation(1), totalW);
		    if (IsTightLepton(1,_MuonID) && !IsTightLepton(0,_MuonID))
		    hLooseIso -> Fill(ElectronIsolation(0), totalW);
	     */
	     //if (IsIsolatedLepton(0))
	     //if (IsIsolatedLepton(1))
	     //if (IsTightLepton(0,_MuonID))
	     //if (IsTightLepton(1,_MuonID))
	     {
		     
		     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		     //
		     // Main analisis
		     //
		     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		     
		     hWTrigger   ->Fill(1, totalW); 
		     hWeffTrigger->Fill(1, efficiencyW);
		     
		    hPtLepton1TwoLeptonsLevel      ->Fill(pt1,       totalW);
		    hPtLepton2TwoLeptonsLevel      ->Fill(pt2,       totalW);
		    hPtDiLeptonTwoLeptonsLevel     ->Fill(ptll,      totalW);
		    hMinvTwoLeptonsLevel           ->Fill(mll,       totalW);
		    hMtTwoLeptonsLevel             ->Fill(mth,       totalW);
		    hMt1TwoLeptonsLevel            ->Fill(Mt1,       totalW);
		    hMt2TwoLeptonsLevel            ->Fill(Mt2,       totalW);
		    hpfMetTwoLeptonsLevel          ->Fill(pfType1Met,totalW);
		    hpminMetTwoLeptonsLevel        ->Fill(mpmet,     totalW);
		    hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,      totalW);
		    hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,    totalW);
		    hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet, totalW);
		    hDPhillMetTwoLeptonsLevel      ->Fill(dphillmet, totalW);		  
		    hMcTwoLeptonsLevel             ->Fill(Mc,        totalW);		  
		    hPtWWTwoLeptonsLevel           ->Fill(ptWW,      totalW);		  
		    hTrkMetTwoLeptonsLevel         ->Fill(trkMet,    totalW);
		    hHtTwoLeptonsLevel             ->Fill(Ht,        totalW);			
		    hdphijetjetTwoLeptonsLevel     ->Fill(dphijetjet,totalW);
		    hdetajetjetTwoLeptonsLevel     ->Fill(detajetjet,totalW);
		     
		    //WW Level Selections
		    if (nextra == 0) {
		      
		      hWExtraLepton->Fill(1, totalW);
		      hWeffExtraLepton->Fill(1, efficiencyW);
		      
		      if (pfType1Met > 20 ) {
			
			hWMetCut->Fill(1, totalW);
			hWeffMetCut->Fill(1, efficiencyW);
			
			if (mll > 12) {
			  
			  hWLowMinv->Fill(1, totalW);
			  hWeffLowMinv->Fill(1, efficiencyW);
			  
			   //zveto (in case of same flavour)
			   if ( (fabs(ZMASS - mll) > 15 && 
				 ( metvar > 45 ) )      ||
				channel == 2            ||
				channel == 3            ){
			     
			     hWZVeto->Fill(1, totalW);
			     hWeffZVeto->Fill(1, efficiencyW);
			     
			     if (mpmet > 20){
			       
			       hWpMetCut->Fill(1, totalW);
			       hWeffpMetCut->Fill(1, efficiencyW);
			       
			       if (dphiv || channel == 2 || channel == 3) {
				 
				 hWDeltaPhiJet->Fill(1, totalW);
				 hWeffDeltaPhiJet->Fill(1, efficiencyW);
				 
				 if ( ptll>30 && (channel == 2 || channel == 3 || ptll>45) ) {
				   
				   hWPtll->Fill(1, totalW); 	   
				   hWeffPtll->Fill(1, efficiencyW);			    				   
				   hWnJets->Fill(njet, totalW);
				   hWeffnJets->Fill(njet, efficiencyW);
				   
				   hWnBtaggedJets->Fill(nbjet, totalW);
				   hWeffnBtaggedJets->Fill(nbjet, efficiencyW);
				   
				  //b-veto
				  if (bveto_ip == 1 && nbjettche == 0) {
				    
				    hWTopTagging->Fill(1, totalW);
				    hWeffTopTagging->Fill(1, efficiencyW);
				    
				    hWSoftMuVeto->Fill(1, totalW);
				    hWeffSoftMuVeto->Fill(1,efficiencyW);
				    
				    hPtLepton1WWLevel      ->Fill(pt1,       totalW);
				    hPtLepton2WWLevel      ->Fill(pt2,       totalW);
				    hPtDiLeptonWWLevel     ->Fill(ptll,      totalW);
				    hMinvWWLevel           ->Fill(mll,       totalW);
				    hMtWWLevel             ->Fill(mth,       totalW);
				    hMt1WWLevel            ->Fill(Mt1,       totalW);
				    hMt2WWLevel            ->Fill(Mt2,       totalW);
				    hpfMetWWLevel          ->Fill(pfType1Met,totalW);
				    hpminMetWWLevel        ->Fill(mpmet,     totalW);
				    hDeltaRLeptonsWWLevel  ->Fill(drll,      totalW);
				    hDeltaPhiLeptonsWWLevel->Fill(dphill,    totalW);
				    hDPhiPtllJetWWLevel    ->Fill(dphilljet, totalW);
				    hDPhillMetWWLevel      ->Fill(dphillmet, totalW);		  
				    hMcWWLevel             ->Fill(Mc,        totalW);		  
				    hPtWWWWLevel           ->Fill(ptWW,      totalW);		  
				    hTrkMetWWLevel         ->Fill(trkMet,    totalW);
				    hHtWWLevel             ->Fill(Ht,        totalW);			
				    hdphijetjetWWLevel     ->Fill(dphijetjet,totalW);
				    hdetajetjetWWLevel     ->Fill(detajetjet,totalW);
				    
				    //Control Region Selections - on Top of WW Cuts
				    if (drll > 1.5){
				      hpfMetCR  -> Fill (pfType1Met, totalW);
				      htrkMetCR -> Fill (trkMet,     totalW);
				      hmpMetCR  -> Fill (mpmet,      totalW);
				    }
				    
				    //monoHiggs Selections
				    if (Mc < 100){
				      
				      hPtLepton1monoHLevel[0]      ->Fill(pt1,       totalW);
				      hPtLepton2monoHLevel[0]      ->Fill(pt2,       totalW);
				      hPtDiLeptonmonoHLevel[0]     ->Fill(ptll,      totalW);
				      hMinvmonoHLevel[0]           ->Fill(mll,       totalW);
				      hMtmonoHLevel[0]             ->Fill(mth,       totalW);
				      hMt1monoHLevel[0]            ->Fill(Mt1,       totalW);
				      hMt2monoHLevel[0]            ->Fill(Mt2,       totalW);
				      hpfMetmonoHLevel[0]          ->Fill(pfType1Met,totalW);
				      hpminMetmonoHLevel[0]        ->Fill(mpmet,     totalW);
				      hDeltaRLeptonsmonoHLevel[0]  ->Fill(drll,      totalW);
				      hDeltaPhiLeptonsmonoHLevel[0]->Fill(dphill,    totalW);
				      hDPhiPtllJetmonoHLevel[0]    ->Fill(dphilljet, totalW);
				      hDPhillMetmonoHLevel[0]      ->Fill(dphillmet, totalW);
				      hPtWWmonoHLevel[0]           ->Fill(ptWW,      totalW);
				      hMcmonoHLevel[0]             ->Fill(Mc,        totalW);		  
				      hTrkMetmonoHLevel[0]         ->Fill(trkMet,    totalW);
				      hHtmonoHLevel[0]             ->Fill(Ht,        totalW);			
				      hdphijetjetmonoHLevel[0]     ->Fill(dphijetjet,totalW);
				      hdetajetjetmonoHLevel[0]     ->Fill(detajetjet,totalW);
				      
				      if (drll < 1.5){
					
					hPtLepton1monoHLevel[1]      ->Fill(pt1,       totalW);
					hPtLepton2monoHLevel[1]      ->Fill(pt2,       totalW);
					hPtDiLeptonmonoHLevel[1]     ->Fill(ptll,      totalW);
					hMinvmonoHLevel[1]           ->Fill(mll,       totalW);
					hMtmonoHLevel[1]             ->Fill(mth,       totalW);
					hMt1monoHLevel[1]            ->Fill(Mt1,       totalW);
					hMt2monoHLevel[1]            ->Fill(Mt2,       totalW);
					hpfMetmonoHLevel[1]          ->Fill(pfType1Met,totalW);
					hpminMetmonoHLevel[1]        ->Fill(mpmet,     totalW);
					hDeltaRLeptonsmonoHLevel[1]  ->Fill(drll,      totalW);
					hDeltaPhiLeptonsmonoHLevel[1]->Fill(dphill,    totalW);
					hDPhiPtllJetmonoHLevel[1]    ->Fill(dphilljet, totalW);
					hDPhillMetmonoHLevel[1]      ->Fill(dphillmet, totalW);
					hPtWWmonoHLevel[1]           ->Fill(ptWW,      totalW);
					hMcmonoHLevel[1]             ->Fill(Mc,        totalW);		  
					hTrkMetmonoHLevel[1]         ->Fill(trkMet,    totalW);
					hHtmonoHLevel[1]             ->Fill(Ht,        totalW);			
					hdphijetjetmonoHLevel[1]     ->Fill(dphijetjet,totalW);
					hdetajetjetmonoHLevel[1]     ->Fill(detajetjet,totalW);
					
					if ( mpmet > 60 ){
					  
					  hPtLepton1monoHLevel[2]      ->Fill(pt1,       totalW);
					  hPtLepton2monoHLevel[2]      ->Fill(pt2,       totalW);
					  hPtDiLeptonmonoHLevel[2]     ->Fill(ptll,      totalW);
					  hMinvmonoHLevel[2]           ->Fill(mll,       totalW);
					  hMtmonoHLevel[2]             ->Fill(mth,       totalW);
					  hMt1monoHLevel[2]            ->Fill(Mt1,       totalW);
					  hMt2monoHLevel[2]            ->Fill(Mt2,       totalW);
					  hpfMetmonoHLevel[2]          ->Fill(pfType1Met,totalW);
					  hpminMetmonoHLevel[2]        ->Fill(mpmet,     totalW);
					  hDeltaRLeptonsmonoHLevel[2]  ->Fill(drll,      totalW);
					  hDeltaPhiLeptonsmonoHLevel[2]->Fill(dphill,    totalW);
					  hDPhiPtllJetmonoHLevel[2]    ->Fill(dphilljet, totalW);
					  hDPhillMetmonoHLevel[2]      ->Fill(dphillmet, totalW);
					  hPtWWmonoHLevel[2]           ->Fill(ptWW,      totalW);
					  hMcmonoHLevel[2]             ->Fill(Mc,        totalW);		  
					  hTrkMetmonoHLevel[2]         ->Fill(trkMet,    totalW);
					  hHtmonoHLevel[2]             ->Fill(Ht,        totalW);						      
					  hdphijetjetmonoHLevel[2]     ->Fill(dphijetjet,totalW);
					  hdetajetjetmonoHLevel[2]     ->Fill(detajetjet,totalW);
					  
					  mvaTree->Fill();
					  
					  
					}
				      }  					
				    }
				  }
				 }
			       }
			     }
			   }
			}
		      }
		    }
		   }
	   }
   
   
   // Define Normalization Factor for MC samples 
   //------------------------------------------------------------------------------
   
   // Define weights
   //------------------------------------------------------------------------------
   
   float pileupweight = 1;
   
   //  if (!IsDATA)
   // pileupweight = fPUWeight->GetWeight(T_Event_nPU);
   
   double factN = 1;
   
   if (_XSection > 0) factN = _XSection * _Luminosity / _NEvents;
   
   //factN = factN*pileupweight;
   //------------------------------------------------------------------------------
   
   // Init variables
   //------------------------------------------------------------------------------
   
} // end inside Loop

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
//



void monoHiggsSelector::Summary() {
  // Get Data Members at the client-master (after finishing the analysis at the workers nodes)
  // Only data members set here will be accesible at the client-master

  mvaTree = FindOutput<TTree*>("nt");
  ///*** 1D histos ***/// 

  h_n_PV  = FindOutput<TH1F*>("h_N_PV");  
  
  // Counting histograms    
  //---------------------------------------------------------------------------- 
  
  hWExtraLepton     = FindOutput<TH1F*>("hWExtraLepton");         
  hWeffExtraLepton  = FindOutput<TH1F*>("hWeffExtraLepton");  
  hWMetCut          = FindOutput<TH1F*>("hWMetCut");          
  hWeffMetCut       = FindOutput<TH1F*>("hWeffMetCut");       
  hWLowMinv         = FindOutput<TH1F*>("hWLowMinv");         
  hWeffLowMinv      = FindOutput<TH1F*>("hWeffLowMinv");      
  hWZVeto           = FindOutput<TH1F*>("hWZVeto");           
  hWeffZVeto        = FindOutput<TH1F*>("hWeffZVeto");        
  hWpMetCut         = FindOutput<TH1F*>("hWpMetCut");         
  hWeffpMetCut      = FindOutput<TH1F*>("hWeffpMetCut");      
  hWDeltaPhiJet     = FindOutput<TH1F*>("hWDeltaPhiJet");     
  hWeffDeltaPhiJet  = FindOutput<TH1F*>("hWeffDeltaPhiJet");  
  hWPtll            = FindOutput<TH1F*>("hWPtll");            
  hWeffPtll         = FindOutput<TH1F*>("hWeffPtll");         
  hWnJets           = FindOutput<TH1F*>("hWnJets");           
  hWeffnJets        = FindOutput<TH1F*>("hWeffnJets");        
  hWnBtaggedJets    = FindOutput<TH1F*>("hWnBtaggedJets");    
  hWeffnBtaggedJets = FindOutput<TH1F*>("hWeffnBtaggedJets"); 
  hWTopTagging      = FindOutput<TH1F*>("hWTopTagging");      
  hWeffTopTagging   = FindOutput<TH1F*>("hWeffTopTagging");   
  hWSoftMuVeto      = FindOutput<TH1F*>("hWSoftMuVeto");      
  hWeffSoftMuVeto   = FindOutput<TH1F*>("hWeffSoftMuVeto");   

  // TwoLeptons level histograms
  //---------------------------------------------------------------------------

  hPtLepton1TwoLeptonsLevel       = FindOutput<TH1F*>("hPtLepton1TwoLeptonsLevel");       
  hPtLepton2TwoLeptonsLevel       = FindOutput<TH1F*>("hPtLepton2TwoLeptonsLevel");       
  hPtDiLeptonTwoLeptonsLevel      = FindOutput<TH1F*>("hPtDiLeptonTwoLeptonsLevel");      
  hMinvTwoLeptonsLevel            = FindOutput<TH1F*>("hMinvTwoLeptonsLevel");            
  hMtTwoLeptonsLevel              = FindOutput<TH1F*>("hMtTwoLeptonsLevel");              
  hMt1TwoLeptonsLevel             = FindOutput<TH1F*>("hMt1TwoLeptonsLevel");             
  hMt2TwoLeptonsLevel             = FindOutput<TH1F*>("hMt2TwoLeptonsLevel");             
  hpfMetTwoLeptonsLevel           = FindOutput<TH1F*>("hpfMetTwoLeptonsLevel");           
  hpminMetTwoLeptonsLevel         = FindOutput<TH1F*>("hpminMetTwoLeptonsLevel");         
  hDeltaRLeptonsTwoLeptonsLevel   = FindOutput<TH1F*>("hDeltaRLeptonsTwoLeptonsLevel");   
  hDeltaPhiLeptonsTwoLeptonsLevel = FindOutput<TH1F*>("hDeltaPhiLeptonsTwoLeptonsLevel"); 
  hDPhiPtllJetTwoLeptonsLevel     = FindOutput<TH1F*>("hDPhiPtllJetTwoLeptonsLevel");     
  hDPhillMetTwoLeptonsLevel       = FindOutput<TH1F*>("hDPhillMetTwoLeptonsLevel");       
  hMcTwoLeptonsLevel              = FindOutput<TH1F*>("hMcTwoLeptonsLevel");               
  hPtWWTwoLeptonsLevel            = FindOutput<TH1F*>("hPtWWTwoLeptonsLevel");            
  hTrkMetTwoLeptonsLevel          = FindOutput<TH1F*>("hTrkMetTwoLeptonsLevel");          
  hHtTwoLeptonsLevel              = FindOutput<TH1F*>("hHtTwoLeptonsLevel");              
  hdphijetjetTwoLeptonsLevel      = FindOutput<TH1F*>("hdphijetjetTwoLeptonsLevel");      
  hdetajetjetTwoLeptonsLevel      = FindOutput<TH1F*>("hdetajetjetTwoLeptonsLevel");      

  // WW level histograms
  //----------------------------------------------------------------------------                                                    

  hPtLepton1WWLevel       = FindOutput<TH1F*>("hPtLepton1WWLevel");       
  hPtLepton2WWLevel       = FindOutput<TH1F*>("hPtLepton2WWLevel");       
  hPtDiLeptonWWLevel      = FindOutput<TH1F*>("hPtDiLeptonWWLevel");      
  hMinvWWLevel            = FindOutput<TH1F*>("hMinvWWLevel");            
  hMtWWLevel              = FindOutput<TH1F*>("hMtWWLevel");              
  hMt1WWLevel             = FindOutput<TH1F*>("hMt1WWLevel");             
  hMt2WWLevel             = FindOutput<TH1F*>("hMt2WWLevel");             
  hpfMetWWLevel           = FindOutput<TH1F*>("hpfMetWWLevel");           
  hpminMetWWLevel         = FindOutput<TH1F*>("hpminMetWWLevel");         
  hDeltaRLeptonsWWLevel   = FindOutput<TH1F*>("hDeltaRLeptonsWWLevel");   
  hDeltaPhiLeptonsWWLevel = FindOutput<TH1F*>("hDeltaPhiLeptonsWWLevel"); 
  hDPhiPtllJetWWLevel     = FindOutput<TH1F*>("hDPhiPtllJetWWLevel");     
  hDPhillMetWWLevel       = FindOutput<TH1F*>("hDPhillMetWWLevel");       
  hMcWWLevel              = FindOutput<TH1F*>("hMcWWLevel");                
  hPtWWWWLevel            = FindOutput<TH1F*>("hPtWWWWLevel");            
  hTrkMetWWLevel          = FindOutput<TH1F*>("hTrkMetWWLevel");          
  hHtWWLevel              = FindOutput<TH1F*>("hHtWWLevel");              
  hdphijetjetWWLevel      = FindOutput<TH1F*>("hdphijetjetWWLevel");      
  hdetajetjetWWLevel      = FindOutput<TH1F*>("hdetajetjetWWLevel");      
    
  //Control Region Plots
  //----------------------------------------------------------------------------                

  hpfMetCR  = FindOutput<TH1F*>("hpfMetCR");           
  htrkMetCR = FindOutput<TH1F*>("htrkMetCR");          
  hmpMetCR  = FindOutput<TH1F*>("hmpMetCR");           

  // monoH level histograms     
  //----------------------------------------------------------------------------  

  char name[80];

  for (Int_t qq=0; qq<4; qq++) {

    sprintf(name,"hPtLepton1monoHLevel%.1i",qq);
    hPtLepton1monoHLevel[qq]       = FindOutput<TH1F*>(name); 
    sprintf(name,"hPtLepton2monoHLevel%.1i",qq);
    hPtLepton2monoHLevel[qq]       = FindOutput<TH1F*>(name); 
    sprintf(name,"hPtDiLeptonmonoHLevel%.1i",qq);
    hPtDiLeptonmonoHLevel[qq]      = FindOutput<TH1F*>(name); 
    sprintf(name,"hMinvmonoHLevel%.1i",qq);
    hMinvmonoHLevel[qq]            = FindOutput<TH1F*>(name); 
    sprintf(name,"hMtmonoHLevel%.1i",qq);
    hMtmonoHLevel[qq]              = FindOutput<TH1F*>(name); 
    sprintf(name,"hMt1monoHLevel%.1i",qq);
    hMt1monoHLevel[qq]             = FindOutput<TH1F*>(name); 
    sprintf(name,"hMt2monoHLevel%.1i",qq);
    hMt2monoHLevel[qq]             = FindOutput<TH1F*>(name); 
    sprintf(name,"hpfMetmonoHLevel%.1i",qq);
    hpfMetmonoHLevel[qq]           = FindOutput<TH1F*>(name); 
    sprintf(name,"hpminMetmonoHLevel%.1i",qq);
    hpminMetmonoHLevel[qq]         = FindOutput<TH1F*>(name); 
    sprintf(name,"hDeltaRLeptonsmonoHLevel%.1i",qq);
    hDeltaRLeptonsmonoHLevel[qq]   = FindOutput<TH1F*>(name); 
    sprintf(name,"hDeltaPhiLeptonsmonoHLevel%.1i",qq);
    hDeltaPhiLeptonsmonoHLevel[qq] = FindOutput<TH1F*>(name); 
    sprintf(name,"hDPhiPtllJetmonoHLevel%.1i",qq);
    hDPhiPtllJetmonoHLevel[qq]     = FindOutput<TH1F*>(name); 
    sprintf(name,"hDPhillMetmonoHLevel%.1i",qq);
    hDPhillMetmonoHLevel[qq]       = FindOutput<TH1F*>(name); 
    sprintf(name,"hPtWWmonoHLevel%.1i",qq);
    hPtWWmonoHLevel[qq]            = FindOutput<TH1F*>(name); 
    sprintf(name,"hMcmonoHLevel%.1i",qq);
    hMcmonoHLevel[qq]              = FindOutput<TH1F*>(name); 
    sprintf(name,"hTrkMetmonoHLevel%.1i",qq);
    hTrkMetmonoHLevel[qq]          = FindOutput<TH1F*>(name); 
    sprintf(name,"hHtmonoHLevel%.1i",qq);
    hHtmonoHLevel[qq]              = FindOutput<TH1F*>(name); 
    sprintf(name,"hdphijetjetmonoHLevel%.1i",qq);
    hdphijetjetmonoHLevel[qq]      = FindOutput<TH1F*>(name); 
    sprintf(name,"hdetajetjetmonoHLevel%.1i",qq);
    hdetajetjetmonoHLevel[qq]      = FindOutput<TH1F*>(name); 
  }
}

//------------------------------------------------------------------------------
// IsTightLepton
//
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
// egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium
//------------------------------------------------------------------------------
bool monoHiggsSelector::IsTightLepton(int k, TString _MuonID_)
{
  bool is_tight_lepton = false;

  // Muon ID
  if (fabs(std_vector_lepton_flavour.at(k)) == 13){
    
    if (_MuonID_ == "MediumID"){
      if (std_vector_lepton_isMediumMuon.at(k) == 1)
	is_tight_lepton = true;
    }
    
    else if(_MuonID_ == "MediumIDTighterIP"){
      if (std_vector_lepton_isMediumMuon.at(k) == 1)
	if (fabs(std_vector_lepton_BestTrackdxy.at(k)) < 0.02)
	  if (fabs(std_vector_lepton_BestTrackdz.at(k)) < 0.1)
	    is_tight_lepton = true;
    }
    
    else if (_MuonID_ == "TightID" ){
      if (std_vector_lepton_isTightMuon.at(k) == 1)
	is_tight_lepton = true;
    }
    
    else if(_MuonID_ == "TightIDTighterIP"){
      if (std_vector_lepton_isTightMuon.at(k) == 1)
	if (fabs(std_vector_lepton_BestTrackdxy.at(k)) < 0.02)
	  if (fabs(std_vector_lepton_BestTrackdz.at(k)) < 0.1)
	    is_tight_lepton = true;
    }
  }
  
  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_flavour.at(k)) == 11)
    {
      is_tight_lepton = std_vector_lepton_eleIdMedium.at(k);
    }
  /*   
  float aeta = fabs(std_vector_electron_scEta.at(k));
      
      if (aeta <= 1.479)
	{
	  if (fabs(std_vector_electron_deltaEtaIn.at(k)) < 0.008925 &&
	      fabs(std_vector_electron_deltaPhiIn.at(k)) < 0.035973 &&
	      std_vector_electron_sigmaIetaIeta.at(k)    < 0.009996 &&
	      std_vector_electron_HoE.at(k)              < 0.050537 &&
	      fabs(std_vector_electron_d0.at(k))         < 0.012235 &&
	      fabs(std_vector_electron_dz.at(k))         < 0.042020 &&
	      fabs(std_vector_electron_ooEooP.at(k))     < 0.091942 &&
	      ElectronIsolation(k)                        < 0.107587 &&
	      !std_vector_electron_passConversion.at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
      else if (aeta > 1.479 && aeta < 2.5)
	{
	  if (fabs(std_vector_electron_deltaEtaIn.at(k)) < 0.007429 &&
	      fabs(std_vector_electron_deltaPhiIn.at(k)) < 0.067879 &&
	      std_vector_electron_sigmaIetaIeta.at(k)    < 0.030135 &&
	      std_vector_electron_HoE.at(k)              < 0.086782 &&
	      fabs(std_vector_electron_d0.at(k))         < 0.036719 &&
	      fabs(std_vector_electron_dz.at(k))         < 0.138142 &&
	      fabs(std_vector_electron_ooEooP.at(k))     < 0.100683 &&
	      ElectronIsolation(k)                        < 0.113254 &&
	      !std_vector_electron_passConversion.at(k))  // Includes expectedMissingInnerHits
	    {
	      is_tight_lepton = true;
	    }
	}
    }
  */
  return is_tight_lepton;

}


//------------------------------------------------------------------------------
// MuonIsolation
//------------------------------------------------------------------------------
float monoHiggsSelector::MuonIsolation(int k)
{
  float pt = std_vector_lepton_pt.at(k);
  float id = std_vector_lepton_flavour.at(k);

  float relative_isolation = -999;

  if (fabs(id) != 13) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso.at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso.at(k) +
	      std_vector_lepton_neutralHadronIso.at(k) -
	      0.5*std_vector_lepton_sumPUPt.at(k)));

  relative_isolation /= pt;

  return relative_isolation;
}


//------------------------------------------------------------------------------
// ElectronIsolation
//------------------------------------------------------------------------------
float monoHiggsSelector::ElectronIsolation(int k)
{
  float pt = std_vector_lepton_pt.at(k);
  float id = std_vector_lepton_flavour.at(k);

  float relative_isolation = -999;

  if (fabs(id) != 11) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso.at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso.at(k) +
	      std_vector_lepton_neutralHadronIso.at(k) -
	      jetRho*std_vector_electron_effectiveArea.at(k)));
  
  relative_isolation /= pt;
  
  return relative_isolation;
}


//------------------------------------------------------------------------------
// IsIsolatedLepton
//------------------------------------------------------------------------------
bool monoHiggsSelector::IsIsolatedLepton(int k)
{
  float id = std_vector_lepton_flavour.at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.12);
  
  return is_isolated_lepton;
}
