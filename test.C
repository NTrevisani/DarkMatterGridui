////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////           DUMMY                                      /////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


#include "test.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"

#include "TLorentzVector.h"
#include <vector>
#include "TROOT.h"
#include <iostream>

#include "TDatabasePDG.h"


test::test(TTree* tree):
  PAFAnalysis(tree) {
}


void test::Initialise() {

  GetInputParameters()->TheNamedBool("IsDATA", IsDATA);
  GetInputParameters()->TheNamedInt("NEvents", NEvents);
  GetInputParameters()->TheNamedDouble("Luminosity", Luminosity);
  GetInputParameters()->TheNamedDouble("XSection", XSection);
  GetInputParameters()->TheNamedInt("WhichRun", WhichRun);
  GetInputParameters()->TheNamedInt("jetChannel", jetChannel);
  
  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  TString flavorChannel = GetInputParameters()->TheNamedString("FlavorChannel");
  
  cout<<"la stringa che voglio passare si chiama: "<<TheSample<<endl;
  
  
  // To do only once
  float LuminosityPU = 0;
  //fInputParameters->TheNamedFloat("LuminosityPU",LuminosityPU);
  //fPUWeight = new PUWeight(LuminosityPU, Spring11);//Summer11InTime);
  
  // Set the channel                                                                                                                                                               
  //----------------------------------------------------------------------------                                                                                                   
  SelectedChannel = -999;
  
  if      (flavorChannel == "MuMu") SelectedChannel =  0;
  else if (flavorChannel == "EE"  ) SelectedChannel =  1;
  else if (flavorChannel == "EMu" ) SelectedChannel =  2;
  else if (flavorChannel == "MuE" ) SelectedChannel =  3;
  else if (flavorChannel == "OF" )  SelectedChannel =  4;
  else if (flavorChannel == "SF" )  SelectedChannel =  5;
  
  else if (flavorChannel == "All" ) SelectedChannel = -1;
  
  G_Debug_DefineAnalysisVariables = false;
  
  //------------------------------------------------------------------------------
  // Create mvaTree branches
  //------------------------------------------------------------------------------

  mvaTree = CreateTree("mvaTree","mvaTree");

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

// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {
  
  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
  
  h_n_PV -> Fill(1,efficiencyW);  
  
  Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
  Float_t jetbin = njet;
  
  //building mpmet
  Float_t dphimin = (min(dphilmet1,dphilmet2));
  Float_t fullpmet = 0;
  Float_t trkpmet = 0;
  
  if (dphimin < TMath::Pi() / 2)
    fullpmet = fabs(pfType1Met * sin(dphimin));
  else 
    fullpmet = pfType1Met;
  
  if (dphimin < TMath::Pi() / 2)
    trkpmet = fabs(trkMet * sin(dphimin));
  else
    trkpmet = trkMet;
  
  mpmet = min(trkpmet,fullpmet);
  
  Float_t metvar = (njet <= 1) ? mpmet : pfType1Met;
  
  //building Ht
  Ht = std_vector_lepton_pt->at(0) + std_vector_lepton_pt->at(1) + pfType1Met;
  
  if(njet > 10) njet = 10;
  for (int i = 0; i < njet; ++i)
    if(std_vector_jet_pt->at(i) > 0)
      Ht += std_vector_jet_pt->at(i);
  
  //building ptWW
  TLorentzVector L1,L2;
  TLorentzVector MET;
  
  if (pt1>0 && pt2>0) {  
    L1.SetPtEtaPhiM(pt1, 0, phi1, 0.);
    L2.SetPtEtaPhiM(pt2, 0, phi2, 0.);
    MET.SetPtEtaPhiM(pfType1Met, 0, pfType1Metphi, 0.);
    ptWW = (L1+L2+MET).Pt();
  }
  
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

  // The selection begins here
  //--------------------------------------------------------------------------
  if (std_vector_lepton_pt->at(0) > 20)
    if (std_vector_lepton_pt->at(1) > 20) 
      if (ch1*ch2 < 0)
	//if (std_vector_lepton_BestTrackdxy -> at(0) < 0.01 && std_vector_lepton_BestTrackdz -> at(0) < 0.1)
	//if (std_vector_lepton_BestTrackdxy -> at(1) < 0.01 && std_vector_lepton_BestTrackdz -> at(1) < 0.1)
	//if (prompt == 1)
	//if (std_vector_lepton_isTightMuon->at(0) == 1)// && std_vector_lepton_isTightMuon->at(1) == 0 || 
	//if (std_vector_lepton_isTightMuon->at(1) == 1)// && std_vector_lepton_isTightMuon->at(1) == 1)
	//if (isoZero < 0.12)
	//if (isoOne < 0.12)
	if ( (SelectedChannel == -1)                                   || 
	     (channel == SelectedChannel)                              || 
	     (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
	     (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
	     ){
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //
	  // Main analisis
	  //
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  
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
	  
	  cout<<"hola"<<endl;

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

  // Define Normalization Factor for MC samples 
  
  //------------------------------------------------------------------------------
  // Define weights
  //------------------------------------------------------------------------------
  
  float pileupweight = 1;
  
  //  if (!IsDATA)
  // pileupweight = fPUWeight->GetWeight(T_Event_nPU);
  
  double factN = 1;
  
  if (XSection > 0) factN = XSection * Luminosity / NEvents;
  
  
  //factN = factN*pileupweight;
  
  
  //------------------------------------------------------------------------------
  // Init variables
  //------------------------------------------------------------------------------
  
}// end inside Loop

void test::SetDataMembersAtTermination() {

  mvaTree = ((TTree*) FindOutput("mvaTree")); 
  
  h_n_PV = ((TH1F*) FindOutput("h_n_PV")); 
  
  // Counting histograms                                                    
  //----------------------------------------------------------------------------                         
  
  hWExtraLepton     = ((TH1F*) FindOutput("hWExtraLepton"));   
  hWeffExtraLepton  = ((TH1F*) FindOutput("hWeffExtraLepton"));
  hWMetCut          = ((TH1F*) FindOutput("hWMetCut"));
  hWeffMetCut       = ((TH1F*) FindOutput("hWeffMetCut"));
  hWLowMinv         = ((TH1F*) FindOutput("hWLowMinv"));
  hWeffLowMinv      = ((TH1F*) FindOutput("hWeffLowMinv"));
  hWZVeto           = ((TH1F*) FindOutput("hWZVeto"));
  hWeffZVeto        = ((TH1F*) FindOutput("hWeffZVeto"));
  hWpMetCut         = ((TH1F*) FindOutput("hWpMetCut"));
  hWeffpMetCut      = ((TH1F*) FindOutput("hWeffpMetCut"));
  hWDeltaPhiJet     = ((TH1F*) FindOutput("hWDeltaPhiJet"));
  hWeffDeltaPhiJet  = ((TH1F*) FindOutput("hWeffDeltaPhiJet"));
  hWPtll            = ((TH1F*) FindOutput("hWPtll"));
  hWeffPtll         = ((TH1F*) FindOutput("hWeffPtll"));
  hWnJets           = ((TH1F*) FindOutput("hWnJets"));
  hWeffnJets        = ((TH1F*) FindOutput("hWeffnJets"));
  hWnBtaggedJets    = ((TH1F*) FindOutput("hWnBtaggedJets"));
  hWeffnBtaggedJets = ((TH1F*) FindOutput("hWeffnBtaggedJets"));
  hWTopTagging      = ((TH1F*) FindOutput("hWTopTagging"));
  hWeffTopTagging   = ((TH1F*) FindOutput("hWeffTopTagging"));
  hWSoftMuVeto      = ((TH1F*) FindOutput("hWSoftMuVeto"));
  hWeffSoftMuVeto   = ((TH1F*) FindOutput("hWeffSoftMuVeto"));
  
  // TwoLeptons level histograms                                                    
  //----------------------------------------------------------------------------                         
  
  hPtLepton1TwoLeptonsLevel       = ((TH1F*) FindOutput("hPtLepton1TwoLeptonsLevel"));
  hPtLepton2TwoLeptonsLevel       = ((TH1F*) FindOutput("hPtLepton2TwoLeptonsLevel"));
  hPtDiLeptonTwoLeptonsLevel      = ((TH1F*) FindOutput("hPtDiLeptonTwoLeptonsLevel"));
  hMinvTwoLeptonsLevel            = ((TH1F*) FindOutput("hMinvTwoLeptonsLevel"));
  hMtTwoLeptonsLevel              = ((TH1F*) FindOutput("hMtTwoLeptonsLevel"));
  hMt1TwoLeptonsLevel             = ((TH1F*) FindOutput("hMt1TwoLeptonsLevel"));
  hMt2TwoLeptonsLevel             = ((TH1F*) FindOutput("hMt2TwoLeptonsLevel"));
  hpfMetTwoLeptonsLevel           = ((TH1F*) FindOutput("hpfMetTwoLeptonsLevel"));
  hpminMetTwoLeptonsLevel         = ((TH1F*) FindOutput("hpminMetTwoLeptonsLevel"));
  hDeltaRLeptonsTwoLeptonsLevel   = ((TH1F*) FindOutput("hDeltaRLeptonsTwoLeptonsLevel"));
  hDeltaPhiLeptonsTwoLeptonsLevel = ((TH1F*) FindOutput("hDeltaPhiLeptonsTwoLeptonsLevel"));
  hDPhiPtllJetTwoLeptonsLevel     = ((TH1F*) FindOutput("hDPhiPtllJetTwoLeptonsLevel"));
  hDPhillMetTwoLeptonsLevel       = ((TH1F*) FindOutput("hDPhillMetTwoLeptonsLevel"));
  hMcTwoLeptonsLevel              = ((TH1F*) FindOutput("hMcTwoLeptonsLevel"));
  hTrkMetTwoLeptonsLevel          = ((TH1F*) FindOutput("hTrkMetTwoLeptonsLevel"));
  hHtTwoLeptonsLevel              = ((TH1F*) FindOutput("hHtTwoLeptonsLevel"));
  hdphijetjetTwoLeptonsLevel      = ((TH1F*) FindOutput("hdphijetjetTwoLeptonsLevel"));
  hdetajetjetTwoLeptonsLevel      = ((TH1F*) FindOutput("hdetajetjetTwoLeptonsLevel"));

  // WW level histograms                                      
  //----------------------------------------------------------------------------                
                                  
  hPtLepton1WWLevel       = ((TH1F*) FindOutput("hPtLepton1WWLevel"));
  hPtLepton2WWLevel       = ((TH1F*) FindOutput("hPtLepton2WWLevel"));
  hPtDiLeptonWWLevel      = ((TH1F*) FindOutput("hPtDiLeptonWWLevel"));
  hMinvWWLevel            = ((TH1F*) FindOutput("hMinvWWLevel"));
  hMtWWLevel              = ((TH1F*) FindOutput("hMtWWLevel"));
  hMt1WWLevel             = ((TH1F*) FindOutput("hMt1WWLevel"));
  hMt2WWLevel             = ((TH1F*) FindOutput("hMt2WWLevel"));
  hpfMetWWLevel           = ((TH1F*) FindOutput("hpfMetWWLevel"));
  hpminMetWWLevel         = ((TH1F*) FindOutput("hpminMetWWLevel"));
  hDeltaRLeptonsWWLevel   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevel"));
  hDeltaPhiLeptonsWWLevel = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevel"));
  hDPhiPtllJetWWLevel     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevel"));
  hDPhillMetWWLevel       = ((TH1F*) FindOutput("hDPhillMetWWLevel"));
  hMcWWLevel              = ((TH1F*) FindOutput("hMcWWLevel"));
  hTrkMetWWLevel          = ((TH1F*) FindOutput("hTrkMetWWLevel"));
  hHtWWLevel              = ((TH1F*) FindOutput("hHtWWLevel"));
  hdphijetjetWWLevel      = ((TH1F*) FindOutput("hdphijetjetWWLevel"));
  hdetajetjetWWLevel      = ((TH1F*) FindOutput("hdetajetjetWWLevel"));

  //Control Region histograms
  //----------------------------------------------------------------------------                
  
  hpfMetCR  = ((TH1F*) FindOutput("hpfMetCR"));
  htrkMetCR = ((TH1F*) FindOutput("htrkMetCR"));
  hmpMetCR  = ((TH1F*) FindOutput("hmpMetCR"));

  // monoH level histograms                                             
  //----------------------------------------------------------------------------  

  for (Int_t qq = 0; qq < 4; ++qq){
    hPtLepton1monoHLevel[qq]       = ((TH1F*) FindOutput("hPtLepton1monoHLevel%.1i",qq));
    hPtLepton2monoHLevel[qq]       = ((TH1F*) FindOutput("hPtLepton2monoHLevel%.1i",qq));
    hPtDiLeptonmonoHLevel[qq]      = ((TH1F*) FindOutput("hPtDiLeptonmonoHLevel%.1i",qq));
    hMinvmonoHLevel[qq]            = ((TH1F*) FindOutput("hMinvmonoHLevel%.1i",qq));
    hMtmonoHLevel[qq]              = ((TH1F*) FindOutput("hMtmonoHLevel%.1i",qq));
    hMt1monoHLevel[qq]             = ((TH1F*) FindOutput("hMt1monoHLevel%.1i",qq));
    hMt2monoHLevel[qq]             = ((TH1F*) FindOutput("hMt2monoHLevel%.1i",qq));
    hpfMetmonoHLevel[qq]           = ((TH1F*) FindOutput("hpfMetmonoHLevel%.1i",qq));
    hpminMetmonoHLevel[qq]         = ((TH1F*) FindOutput("hpminMetmonoHLevel%.1i",qq));
    hDeltaRLeptonsmonoHLevel[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsmonoHLevel%.1i",qq));
    hDeltaPhiLeptonsmonoHLevel[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsmonoHLevel%.1i",qq));
    hDPhiPtllJetmonoHLevel[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetmonoHLevel%.1i",qq));
    hDPhillMetmonoHLevel[qq]       = ((TH1F*) FindOutput("hDPhillMetmonoHLevel%.1i",qq));
    hPtWWmonoHLevel[qq]            = ((TH1F*) FindOutput("hPtWWmonoHLevel%.1i",qq));
    hMcmonoHLevel[qq]              = ((TH1F*) FindOutput("hMcmonoHLevel%.1i",qq));
    hTrkMetmonoHLevel[qq]          = ((TH1F*) FindOutput("hTrkMetmonoHLevel%.1i",qq));
    hHtmonoHLevel[qq]              = ((TH1F*) FindOutput("hHtmonoHLevel%.1i",qq));
    hdphijetjetmonoHLevel[qq]      = ((TH1F*) FindOutput("hdphijetjetmonoHLevel%.1i",qq));
    hdetajetjetmonoHLevel[qq]      = ((TH1F*) FindOutput("hdetajetjetmonoHLevel%.1i",qq));
  }

}

void test::Summary() {

  cout << " ---------------------------------------------------" << endl;
  cout << " " << endl;
  
  InitialiseParameters();

  cout << " Number of Events::  " << NEvents  << endl;

  double factN = 1.; 
  if (XSection > 0) factN = XSection * Luminosity / NEvents; //fractionoftotalevents;
  
  cout << " Normalization factor: " << factN << endl;
  cout << endl;
  cout << " ---------------------------------------------------" << endl;

}
