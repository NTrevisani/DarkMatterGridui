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
  GetInputParameters()->TheNamedInt("MultiplicateXS", MultiplicateXS);
  
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
  // Create histos
  //------------------------------------------------------------------------------
  
  h_n_PV = CreateH1F("h_n_PV","h_n_PV",50,0.,10.);

  // WW level histograms     
  //----------------------------------------------------------------------------  
  
  for (Int_t nC=0; nC<4; nC++) {

    hPtLepton1WWLevel[nC]       = CreateH1F(Form("hPtLepton1WWLevel%.1i", nC),       "", 1000, 0, 1000);
    hPtLepton2WWLevel[nC]       = CreateH1F(Form("hPtLepton2WWLevel%.1i", nC),       "", 1000, 0, 1000);
    hPtDiLeptonWWLevel[nC]      = CreateH1F(Form("hPtDiLeptonWWLevel%.1i", nC),      "", 1000, 0, 1000);
    hMinvWWLevel[nC]            = CreateH1F(Form("hMinvWWLevel%.1i", nC),            "", 1000, 0, 1000);
    hMtWWLevel[nC]              = CreateH1F(Form("hMtWWLevel%.1i", nC),              "", 1000, 0, 1000);
    hMt1WWLevel[nC]             = CreateH1F(Form("hMt1WWLevel%.1i", nC),             "", 1000, 0, 1000);
    hMt2WWLevel[nC]             = CreateH1F(Form("hMt2WWLevel%.1i", nC),             "", 1000, 0, 1000);
    hpfMetWWLevel[nC]           = CreateH1F(Form("hpfMetWWLevel%.1i", nC),           "", 1000, 0, 1000);
    hpminMetWWLevel[nC]         = CreateH1F(Form("hpminMetWWLevel%.1i", nC),         "", 1000, 0, 1000);
    hDeltaRLeptonsWWLevel[nC]   = CreateH1F(Form("hDeltaRLeptonsWWLevel%.1i", nC),   "",  50, 0,   5);
    hDeltaPhiLeptonsWWLevel[nC] = CreateH1F(Form("hDeltaPhiLeptonsWWLevel%.1i", nC), "",  32, 0, 3.2);
    hDPhiPtllJetWWLevel[nC]     = CreateH1F(Form("hDPhiPtllJetWWLevel%.1i", nC),     "",  32, 0, 3.2);
    hPtWWWWLevel[nC]            = CreateH1F(Form("hPtWWWWLevel%.1i", nC),                   "", 1000, 0, 1000);
    hMcWWLevel[nC]              = CreateH1F(Form("hMcWWLevel%.1i", nC),              "", 1000, 0, 1000);
    hTrkMetWWLevel[nC]          = CreateH1F(Form("hTrkMetWWLevel%.1i", nC),          "", 1000, 0, 1000);
    hHtWWLevel[nC]              = CreateH1F(Form("hHtWWLevel%.1i", nC),              "", 1000, 0, 1000);
  }
  
  // TwoLeptons level histograms   
  //----------------------------------------------------------------------------
  
  hPtLepton1TwoLeptonsLevel       = CreateH1F("hPtLepton1TwoLeptonsLevel",       "", 1000, 0, 1000);
  hPtLepton2TwoLeptonsLevel       = CreateH1F("hPtLepton2TwoLeptonsLevel",       "", 1000, 0, 1000);
  hPtDiLeptonTwoLeptonsLevel      = CreateH1F("hPtDiLeptonTwoLeptonsLevel",      "", 1000, 0, 1000);
  hMinvTwoLeptonsLevel            = CreateH1F("hMinvTwoLeptonsLevel",            "", 1000, 0, 1000);
  hMtTwoLeptonsLevel              = CreateH1F("hMtTwoLeptonsLevel",              "", 1000, 0, 1000);
  hMt1TwoLeptonsLevel             = CreateH1F("hMt1TwoLeptonsLevel",             "", 1000, 0, 1000);
  hMt2TwoLeptonsLevel             = CreateH1F("hMt2TwoLeptonsLevel",             "", 1000, 0, 1000);
  hpfMetTwoLeptonsLevel           = CreateH1F("hpfMetTwoLeptonsLevel",           "", 1000, 0, 1000);
  hpminMetTwoLeptonsLevel         = CreateH1F("hpminMetTwoLeptonsLevel",         "", 1000, 0, 1000);
  hDeltaRLeptonsTwoLeptonsLevel   = CreateH1F("hDeltaRLeptonsTwoLeptonsLevel",   "",  50, 0,   5);
  hDeltaPhiLeptonsTwoLeptonsLevel = CreateH1F("hDeltaPhiLeptonsTwoLeptonsLevel", "",  32, 0, 3.2);
  hDPhiPtllJetTwoLeptonsLevel     = CreateH1F("hDPhiPtllJetTwoLeptonsLevel",     "",  32, 0, 3.2);
  hDPhillMetTwoLeptonsLevel       = CreateH1F("hDPhillMetTwoLeptonsLevel",       "",  32, 0, 3.2);
  hMcTwoLeptonsLevel              = CreateH1F("hMcTwoLeptonsLevel",              "", 1000, 0, 1000);  
  hPtWWTwoLeptonsLevel            = CreateH1F("hPtWWTwoLeptonsLevel",            "", 1000, 0, 1000);
  hTrkMetTwoLeptonsLevel          = CreateH1F("hTrkMetTwoLeptonsLevel",          "", 1000, 0, 1000);
  hHtTwoLeptonsLevel              = CreateH1F("hHtTwoLeptonsLevel",              "", 1000, 0, 1000);
}

// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {
  
  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
  
  if (TheSample.Contains("Dark")){
    efficiencyW = MultiplicateXS * efficiencyW;
    totalW = MultiplicateXS * totalW;
  }
  
  h_n_PV -> Fill(1,efficiencyW);  
  
  Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
  
  Float_t jetbin = njet;
  
  //building mpmet
  Float_t dphimin = (min(dphilmet1,dphilmet2));
  Float_t fullpmet = 0;
  Float_t trkpmet = 0;
  
  if (dphimin < TMath::Pi() / 2)
    fullpmet = pfType1Met * sin(dphimin);
  else 
    fullpmet = pfType1Met;
  
  if (dphimin < TMath::Pi() / 2)
    trkpmet = trkMet * sin(dphimin);
  else
    trkpmet = trkMet;
  
  Float_t mpmet = min(trkpmet,fullpmet);
  
  Float_t metvar = (njet <= 1) ? mpmet : pfType1Met;
  
  //building Ht
  Float_t Ht = std_vector_lepton_pt->at(0) + std_vector_lepton_pt->at(1) + pfType1Met;
  
  if(njet > 10) njet = 10;
  for (int i = 0; i < njet; ++i)
    if(std_vector_jet_pt->at(i) > 0)
      Ht += std_vector_jet_pt->at(i);
  
  //reject background prompt-prompt events
  Float_t prompt = 1;
  if(TheSample == "TTJets" || TheSample == "Top" || TheSample == "QCD" || TheSample == "WJets")
    if(fabs(std_vector_leptonGen_mpid -> at(0)) == 24 && 
       fabs(std_vector_leptonGen_mstatus -> at(0)) > 20 && 
      fabs(std_vector_leptonGen_mstatus -> at(0)) < 30)
      if(fabs(std_vector_leptonGen_mpid -> at(1)) == 24 && 
	 fabs(std_vector_leptonGen_mstatus -> at(1)) > 20 && 
	 fabs(std_vector_leptonGen_mstatus -> at(1)) < 30)
	prompt = 0;
  

  //building ptWW
  TLorentzVector L1,L2;
  TLorentzVector MET;
  Float_t ptWW = 0;
  
  if (pt1>0 && pt2>0) {  
    L1.SetPtEtaPhiM(pt1, 0, phi1, 0.);
    L2.SetPtEtaPhiM(pt2, 0, phi2, 0.);
    MET.SetPtEtaPhiM(pfType1Met, 0, pfType1Metphi, 0.);
    ptWW = (L1+L2+MET).Pt();
  }
  
  //building Mt of the two Ws
  Float_t Mt1 = 0;
  Float_t Mt2 = 0;
   
  if( pfType1Met > 0){
    if( pt1 > 0 )
      Mt1 = sqrt(2*pt1*pfType1Met*(1 - cos(dphilmet1)));  
    if( pt2 > 0 )
      Mt2 = sqrt(2*pt2*pfType1Met*(1 - cos(dphilmet2)));  
  }
  
  //building Mc
  Float_t Mc = 0;
  
  if (ptll > 0 && mll > 0 && pfType1Met > 0)
    Mc = sqrt( pow(sqrt(ptll*ptll + mll*mll) + pfType1Met,2) - pow(ptll + pfType1Met,2) );
  
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
	     )
	  {

	    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    //
	    // Main analisis
	    //
	    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    
	    /*if( mth > 260)*/
	    /*if( mpmet > 132)*/
	    /*if ( drll < 0.5)*/ 
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
	      
	      if (nextra == 0) {
		
		if (dphiv || channel == 2 || channel == 3) {
		  
		  if ( mth > 260){
		    
		    hPtLepton1WWLevel[0]      ->Fill(pt1,       totalW);
		    hPtLepton2WWLevel[0]      ->Fill(pt2,       totalW);
		    hPtDiLeptonWWLevel[0]     ->Fill(ptll,      totalW);
		    hMinvWWLevel[0]           ->Fill(mll,       totalW);
		    hMtWWLevel[0]             ->Fill(mth,       totalW);
		    hMt1WWLevel[0]            ->Fill(Mt1,       totalW);
		    hMt2WWLevel[0]            ->Fill(Mt2,       totalW);
		    hpfMetWWLevel[0]          ->Fill(pfType1Met,totalW);
		    hpminMetWWLevel[0]        ->Fill(mpmet,     totalW);
		    hDeltaRLeptonsWWLevel[0]  ->Fill(drll,      totalW);
		    hDeltaPhiLeptonsWWLevel[0]->Fill(dphill,    totalW);
		    hDPhiPtllJetWWLevel[0]    ->Fill(dphilljet, totalW);
		    hPtWWWWLevel[0]           ->Fill(ptWW,      totalW);
		    hMcWWLevel[0]             ->Fill(Mc,        totalW);		  
		    hTrkMetWWLevel[0]         ->Fill(trkMet,    totalW);
		    hHtWWLevel[0]             ->Fill(Ht,        totalW);			
		    
		    if (drll < 0.5){
		      
		      hPtLepton1WWLevel[1]      ->Fill(pt1,       totalW);
		      hPtLepton2WWLevel[1]      ->Fill(pt2,       totalW);
		      hPtDiLeptonWWLevel[1]     ->Fill(ptll,      totalW);
		      hMinvWWLevel[1]           ->Fill(mll,       totalW);
		      hMtWWLevel[1]             ->Fill(mth,       totalW);
		      hMt1WWLevel[1]            ->Fill(Mt1,       totalW);
		      hMt2WWLevel[1]            ->Fill(Mt2,       totalW);
		      hpfMetWWLevel[1]          ->Fill(pfType1Met,totalW);
		      hpminMetWWLevel[1]        ->Fill(mpmet,     totalW);
		      hDeltaRLeptonsWWLevel[1]  ->Fill(drll,      totalW);
		      hDeltaPhiLeptonsWWLevel[1]->Fill(dphill,    totalW);
		      hDPhiPtllJetWWLevel[1]    ->Fill(dphilljet, totalW);
		      hPtWWWWLevel[1]           ->Fill(ptWW,      totalW);
		      hMcWWLevel[1]             ->Fill(Mc,        totalW);		  
		      hTrkMetWWLevel[1]         ->Fill(trkMet,    totalW);
		      hHtWWLevel[1]             ->Fill(Ht,        totalW);			
		      
		      if ( Mt2 > 145 ){
			
			hPtLepton1WWLevel[2]      ->Fill(pt1,       totalW);
			hPtLepton2WWLevel[2]      ->Fill(pt2,       totalW);
			hPtDiLeptonWWLevel[2]     ->Fill(ptll,      totalW);
			hMinvWWLevel[2]           ->Fill(mll,       totalW);
			hMtWWLevel[2]             ->Fill(mth,       totalW);
			hMt1WWLevel[2]            ->Fill(Mt1,       totalW);
			hMt2WWLevel[2]            ->Fill(Mt2,       totalW);
			hpfMetWWLevel[2]          ->Fill(pfType1Met,totalW);
			hpminMetWWLevel[2]        ->Fill(mpmet,     totalW);
			hDeltaRLeptonsWWLevel[2]  ->Fill(drll,      totalW);
			hDeltaPhiLeptonsWWLevel[2]->Fill(dphill,    totalW);
			hDPhiPtllJetWWLevel[2]    ->Fill(dphilljet, totalW);
			hPtWWWWLevel[2]           ->Fill(ptWW,      totalW);
			hMcWWLevel[2]             ->Fill(Mc,        totalW);		  
			hTrkMetWWLevel[2]         ->Fill(trkMet,    totalW);
			hHtWWLevel[2]             ->Fill(Ht,        totalW);						      
			
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

  h_n_PV = ((TH1F*) FindOutput("h_n_PV")); 
  
  // WW level histograms                                                                                                                         //----------------------------------------------------------------------------                                                                                                   
  for (Int_t qq = 0; qq < 4; ++qq){
  hPtLepton1WWLevel[qq]       = ((TH1F*) FindOutput("hPtLepton1WWLevel%.1i",qq));
  hPtLepton2WWLevel[qq]       = ((TH1F*) FindOutput("hPtLepton2WWLevel%.1i",qq));
  hPtDiLeptonWWLevel[qq]      = ((TH1F*) FindOutput("hPtDiLeptonWWLevel%.1i",qq));
  hMinvWWLevel[qq]            = ((TH1F*) FindOutput("hMinvWWLevel%.1i",qq));     // %.1i",qq 
  hMtWWLevel[qq]              = ((TH1F*) FindOutput("hMtWWLevel%.1i",qq));
  hMt1WWLevel[qq]             = ((TH1F*) FindOutput("hMt1WWLevel%.1i",qq));
  hMt2WWLevel[qq]             = ((TH1F*) FindOutput("hMt2WWLevel%.1i",qq));
  hpfMetWWLevel[qq]           = ((TH1F*) FindOutput("hpfMetWWLevel%.1i",qq));
  hpminMetWWLevel[qq]         = ((TH1F*) FindOutput("hpminMetWWLevel%.1i",qq));
  hDeltaRLeptonsWWLevel[qq]   = ((TH1F*) FindOutput("hDeltaRLeptonsWWLevel%.1i",qq));
  hDeltaPhiLeptonsWWLevel[qq] = ((TH1F*) FindOutput("hDeltaPhiLeptonsWWLevel%.1i",qq));
  hDPhiPtllJetWWLevel[qq]     = ((TH1F*) FindOutput("hDPhiPtllJetWWLevel%.1i",qq));
  hPtWWWWLevel[qq]            = ((TH1F*) FindOutput("hPtWWWWLevel%.1i",qq));
  hMcWWLevel[qq]              = ((TH1F*) FindOutput("hMcWWLevel%.1i",qq));
  hTrkMetWWLevel[qq]          = ((TH1F*) FindOutput("hTrkMetWWLevel%.1i",qq));
  hHtWWLevel[qq]              = ((TH1F*) FindOutput("hHtWWLevel%.1i",qq));
  }

  // TwoLeptons level histograms                                                                                                                 //----------------------------------------------------------------------------                                                                                                   
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
