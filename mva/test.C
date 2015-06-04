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
  // Create tree and branches
  //------------------------------------------------------------------------------



  tree = CreateTree("nt","nt");

  tree->Branch("hPtLepton1",&hPtLepton1,"hPtLepton1");
  tree->Branch("hPtLepton2",&hPtLepton2,"hPtLepton2");
  tree->Branch("hPtDiLepton",&hPtDiLepton,"hPtDiLepton");
  tree->Branch("hMinv",&hMinv,"hMinv");
  tree->Branch("hMth",&hMth,"hMth");
  tree->Branch("hMt1",&hMt1,"hMt1");
  tree->Branch("hMt2",&hMt2,"hMt2");
  tree->Branch("hNJets30",&hNJets30,"hNJets30");
  tree->Branch("hpfMet",&hpfMet,"hpfMet");
  tree->Branch("hmpMet",&hmpMet,"hmpMet");
  tree->Branch("hDeltaRLeptons",&hDeltaRLeptons,"hDeltaRLeptons");
  tree->Branch("hDeltaPhiLeptons",&hDeltaPhiLeptons,"hDeltaPhiLeptons");
  tree->Branch("hDPhiPtllMet",&hDPhiPtllMet,"hDPhiPtllMet");
  tree->Branch("hDPhiPtllJet",&hDPhiPtllJet,"hDPhiPtllJet");
  tree->Branch("hPtWW",&hPtWW,"hPtWW");
  tree->Branch("hMc",&hMc,"hMc");
  tree->Branch("hHt",&hHt,"hHt");
  tree->Branch("hTrkMet",&hTrkMet,"hTrkMet");
  tree->Branch("htotalW",&htotalW,"htotalW");
}


// The InsideLoop() function is called for each entry in the tree to be processed  
void test::InsideLoop() {

  TString TheSample = GetInputParameters()->TheNamedString("theSample");
  
  Double_t efficiencyW = effW * triggW;
  Double_t totalW      = -999;
  
  efficiencyW = puW * effW * triggW ;
  totalW      = (1 + 0.5 * (dataset >= 82 && dataset <= 84)) * baseW * efficiencyW * Luminosity;
  
  cout<<"ciao"<<endl;
  
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
  
   Float_t mpmet = min(trkpmet,fullpmet);

   Float_t metvar = (njet <= 1) ? mpmet : pfType1Met;
  

   //Defining Isolation
   //--------------------------------------------------------------------------

   Float_t isoOne = 1.;
   Float_t isoZero = 1.;
   
   if(std_vector_lepton_pt->size() > 1){
     if(std_vector_lepton_neutralHadronIso->at(0) + std_vector_lepton_photonIso->at(0) - 0.5 * std_vector_lepton_sumPUPt->at(0) > 0)
       isoZero = 
	 std_vector_lepton_chargedHadronIso->at(0) + 
	 std_vector_lepton_neutralHadronIso->at(0) + 
	 std_vector_lepton_photonIso->at(0) - 
	 0.5 * std_vector_lepton_sumPUPt->at(0);
     else
       isoZero = std_vector_lepton_chargedHadronIso->at(0);
     isoZero = isoZero / std_vector_lepton_pt->at(0);
     cout<<"isoZero: "<<isoZero<<endl;
     
     if(std_vector_lepton_neutralHadronIso->at(1) + std_vector_lepton_photonIso->at(1) - 0.5 * std_vector_lepton_sumPUPt->at(1) > 0)
       isoOne = 
	 std_vector_lepton_chargedHadronIso->at(1) + 
	 std_vector_lepton_neutralHadronIso->at(1) +
	 std_vector_lepton_photonIso->at(1) -
	 0.5 * std_vector_lepton_sumPUPt->at(1);
     else
       isoOne = std_vector_lepton_chargedHadronIso->at(1);
     isoOne = isoOne / std_vector_lepton_pt->at(1);
     cout<<"isoOne :"<<isoOne<<endl;
     cout<<" "<<endl;
   }
   /*   
   if(fabs(std_vector_lepton_id->at(0)) == 13){
     //hIsoMu->Fill(isoZero);
     hchHadronMu->Fill(std_vector_lepton_chargedHadronIso->at(0));
     hneuHadronMu->Fill(std_vector_lepton_neutralHadronIso->at(0));
     hphotonMu->Fill(std_vector_lepton_photonIso->at(0));
     hPUMu->Fill(std_vector_lepton_sumPUPt->at(0));
   }

   if(fabs(std_vector_lepton_id->at(1)) == 13){
     hIsoMu->Fill(isoOne);
     hchHadronMu->Fill(std_vector_lepton_chargedHadronIso->at(1));
     hneuHadronMu->Fill(std_vector_lepton_neutralHadronIso->at(1));
     hphotonMu->Fill(std_vector_lepton_photonIso->at(1));
     hPUMu->Fill(std_vector_lepton_sumPUPt->at(1));
   }
   
   if(fabs(std_vector_lepton_id->at(0)) == 11){
     hIsoEl->Fill(isoZero);
     hchHadronEl->Fill(std_vector_lepton_chargedHadronIso->at(0));
     hneuHadronEl->Fill(std_vector_lepton_neutralHadronIso->at(0));
     hphotonEl->Fill(std_vector_lepton_photonIso->at(0));
     hPUEl->Fill(std_vector_lepton_sumPUPt->at(0));
   }

   if(fabs(std_vector_lepton_id->at(1)) == 11){
     hIsoEl->Fill(isoOne);
     hchHadronEl->Fill(std_vector_lepton_chargedHadronIso->at(1));
     hneuHadronEl->Fill(std_vector_lepton_neutralHadronIso->at(1));
     hphotonEl->Fill(std_vector_lepton_photonIso->at(1));
     hPUEl->Fill(std_vector_lepton_sumPUPt->at(1));
   }
   */
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

   //building Ht
   Float_t Ht = 0;
   if(pt1 > 0 && pt2 > 0 && pfType1Met > 0 )
     std_vector_lepton_pt->at(0) + std_vector_lepton_pt->at(1) + pfType1Met;

   if(njet > 10) njet = 10;
   for (int i = 0; i < njet; ++i)
     if(std_vector_jet_pt->at(i) > 0)
       Ht += std_vector_jet_pt->at(i);

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
	     
	     
	     //	     if (nextra == 0) {
	       
	     //if (pfType1Met > 20 ) { 
		 
	     //	 if (mll > 12 && mll < 80) {
		   
	     //	   if (mpmet > 80){
		     
	     //	     if (dphiv || channel == 2 || channel == 3) {
		       
	     //	       if ( ptll > 70 ) {
			 
	     		 if ( mth > 0 && mth < 10000){
			   
	     //		   if ( dphill < 3.14159265 / 2.){
			     
	     //		     if(Ht < 250){
			       
			       hPtLepton1       = pt1;       
			       hPtLepton2       = pt2;       
			       hPtDiLepton      = ptll;      
			       hMinv            = mll;       
			       hMth             = mth;       
			       hNJets30         = jetbin;    
			       hpfMet           = pfType1Met;
			       hmpMet           = mpmet;     
			       hDeltaRLeptons   = drll;      
			       hDeltaPhiLeptons = dphill;    
			       hDPhiPtllJet     = dphilljet; 
			       hDPhiPtllMet     = dphillmet; 
			       hPtWW            = ptWW;      
			       hMc              = Mc;        		  
			       hTrkMet          = trkMet;    
			       hMt1             = Mt1;
			       hMt2             = Mt2;
			       hHt              = Ht;
			       htotalW          = totalW;
			       //}
			       // }
			 }  					
			       // }
			       //}
			       //}
			       //}
			       //}
			       //}
	   }
   
   	
   tree->Fill();		        
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

  tree = ((TTree*) FindOutput("nt"));

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
