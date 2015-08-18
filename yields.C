//root -l -q 'yields.C("MuMu")'

#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TMath.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iomanip>

std::ofstream inFile("yields.txt",std::ios::out);

const int nProcesses = 4;

enum {iWW, iHWW, iZH, iDark10};

TFile *input[nProcesses];
TH1F  *histo[nProcesses];

TString process[nProcesses];

process[iWW]     = "WW";
process[iHWW]    = "HWW";
process[iZH]     = "ZH";
process[iDark10] = "Dark10";

Float_t yield[nProcesses];

const int nCuts = 3;

//defining cut levels
TString cutLevel[nCuts];
cutLevel[0] = "WW Level";
cutLevel[1] = "mono Higgs Level";
cutLevel[2] = "Control Region";

TString cutHisto[nCuts];
cutHisto[0] = "hpfMetWWLevel";
cutHisto[1] = "hpfMetmonoHLevel2";
cutHisto[2] = "hpfMetCR";
/*
TString cutLevel[nCuts];
cutLevel[0] = "Two Leptons Level";
cutLevel[1] = "WW Level";
cutLevel[2] = "mono Higgs Level";

TString cutHisto[nCuts];
cutHisto[0] = "hDeltaPhiLeptonsTwoLeptonsLevel";
cutHisto[1] = "hDeltaPhiLeptonsWWLevel";
cutHisto[2] = "hDeltaPhiLeptonsmonoHLevel1";
*/
void yields(TString flavourChannel = "",
	    TString muonID         = ""
	    ){

  if( muonID != "MediumID" && muonID != "TightID" && muonID != "TightIDTighterIP")  
  if (flavourChannel != "All" && flavourChannel !=  "OF" && flavourChannel !=  "SF" && flavourChannel !=  "MuMu" && flavourChannel !=  "EE" && flavourChannel !=  "EMu" && flavourChannel !=  "MuE"){

    cout<<"**************************************************************************************************************"<<endl;
    cout<<"Please select a valid flavour channel to look at: 'All' or 'OF' or 'SF' or 'MuMu' or 'EE' or 'EMu' or 'MuE'..."<<endl;
    cout<<"...And a valid muonID: 'MediumID', 'TightID' or 'TightIDTighterIP'"<<endl;
    cout<<"For example: root -l -q 'yields.C(\"OF\",\"MediumIDTighterIP\")'"<<endl;
    cout<<"**************************************************************************************************************"<<endl;
    return;
  }
  /*
  if (flavourChannel == "OF" || flavourChannel == "EMu" || flavourChannel == "MuE")
    cutLevel[7] = cutLevel[7] + " 30~GeV";

  else if (flavourChannel == "SF" || flavourChannel == "MuMu" || flavourChannel == "EE")
    cutLevel[7] = cutLevel[7] + " 45~GeV";

  else if (flavourChannel == "All"){
    cutLevel[6] = cutLevel[6] + " (SF only)";
    cutLevel[7] = "p$_T^{ll}$ > 30~GeV (45~GeV) in OF (SF)";
  }
  */
  TString path = "rootFiles/AllJet/" + flavourChannel + "/";
  
  inFile<<"\\begin{tabular}{cSSS}"<<endl;
  cout<<"\\begin{tabular}{cSSS}"<<endl;
  inFile<<"\\toprule"<<endl;
  cout<<"\\toprule"<<endl;
  inFile<<"Cut Level";
  cout<<"Cut Level";
  for (int q = 0; q < nProcesses; ++q){
    inFile<<" & {"<<process[q]<<"}";
    cout<<" & {"<<process[q]<<"}";
  }			
  inFile<<" \\\\ "<<endl;	
  cout<<" \\\\ "<<endl;
  inFile<<"\\midrule"<<endl;
  cout<<"\\midrule"<<endl;

  for(int c = 0; c < nCuts; ++c){

    if (flavourChannel == "OF" || flavourChannel == "EMu" || flavourChannel == "MuE") 
      if (c == 6 || c == 4) continue;
    
    //signal and background files and histograms
    for (int ip = 0; ip < nProcesses; ++ip){
      TString fileName = path + process[ip] + ".root";
      input[ip] = new TFile(fileName, "read");
      
      if (input[ip] -> GetListOfKeys() -> Contains(cutHisto[c]) == 0){
	cout<<"I cannot find the histogram named "<<cutHisto[c]<<". I'll skip it."<<endl;
	return 0;
      }
      
      histo[ip]  = (TH1F*) input[ip] -> Get(cutHisto[c]);
    }
    
    //begin the table
    inFile<<cutLevel[c];
    cout<<cutLevel[c];
    //number of events
    for(int i = 0; i < nProcesses; ++i){
      
      yield[i] = /*TMath::Nint*/(histo[i] -> Integral()); 
      std::cout << std::fixed;
      inFile<<" & "<<setprecision(2)<<yield[i];
      cout<<" & "<<setprecision(2)<<yield[i];
    }
    inFile<<"\\\\"<<endl;
    cout<<"\\\\"<<endl;
  }
  
  //end of the table
  inFile<<"\\bottomrule"<<endl;
  cout<<"\\bottomrule"<<endl;
  //inFile<<"Percentage"<<" & "<<DarkNoCuts<<" & "<<ZHNoCuts<<" & "<<HWWNoCuts<<" & "<<WWNoCuts<<"\\"<<endl;
  inFile<<"\\end{tabular}"<<endl;
  cout<<"\\end{tabular}"<<endl;
  inFile.close();

  gSystem -> Exec("mv yields.txt distributions/" + flavourChannel + "/" + muonID);
}
