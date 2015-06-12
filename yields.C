//root -l -q 'yields.C("Dark1","1")'

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
#include <stdio.h>
#include <stdlib.h>
#include <string>

void yields(TString nameF = "Dark1", TString scaleXS = "1"){

  //initial variables
  TString variable = "";
  TH1F *hbkg[3];

  //signal files
  TFile *fsig = new TFile("rootFiles/AllJet/OF/" + scaleXS + "/" + nameF + ".root","read");

  //defining background names
  TString name[3];
  name[0] = "ZH"; 
  name[1] = "HWW";  
  name[2] = "WW";

  //defining cut level names
  TString cutLevel[3];
  cutLevel[0] = "$m_T Higgs$ > 260~GeV"; 
  cutLevel[1] = "\\Delta R_{ll} < 0.5";  
  cutLevel[2] = "M$_T$(W$_2$) > 145~GeV";

  //background files  
  TFile *fbkg[3];
  fbkg[0] = new TFile("rootFiles/AllJet/OF/" + scaleXS + "/"  + name[0]  + ".root","read");    
  fbkg[1] = new TFile("rootFiles/AllJet/OF/" + scaleXS + "/"  + name[1]  + ".root","read");
  fbkg[2] = new TFile("rootFiles/AllJet/OF/" + scaleXS + "/"  + name[2]  + ".root","read"); 

  //begin the table
  std::ofstream inFile("yields.txt",std::ios::out);

  inFile<<"\\begin{tabular}{cSSSS}"<<endl;
  inFile<<"\\toprule"<<endl;
  inFile<<"Cut Level & {Dark Matter} & {ZH} & {HWW} & {WW}\\"<<endl;
  inFile<<"\\midrule"<<endl;

  //******************************************
  variable = "hPtLepton1TwoLeptonsLevel";

  //defining signal histogram (no cuts applied)
  TH1F* hsig = (TH1F*)fsig -> Get(variable);

  //defining background histograms (no cuts applied)
  hbkg[0]  = (TH1F*)fbkg[0] -> Get(variable);
  hbkg[1]  = (TH1F*)fbkg[1] -> Get(variable);
  hbkg[2]  = (TH1F*)fbkg[2] -> Get(variable);
  
  //number of events (no cuts applied)
  Float_t DarkNoCuts = hsig    -> Integral();
  Float_t ZHNoCuts   = hbkg[0] -> Integral();
  Float_t HWWNoCuts  = hbkg[1] -> Integral();
  Float_t WWNoCuts   = hbkg[2] -> Integral();

  //writing the table (no cuts line)
  inFile<<"None & "<<DarkNoCuts<<" & "<<ZHNoCuts<<" & "<<HWWNoCuts<<" & "<<WWNoCuts<<"\\"<<endl;

  //******************************************
  variable = "hPtLepton1WWLevel";
  //std::string buffer;
  char buffer[10];

  for (int i = 0; i < 3; ++i){
    //defining signal histogram (subsequential cuts)
    sprintf(buffer,"%d",i);
    TString addNumber = buffer;
    hsig = (TH1F*)fsig -> Get(variable + addNumber);

    //defining background histograms (subsequential cuts)
    hbkg[0]  = (TH1F*)fbkg[0] -> Get(variable + addNumber);
    hbkg[1]  = (TH1F*)fbkg[1] -> Get(variable + addNumber);
    hbkg[2]  = (TH1F*)fbkg[2] -> Get(variable + addNumber);
    
    //number of events (subsequential cuts)
    Float_t DarkNoCuts = hsig    -> Integral();
    Float_t ZHNoCuts   = hbkg[0] -> Integral();
    Float_t HWWNoCuts  = hbkg[1] -> Integral();
    Float_t WWNoCuts   = hbkg[2] -> Integral();
    
    //writing the table (subsequential cuts)
    inFile<<cutLevel[i]<<" & "<<DarkNoCuts<<" & "<<ZHNoCuts<<" & "<<HWWNoCuts<<" & "<<WWNoCuts<<"\\"<<endl;
  }
  //******************************************
  /*
  variable = "hPtLepton1WWLevel3";

  //defining signal histogram (cuts with correlation)
  hsig = (TH1F*)fsig -> Get(variable);

  //defining background histograms (cuts with correlation)
  hbkg[0]  = (TH1F*)fbkg[0] -> Get(variable);
  hbkg[1]  = (TH1F*)fbkg[1] -> Get(variable);
  hbkg[2]  = (TH1F*)fbkg[2] -> Get(variable);
  
  //number of events (cuts with correlation)
  Float_t DarkNoCuts = hsig    -> Integral();
  Float_t ZHNoCuts   = hbkg[0] -> Integral();
  Float_t HWWNoCuts  = hbkg[1] -> Integral();
  Float_t WWNoCuts   = hbkg[2] -> Integral();

  //writing the table (cuts with correlation)
  inFile<<"Correlation & "<<DarkNoCuts<<" & "<<ZHNoCuts<<" & "<<HWWNoCuts<<" & "<<WWNoCuts<<"\\"<<endl;
  */
  //end of the table
  inFile<<"\\bottomrule"<<endl;
  inFile<<"\\end{tabular}"<<endl;

  inFile.close();

  //cout<<"Data = "<<hdata -> Integral()<<" -> "<<100.*hsig1 -> Integral()/total<<"%"<<endl;
  //cout<<"total = "<<total<<endl;

  //  return 0;//check;
}
