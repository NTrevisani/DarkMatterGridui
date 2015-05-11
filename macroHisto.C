//compares variables value from latinos H->WW->lvlv with HXX->XXWW->XXlvlv and HZ->WWvv->lvlvvv
//run typing:  root -l -b -q 'macroTesto.C("pdf")'

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h
#include "TString.h"
#include "TLegend.h"
#include "TObjArray.h"
#include <algorithm>
#include "TCut.h"
#include "TStyle.h"
#include <iostream.h>
#include "latinoTree.h"

TString latinoVar = "";
TString darkVar = "";
Float_t range = 1000.;

//drawing instructions
void drawPlots(TString variable,
	       Float_t left, 
	       Float_t right, 
	       Int_t nrebin,
	       TString units = "",
	       TFile *fDark_,
	       TFile *fDark100_,
	       TFile *fZH_,
	       TFile *fHWW_,
	       TString format
	       ){

  TH1F *hDark    = (TH1F*)fDark_->Get(variable);
  TH1F *hDark100 = (TH1F*)fDark100_->Get(variable);
  TH1F *hZH      = (TH1F*)fZH_->Get(variable);
  TH1F *hHWW     = (TH1F*)fHWW_->Get(variable);
  
  hDark    -> Rebin(nrebin);
  hDark100 -> Rebin(nrebin);
  hZH      -> Rebin(nrebin);
  hHWW     -> Rebin(nrebin);

  hDark    -> GetXaxis() -> SetRangeUser(left,right);
  hDark100 -> GetXaxis() -> SetRangeUser(left,right);
  hZH      -> GetXaxis() -> SetRangeUser(left,right);
  hHWW     -> GetXaxis() -> SetRangeUser(left,right);

  /*
  //histograms normalization
  hDark   ->Scale(1./hDark->Integral());
  hDark100->Scale(1./hDark100->Integral());
  hZH     ->Scale(1./hZH->Integral());
  hHWW    ->Scale(1./hHWW->Integral());
  */
  TCanvas *c1 = new TCanvas("variable","variable",600,600);
  c1->cd();
  
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetBottomMargin(0.18);
  pad1->Draw();
  pad1->cd();
  Float_t rangeY = max(hDark->GetMaximum(),hZH->GetMaximum());
  rangeY = max(rangeY,hHWW->GetMaximum());
  rangeY = 1.33*rangeY;
  hZH->SetTitleSize(5.);
  hZH->GetYaxis()->SetRangeUser(0.,rangeY);
  hDark->GetYaxis()->SetRangeUser(0.,rangeY);
  hDark100->GetYaxis()->SetRangeUser(0.,rangeY);

  if(units != "[]")
    hZH->GetXaxis()->SetTitle(units);
  //else
  //hZH->GetXaxis()->SetTitle(variable);
  
  hZH->GetXaxis()->SetTitleSize(0.07);
  hZH->GetXaxis()->SetTitleOffset(0.9);
  hZH->GetXaxis()->SetLabelSize(0.05);
  hZH->GetYaxis()->SetTitleSize(0.05);
  hZH->GetYaxis()->SetLabelSize(0.05);  

  hZH->SetLineWidth(3);
  hDark->SetLineWidth(3);
  hDark100->SetLineWidth(3);
  hHWW->SetLineWidth(3);

  hZH->SetStats(0);
  hDark->SetStats(0);
  hDark100->SetStats(0);
  hHWW->SetStats(0);

  hDark->SetLineColor(kBlue);
  hDark100->SetLineColor(kBlack);
  hZH->SetLineColor(kRed);
  hHWW->SetLineColor(kGreen+2);

  hZH->Draw();
  hDark100->Draw("same");
  hDark->Draw("same");
  hHWW->Draw("same");

  TLegend* leg = new TLegend(0.20,0.75,0.70,0.89);
  leg->AddEntry(hZH,"ZH #rightarrow #nu#nuWW","l");
  leg->AddEntry(hDark,"HXX #rightarrow WWXX, m_{X} = 1GeV","l");
  leg->AddEntry(hDark100,"HXX #rightarrow WWXX, m_{X} = 100GeV","l");
  leg->AddEntry(hHWW,"H #rightarrow WW","l");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->Draw();
  
  c1->Print(variable+"."+format,format);
  delete hDark;
  delete hDark100;
  delete hZH;
  delete hHWW;
  delete c1;
  gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + format);
}

//main function
void macroHisto(TString printMode = "png"){
  if (printMode == "C" || printMode == "png" || printMode == "pdf"){
    gSystem->Exec("mkdir distributions");
    gSystem->Exec("mkdir distributions/" + printMode);
    cout<<"siiiiiiiiiiii"<<endl;
  }
  else{
    cout<<"please print a valid plot format"<<endl;
    return;
  }

  gSystem->Exec("mkdir distributions/" + printMode);
  TFile* fDark    = new TFile("rootFiles/AllJet/OF/Dark1.root");
  TFile* fZH      = new TFile("rootFiles/AllJet/OF/ZH.root");
  TFile* fHWW     = new TFile("rootFiles/AllJet/OF/HWW.root");
  TFile* fDark100 = new TFile("rootFiles/AllJet/OF/Dark100.root");

  Int_t cont = 0;
  TString var;
  Float_t leftBound = 0;
  Float_t rightBound = 0;
  TString units = "";
  Int_t nbin = 10;
  ifstream inFile("input.txt");
  std::string line;
  
  while (getline(inFile,line)){
    
    std::ofstream outFile("out.tmp",std::ios::out);
    outFile<<line<<endl;
    outFile.close();
    std::ifstream input;
    input.open("out.tmp",std::ios::in);
    input >> var >> leftBound >> rightBound >> nbin >> units;
    input.close();
    cout<<var<<" "<<leftBound<<" "<<rightBound<<" "<<nbin<<" "<<units<<" "<<printMode<<endl;
    drawPlots(var, leftBound, rightBound, nbin, units, fDark, fDark100, fZH, fHWW, printMode);
  } 
    
  inFile.close();
  gSystem->Exec("rm out.tmp");
}
