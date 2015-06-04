//compares variables value from latinos H->WW->lvlv with HXX->XXWW->XXlvlv and HZ->WWvv->lvlvvv
//run typing:  root -l -b -q 'macroHisto.C("pdf")'

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
	       TFile *fDark10_,
	       TFile *fDark1000_,	      
	       TFile *fDark500_,
	       TFile *fZH_,
	       TFile *fHWW_,
	       TFile *fWW_,
	       TString format,
	       TString drawLog,
	       TString norm
	       ){

  TH1F *hDark     = (TH1F*)fDark_->Get(variable);
  TH1F *hDark100  = (TH1F*)fDark100_->Get(variable);
  TH1F *hDark10   = (TH1F*)fDark10_->Get(variable);
  TH1F *hDark500  = (TH1F*)fDark500_->Get(variable);
  TH1F *hDark1000 = (TH1F*)fDark1000_->Get(variable);
  TH1F *hZH       = (TH1F*)fZH_->Get(variable);
  TH1F *hHWW      = (TH1F*)fHWW_->Get(variable);
  TH1F *hWW       = (TH1F*)fWW_->Get(variable);
  
  hDark     -> Rebin(nrebin);
  hDark10   -> Rebin(nrebin);
  hDark100  -> Rebin(nrebin);
  hDark500  -> Rebin(nrebin);
  hDark1000 -> Rebin(nrebin);
  hZH       -> Rebin(nrebin);
  hHWW      -> Rebin(nrebin);
  hWW       -> Rebin(nrebin);

  hDark     -> GetXaxis() -> SetRangeUser(left,right);
  hDark10   -> GetXaxis() -> SetRangeUser(left,right);
  hDark100  -> GetXaxis() -> SetRangeUser(left,right);
  hDark500  -> GetXaxis() -> SetRangeUser(left,right);
  hDark1000 -> GetXaxis() -> SetRangeUser(left,right);
  hZH       -> GetXaxis() -> SetRangeUser(left,right);
  hHWW      -> GetXaxis() -> SetRangeUser(left,right);
  hWW       -> GetXaxis() -> SetRangeUser(left,right);

  //histograms normalization
  if (norm == "on"){
    hDark     ->Scale(1./hDark->Integral());
    hDark10   ->Scale(1./hDark10->Integral());
    hDark100  ->Scale(1./hDark100->Integral());
    hDark500  ->Scale(1./hDark500->Integral());
    hDark1000 ->Scale(1./hDark1000->Integral());
    hZH       ->Scale(1./hZH->Integral());
    hHWW      ->Scale(1./hHWW->Integral());
    hWW       ->Scale(1./hWW->Integral());
  }

  TCanvas *c1 = new TCanvas("variable","variable",600,800);
  c1->cd();
  
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetBottomMargin(0.18);
  pad1->Draw();
  pad1->cd();

  Float_t rangeY = max(hDark->GetMaximum(),hZH->GetMaximum());
  rangeY = max(rangeY,hHWW->GetMaximum());
  rangeY = max(rangeY,hDark->GetMaximum());
  rangeY = max(rangeY,hDark10->GetMaximum());
  rangeY = max(rangeY,hDark100->GetMaximum());
  rangeY = max(rangeY,hDark500->GetMaximum());
  rangeY = max(rangeY,hDark1000->GetMaximum());
  rangeY = max(rangeY,hWW->GetMaximum());
  hZH->SetTitleSize(4.);

  Float_t rangeMin = 0.00000001;

  if (drawLog == "on"){
    hZH->GetYaxis()->SetRangeUser(rangeMin,rangeY*pow(10.,log10(rangeY/rangeMin))/10);
    hDark->GetYaxis()->SetRangeUser(rangeMin,rangeY*pow(10.,log10(rangeY/rangeMin))/10);
    hDark100->GetYaxis()->SetRangeUser(rangeMin,rangeY*pow(10.,log10(rangeY/rangeMin))/10);
    pad1->SetLogy();
  }  
  else if (drawLog == "off"){
    rangeY = 1.5*rangeY;
    hZH->GetYaxis()->SetRangeUser(0.,rangeY);
    hDark->GetYaxis()->SetRangeUser(0.,rangeY);
    hDark100->GetYaxis()->SetRangeUser(0.,rangeY);
  }

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
  hDark10->SetLineWidth(3);
  hDark1000->SetLineWidth(3);
  hDark500->SetLineWidth(3);
  hHWW->SetLineWidth(3);
  hWW->SetLineWidth(3);

  hZH->SetStats(0);
  hDark->SetStats(0);
  hDark100->SetStats(0);
  hDark10->SetStats(0);
  hDark1000->SetStats(0);
  hDark500->SetStats(0);
  hHWW->SetStats(0);
  hWW->SetStats(0);

  hDark->SetLineColor(kBlue);
  hDark100->SetLineColor(kBlack);
  hDark10->SetLineColor(kYellow);
  hDark1000->SetLineColor(46);
  hDark500->SetLineColor(40);
  hZH->SetLineColor(kRed);
  hHWW->SetLineColor(kGreen+2);
  hWW->SetLineColor(kGreen+4);

  hZH->Draw();
  hDark1000->Draw("same");
  hDark500->Draw("same");
  hDark100->Draw("same");
  hDark10->Draw("same");
  hDark->Draw("same");
  hHWW->Draw("same");
  hWW->Draw("same");

  TLegend* leg = new TLegend(0.20,0.70,0.70,0.89);
  leg->AddEntry(hZH,"ZH #rightarrow #nu#nuWW","l");
  leg->AddEntry(hDark,"HXX #rightarrow WWXX, m_{X} = 1GeV","l");
  leg->AddEntry(hDark100,"HXX #rightarrow WWXX, m_{X} = 100GeV","l");
  leg->AddEntry(hDark10,"HXX #rightarrow WWXX, m_{X} = 10GeV","l");
  leg->AddEntry(hDark1000,"HXX #rightarrow WWXX, m_{X} = 1000GeV","l");
  leg->AddEntry(hDark500,"HXX #rightarrow WWXX, m_{X} = 500GeV","l");
  leg->AddEntry(hHWW,"H #rightarrow WW","l");
  leg->AddEntry(hWW,"WW","l");
  leg->SetTextSize(0.03);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->Draw();
  
  c1->Print(variable+"."+format,format);
  delete hDark;
  delete hDark100;
  delete hDark10;
  delete hDark500;
  delete hDark1000;
  delete hZH;
  delete hHWW;
  delete hWW;
  delete c1;
  gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + format);
}

//main function
void macroHisto(TString printMode = "png", TString logMode = "on", TString normMode = "off"){
  if (printMode == "C" || printMode == "png" || printMode == "pdf"){
    gSystem->Exec("mkdir distributions");
    gSystem->Exec("mkdir distributions/" + printMode);
    cout<<"siiiiiiiiiiii"<<endl;
  }
  else{
    cout<<"please print a valid plot format"<<endl;
    return;
  }

  if(logMode != "on" && logMode != "off"){
    cout<<"I'd like to know if you want the log mode to be 'on' or 'off'"<<endl;
    return;
  }

  if(normMode != "on" && normMode != "off"){
    cout<<"I'd like to know if you want the normalization mode to be 'on' or 'off'"<<endl;
    return;
  }

  gSystem->Exec("mkdir distributions/" + printMode);
  TFile* fDark     = new TFile("rootFiles/AllJet/OF/Dark1.root");
  TFile* fZH       = new TFile("rootFiles/AllJet/OF/ZH.root");
  TFile* fHWW      = new TFile("rootFiles/AllJet/OF/HWW.root");
  TFile* fDark100  = new TFile("rootFiles/AllJet/OF/Dark100.root");
  TFile* fDark10   = new TFile("rootFiles/AllJet/OF/Dark10.root");
  TFile* fDark1000 = new TFile("rootFiles/AllJet/OF/Dark1000.root");
  TFile* fDark500  = new TFile("rootFiles/AllJet/OF/Dark500.root");
  TFile* fWW       = new TFile("rootFiles/AllJet/OF/WW.root");

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
    drawPlots(var, leftBound, rightBound, nbin, units, fDark, fDark100, fDark10, fDark1000, fDark500, fZH, fHWW, fWW, printMode,logMode,normMode);
  } 
    
  inFile.close();
  gSystem->Exec("rm out.tmp");
}
