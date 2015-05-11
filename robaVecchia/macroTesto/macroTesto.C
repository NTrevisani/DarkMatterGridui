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
	       Int_t nbins,
	       TString units = "",
	       TTree *tDark_,
	       TTree *tDark100_,
	       TTree *tZH_,
	       TTree *tHWW_,
	       TTree *tZHjonat_,
	       TString format
	       ){

  TH1F* hDark    = new TH1F("hDark",variable,nbins,left,right);
  TH1F* hDark100 = new TH1F("hDark100",variable,nbins,left,right);
  TH1F* hZH      = new TH1F("hZH",variable,nbins,left,right);
  TH1F* hHWW     = new TH1F("hHWW",variable,nbins,left,right);
  TH1F* hZHjonat = new TH1F("hZHjonat",variable,nbins,left,right);

  TString selectionCuts = "pt1>20 && pt2>20";

  tZH_      -> Draw(variable+">>hZH", + selectionCuts);

  /*  if(variable == "pfmet"){
    tDark_    -> Draw("pfType1Met>>hDark", + selectionCuts);
    tDark100_ -> Draw("pfType1Met>>hDark100", + selectionCuts);
    tHWW_     -> Draw("pfType1Met>>hHWW", + selectionCuts);
    tZHjonat_ -> Draw("pfType1Met>>hZHjonat", + selectionCuts);
  }
  else{*/
    tZHjonat_ -> Draw(variable+">>hZHjonat", + selectionCuts);
    tDark_    -> Draw(variable+">>hDark", + selectionCuts);
    tDark100_ -> Draw(variable+">>hDark100", + selectionCuts);
    tHWW_     -> Draw(variable+">>hHWW", + selectionCuts);
    //}


  //histograms normalization
  Double_t qqDark    = hDark->Integral();
  Double_t qqDark100 = hDark100->Integral();
  Double_t qqZH      = hZH->Integral();
  Double_t qqHWW     = hHWW->Integral();
  Double_t qqZHjonat = hZHjonat->Integral();

  for (int i = 0; i<nbins+1; ++i) {
    hDark    -> SetBinContent(i,hDark->GetBinContent(i)/qqDark);
    hDark100 -> SetBinContent(i,hDark100->GetBinContent(i)/qqDark100);
    hZH      -> SetBinContent(i,hZH->GetBinContent(i)/qqZH);
    hHWW     -> SetBinContent(i,hHWW->GetBinContent(i)/qqHWW);
    hZHjonat -> SetBinContent(i,hZHjonat->GetBinContent(i)/qqZHjonat);
  } 
  
  TCanvas *c1 = new TCanvas("variable","variable",600,600);
  c1->cd();
  
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetBottomMargin(0.18);
  pad1->Draw();
  pad1->cd();
  Float_t rangeY = max(hDark->GetMaximum(),hZH->GetMaximum());
  rangeY = max(rangeY,hHWW->GetMaximum());
  rangeY = max(rangeY,hDark100->GetMaximum());
  rangeY = max(rangeY,hZHjonat->GetMaximum());
  rangeY = 1.33*rangeY;
  hZH->SetTitleSize(5.);
  hZH->GetYaxis()->SetRangeUser(0.,rangeY);
  hDark->GetYaxis()->SetRangeUser(0.,rangeY);
  hDark100->GetYaxis()->SetRangeUser(0.,rangeY);

  if(units != "[]")
    hZH->GetXaxis()->SetTitle(variable+" "+units);
  else
    hZH->GetXaxis()->SetTitle(variable);
  
  hZH->GetXaxis()->SetTitleSize(0.07);
  hZH->GetXaxis()->SetTitleOffset(0.9);
  hZH->GetXaxis()->SetLabelSize(0.05);
  hZH->GetYaxis()->SetTitleSize(0.05);
  hZH->GetYaxis()->SetLabelSize(0.05);  

  hZH->SetLineWidth(3);
  hDark->SetLineWidth(3);
  hDark100->SetLineWidth(3);
  hHWW->SetLineWidth(3);
  hZHjonat->SetLineWidth(3);

  hZH->SetStats(0);
  hDark->SetStats(0);
  hDark100->SetStats(0);
  hHWW->SetStats(0);
  hZHjonat->SetStats(0);

  hDark->SetLineColor(kBlue);
  hDark100->SetLineColor(kBlack);
  hZH->SetLineColor(kRed);
  hHWW->SetLineColor(kGreen+2);
  hZHjonat->SetLineColor(kViolet);

  hZH->Draw();
  hDark100->Draw("same");
  hDark->Draw("same");
  hHWW->Draw("same");
  //hZHjonat->Draw("same");

  TLegend* leg = new TLegend(0.20,0.75,0.70,0.89);
  leg->AddEntry(hZH,"ZH #rightarrow #nu#nuWW","l");
  leg->AddEntry(hDark,"HXX #rightarrow WWXX, m_{X} = 1GeV","l");
  leg->AddEntry(hDark100,"HXX #rightarrow WWXX, m_{X} = 100GeV","l");
  leg->AddEntry(hHWW,"H #rightarrow WW","l");
  //leg->AddEntry(hZHjonat,"ZH #rightarrow #nu#nuWW - spring14","l");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->Draw();
  
  c1->Print(variable+"."+format,format);
  delete hDark;
  delete hDark100;
  delete hZH;
  delete hHWW;
  delete hZHjonat;
  delete c1;

  gSystem->Exec("mv " + variable + "." + format + " " + "distributions/" + format);
}

//main function
void macroTesto(TString printMode = "png"){
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
  TFile* fDark    = new TFile("/gpfs/csic_projects/cms/trevisanin/newLatino/latino_stepB_Higgs_hzpzp_ww_1GeV.root");
  TFile* fZH      = new TFile("/gpfs/csic_projects/cms/trevisanin/newLatino/stepB_ppTOzh_zTO2v_hTOwwTO2l2v_13TeV.root");
  TFile* fHWW     = new TFile("/gpfs/csic_projects/cms/trevisanin/newLatino/latino_stepB_HWW125.root");
  TFile* fDark100 = new TFile("/gpfs/csic_projects/cms/trevisanin/newLatino/latino_stepB_Higgs_hzpzp_ww_100GeV.root");
  TFile* fZHjonat = new TFile("/gpfs/csic_projects/cms/trevisanin/newLatino/latino_stepB_VH_HToWW_WToLNu.root");

  TTree * tDark = (TTree*) fDark -> Get("latino");
  TTree * tDark100 = (TTree*) fDark100 -> Get("latino");
  TTree * tZH = (TTree*) fZH -> Get("latino");
  TTree * tHWW = (TTree*) fHWW -> Get("latino");
  TTree * tZHjonat = (TTree*) fZHjonat -> Get("latino");

  TObjArray *tl = tDark->GetListOfBranches();
  TString nBranch = tl->First()->GetName();

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
    drawPlots(var, leftBound, rightBound, nbin, units, tDark, tDark100, tZH, tHWW, tZHjonat, printMode);
  } 
  
  inFile.close();
  gSystem->Exec("rm out.tmp");
}
