#include "TROOT.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "math.h"

float GetXofTheMax(TGraph* grafico){
  float value = 0.;
  float xmax = 0.;
  for (int  i = 0; i < 3000; ++i){
    if (grafico->Eval(i) > value){
      value = grafico->Eval(i);
      xmax = i;
    }
  }
  return xmax;
}

float GetMax(TGraph* grafico){
  float value = 0.;
  for (int  i = 0; i < 3000; ++i){
    if (grafico->Eval(i) > value){
      value = grafico->Eval(i);
    }
  }
  return value;
}

float FitBestCut(TGraph* grafico){
  float mean = GetXofTheMax(grafico);
  TF1* curva = new TF1("curva","pol5",0.,3000.);
  grafico->Fit("curva","","",mean/2,1.5*mean);
  float max = -9999;
  float maxX = -9999;
  float passo = mean/10000.;
  for (int i = 0; i < 10000; ++i){
    if(curva->Eval(mean/2 + passo*i) > max){
      max = curva->Eval(mean/2 + passo*i);
      maxX = mean/2 + passo*i;
    }
    //cout<<"max = "<<max<<", maxX = "<<maxX<<", i = "<<i<<", x = "<<mean/2 + passo*i<<", y = "<<curva->Eval(mean/2 + passo*i)<<endl;
  }
  return maxX;//curva->GetMaximum(mean/2,1.5*mean);
}

int doTheROC(TString variable = "hTrkMetTwoLeptonsLevel", 
	     TString sampleName = "Dark1", 
	     TString typeOfCut = ">",
	     Float_t max = 3000.,
	     TString printFormat = "pdf",
	     TString text = "Tracker MET [GeV]"
	     ){

  gSystem->Exec("mkdir orthogonalCuts");
  gSystem->Exec("mkdir orthogonalCuts/" + sampleName);

  if(typeOfCut != ">" && typeOfCut != "<"){
    cout<<"Please tell me if you want to keep the events above the optimal value I'll find ('>') or below it ('<')"<<endl;
    return 0;
  }

  if(printFormat != "png" && printFormat != "pdf" && printFormat != "C"){
    cout<<"Plese select a valid output format for the plots: pdf, png or C"<<endl;
    return 0;
  }

  //signal files
  TFile *fsig = new TFile("rootFiles/AllJet/OF/" + sampleName  + ".root","read");

  //defining signal histogram
  TH1F* hsig = (TH1F*)fsig -> Get(variable);

  //defining background names
  TString name[3];
  name[0] = "ZH"; 
  name[1] = "HWW";  
  name[2] = "WW";

  //background files
  TFile *fbkg[3];
  fbkg[0] = new TFile("rootFiles/AllJet/OF/" + name[0]  +".root","read");
  fbkg[1] = new TFile("rootFiles/AllJet/OF/" + name[1]  +".root","read");
  fbkg[2] = new TFile("rootFiles/AllJet/OF/" + name[2]  +".root","read");

  //defining background histograms
  TH1F *hbkg[3];
  hbkg[0]  = (TH1F*)fbkg[0] -> Get(variable);
  hbkg[1]  = (TH1F*)fbkg[1] -> Get(variable);
  hbkg[2]  = (TH1F*)fbkg[2] -> Get(variable);
  
  //defining graphs
  TGraph *g  = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();

  //total number of signal events (no cuts applied)
  Float_t totSig = hsig->Integral();

  //total number of background events (no cuts applied)
  Float_t totBkg = 0;
  for(int q = 2; q < 3; ++q)
    totBkg = totBkg + hbkg[q] ->Integral();

  //defining variables
  Float_t effSig = 0.;
  Float_t effBkg = 0.;

  //main loop: scanning efficiencies
  for (Int_t i = 0; i < 3000; ++i){
    
    //calculating signal efficiency
    Float_t sum = 0.;
    if (typeOfCut == "<")
      sum = hsig->Integral(0.,i);                               //to use when the selection will be < than the value. ex: drll
    else if (typeOfCut == ">")
      sum = hsig->Integral(3000. - i, 3000.);                     //to use when the selection will be > than the value. ex: trkMet
    effSig = sum / totSig;
    
    //calculating background efficiency (and rejection)
    Float_t bkgInt = 0;
    for(int pp = 2; pp < 3; ++pp)
      if (typeOfCut == "<")
	bkgInt = bkgInt + hbkg[pp]->Integral(0.,i);             //to use when the selection will be < than the value. ex: drll
      else if (typeOfCut == ">")      
	bkgInt = bkgInt + hbkg[pp]->Integral(3000. - i, 3000.);   //to use when the selection will be > than the value. ex: trkMet

    effBkg = bkgInt / totBkg; 
    
    //filling graphs
    g ->SetPoint(i, effSig, 1 - effBkg);

    if (typeOfCut == "<"){
      Float_t den = (bkgInt + hsig->Integral(0.,i));             //to use when the selection will be < than the value. ex: drll
      if (den != 0)
	g2->SetPoint(i, i, (hsig->Integral(0.,i)) / sqrt(den));
      g3->SetPoint(i, i, effSig);                              //to use when the selection will be < than the value. ex: drll
      g4->SetPoint(i, i, 1 - effBkg); 
    }   
 
    else if (typeOfCut == ">"){
      Float_t den = (bkgInt + hsig->Integral(3000. - i, 3000.));   //to use when the selection will be > than the value. ex: trkMet
      if (den != 0)
	g2->SetPoint(i, 3000. - i, (hsig->Integral(3000. - i, 3000.)) / sqrt(den));
      g3->SetPoint(i, 3000. - i, effSig);                        //to use when the selection will be > than the value. ex: trkMet
      g4->SetPoint(i, 3000. - i, 1. - effBkg);
    }
  }

  TCanvas *cb = new TCanvas("cb","cb",600.,600.);
  cb -> cd();

  TPad *padb = new TPad("padb","padb",0.,0.,1.,1.);
  padb->SetLeftMargin(0.15);
  padb->SetRightMargin(0.15);
  padb->SetBottomMargin(0.15);
  padb->Draw();
  padb->cd();

  hsig -> Rebin(10);
  hsig -> Scale (1. / hsig -> Integral());
  hsig -> GetYaxis()->SetRangeUser(0.,0.0180);
  hsig -> SetLineWidth(5);
  hsig -> SetStats(0);

  hsig -> GetXaxis()->SetTitle(text);
  hsig -> SetTitle("Normalized Distributions for the Sum of Leptons p_{T}");
  hsig->GetXaxis()->SetTitleSize(0.05);
  hsig->GetXaxis()->SetTitleOffset(1.2);
  hsig->GetXaxis()->SetLabelSize(0.05);
  hsig->GetYaxis()->SetTitleOffset(1.4);
  hsig->GetYaxis()->SetTitleSize(0.05);
  hsig->GetYaxis()->SetLabelSize(0.05);
  hsig->GetXaxis()->SetNdivisions(5,10,20);

  hsig -> Draw();

  TLegend* legb = new TLegend(0.60,0.55,0.75,0.85);
  legb->SetTextSize(0.04);
  legb->SetFillColor(kWhite);
  legb->SetLineColor(kWhite);

  legb -> AddEntry(hsig,"signal","l");

  for(int b = 0; b < 3; ++b){
    hbkg[b] -> Rebin(10);
    hbkg[b] -> SetStats(0);
    hbkg[b] -> SetLineWidth(5);
    hbkg[b] -> Scale (1. / hbkg[b] -> Integral());
    hbkg[b] -> GetYaxis()->SetRangeUser(0.,0.018);
    hbkg[b] -> SetLineColor(b + 1);
    legb -> AddEntry(hbkg[b], name[b],"l");
    cout<<hbkg[b]->GetTitle()<<endl;
    hbkg[b] -> Draw("same");
  }
  legb->Draw();

  g->SetTitle("ROC");
  g->GetXaxis()->SetTitle("Signal Efficiency");
  g->GetYaxis()->SetTitle("Background Rejection");
  g->SetLineWidth(5);  
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetTitleOffset(1.2);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);

  TGraph *p = new TGraph();
  p->SetPoint(0.,0.84,0.55);
  p->SetMarkerColor(kRed);
  p->SetMarkerStyle(20);
  p->SetMarkerSize(2);
  
  g2->SetTitle("Significance");
  g2->GetXaxis()->SetTitle(text);
  g2->SetLineWidth(5);
  g2->GetYaxis()->SetRangeUser(0.,1.5*GetMax(g2));//0.018);
  g2->GetXaxis()->SetRangeUser(0.,max);
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetXaxis()->SetTitleOffset(1.2);
  g2->GetXaxis()->SetLabelSize(0.05);
  g2->GetYaxis()->SetTitle("S / #sqrt{S + B}");
  g2->GetYaxis()->SetTitleOffset(2.0);
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->SetLabelSize(0.05);

  Float_t drawCut = FitBestCut(g2);
  cout<<"I'd like to cut here: "<<TMath::Nint(drawCut)<<endl;

  TLine *l = new TLine(drawCut,0.,drawCut,1.5*GetMax(g2));
  l->SetLineColor(kViolet);
  l->SetLineWidth(5);

  g3->SetLineColor(kRed);
  g3->SetLineWidth(5);
  g3->SetTitle("Signal Efficiency and Background Rejection");
  g3->GetXaxis()->SetTitle(text);
  g3->GetYaxis()->SetTitle("Signal Efficiency");
  g3->GetXaxis()->SetRangeUser(0.,max);  
  g3->GetYaxis()->SetRangeUser(0.,1.);
  g3->GetXaxis()->SetLabelOffset(0.02);
  g3->GetXaxis()->SetTitleSize(0.05);
  g3->GetXaxis()->SetTitleOffset(1.3);
  g3->GetXaxis()->SetLabelSize(0.05);
  g3->GetYaxis()->SetTitleOffset(1.4);
  g3->GetYaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->SetLabelSize(0.05);
  g4->SetLineColor(kBlack);
  g4->SetLineWidth(5);

  TLegend* leg = new TLegend(0.60,0.50,0.75,0.70);
  leg->AddEntry(g3,"Sig Eff","l");
  leg->AddEntry(g4,"Bkg Rej","l");
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);

  TLine *ll = new TLine(drawCut,0.,drawCut,1.);
  ll->SetLineColor(kBlue);
  ll->SetLineWidth(5);

  TCanvas *c1 = new TCanvas("ROC","ROC",600.,600.);
  c1->cd();

  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1.0, 1.0);
  pad1->SetLeftMargin(0.15);
  pad1->SetBottomMargin(0.15);
  pad1->Draw();
  pad1->cd();

  g->Draw("AL");
  //p->Draw("Psame");

  c1->Print("RocAll" + variable + ".pdf","pdf");  
  gSystem->Exec("mv RocAll" + variable + ".pdf orthogonalCuts/" + sampleName);

  TCanvas *c2 = new TCanvas("Significance","Significance",600.,600.);
  c2->cd();

  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1.0, 1.0);
  pad2->SetLeftMargin(0.20);
  pad2->SetBottomMargin(0.15);
  pad2->Draw();
  pad2->cd();

  g2->Draw("AL");
  l->Draw("same");
  //  g4->Draw("same");

  c2->Print("SignificanceAll" + variable + ".pdf","pdf");  
  gSystem->Exec("mv SignificanceAll" + variable + ".pdf orthogonalCuts/" + sampleName);

  TCanvas *c3 = new TCanvas("Efficiencies","Efficiencies",600.,600.);
  c3->cd();

  TPad* pad3 = new TPad("pad3", "pad3", 0., 0., 1.0, 1.0);
  pad3->SetLeftMargin(0.15);
  pad3->SetRightMargin(0.15);
  pad3->SetBottomMargin(0.15);
  pad3->Draw();
  pad3->cd();

  g3->Draw("AL");
  g4->Draw("Lsame");
  ll->Draw("same");
  pad3->Update();
  TGaxis *axis = new TGaxis(pad3->GetUxmax(),pad3->GetUymin(),pad3->GetUxmax(),pad3->GetUymax(),0,1,510,"+L");
  axis->SetLabelSize(g3->GetXaxis()->GetLabelSize());
  axis->SetLabelFont(g3->GetXaxis()->GetLabelFont());
  axis->SetTitleSize(g3->GetXaxis()->GetTitleSize());
  axis->SetTitleFont(g3->GetXaxis()->GetTitleFont());
  axis->SetTitleOffset(1.4);
  axis->SetTitle("Background Rejection");
  axis->Draw("same");
  leg->Draw("same");

  c3->Print("EfficienciesAll" + variable + "." + printFormat,printFormat);
  gSystem->Exec("mv EfficienciesAll" + variable + "." + printFormat +" orthogonalCuts/" + sampleName);
  
  delete c1;
  delete c2;
  delete c3;
  delete g;
  delete g2;
  delete g3;
  delete g4;

  return 0;
}

void ROC(TString printMode = "png", TString sample = "Dark1"){

  TString var = ""; 
  TString cutType = "";
  Int_t rightBound = 0;
  TString tagText = "";

  ifstream inFile("letsROC.txt");
  std::string line;
  
  while (getline(inFile,line)){
    
    std::ofstream outFile("out.tmp",std::ios::out);
    outFile<<line<<endl;
    outFile.close();
    std::ifstream input;
    input.open("out.tmp",std::ios::in);
    input >> var >> cutType >> rightBound >> tagText;
    input.close();
    cout<<" "<<" "<< var <<" "<< sample <<" "<< cutType <<" "<< rightBound <<" "<< tagText<<" "<<endl;
    doTheROC(var, sample, cutType, rightBound, printMode, tagText);

    //drawPlots(var, leftBound, rightBound, nbin, units, fDark, fDark100, fDark10, fDark1000, fDark500, fZH, fHWW, fWW, printMode,logMode,normMode);
  } 
    
  inFile.close();
  gSystem->Exec("rm out.tmp");
}
