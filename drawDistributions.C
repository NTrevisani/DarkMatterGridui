#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"
#include <algorithm>

const UInt_t nProcesses = 6;//7 ;

enum {iData, iWW, iHWW, /*iQCD,*/ iTop, iTTJets, iVBF/*, iWJets*/};

TFile* input[nProcesses];


TString process[nProcesses];

process[iData]   = "Data2012";
process[iWW]     = "WW";
process[iHWW]    = "HWW";
//process[iQCD]    = "QCD";
process[iTop]    = "Top";
process[iTTJets] = "TTJets";
process[iVBF]    = "VBF";
//process[iWJets]  = "WJets";

Color_t color[nProcesses];

color[iData]   = kBlack;
color[iWW]     = kAzure-9;
color[iHWW]    = kRed;
//color[iQCD]    = kBlue;
color[iTop]    = kYellow+1;
color[iTTJets] = kYellow;
color[iVBF]    = kBlack;
//color[iWJets]  = kGreen;

// Relative errors directly from data cards (in %) 
//------------------------------------------------------------------------------

enum backgrounds {
   WW,
   HWW,
   QCD,
   Top,
   TTJets,
   VBF//,
   //   WJets
};
 

enum myChannels {
  jet0_SF=0*9,
  jet0_OF=1*9,
  jet1_SF=2*9,
  jet1_OF=3*9
};

const UInt_t Nsyst = 6 * 4;

//[*] WW / ggWW/  Vg/ Wjets/ VVV/ DY/ VV/ ggH/ VgS/ ?VBF?
//    WW / HWW / QCD / Top / VBF / WJets  

Double_t Systematics [Nsyst] = { 6.9,  30.4,  31.3,  27.4,  50.6,   2.6, //  5.9,  32.5,  34.1,  //0.1,   // jet0_SF
				 6.5,  30.4,  30.0,  28.3,  50.4,  30.7, //  6.2,  28.7,  30.3,  //0.1, // jet0_OF
				 13.7,  32.1,  33.3,  23.9,  50.3,   2.0,//   6.4,  41.0,  58.3, //0.1, // jet1_SF 
				 13.3,  31.8,  30.1,  26.4,  50.3,  30.4};//,   5.6,  35.5,  31.2};//, 0.1}; // jet1_OF
				 //9.5,	43.0,  43.4,  39.4,  71.4,  30.8,   8.6,  43.4,	45.6,   // jet0
				 //19.1,	45.2,  44.9,  35.6,  71.1,  30.5,   8.5,  54.2,	66.1 }; // jet1


Double_t systError[nProcesses];

systError[iData] = 0.0;





// Settings
//------------------------------------------------------------------------------
TString  _channel;
Int_t    _njet;
Double_t _luminosity;
TString  _format;
Bool_t   _drawRatio;
Bool_t   _dataDriven;
Bool_t   _savePlots;
Bool_t   _setLogy;


// Scale factors
//------------------------------------------------------------------------------
Double_t ttScale[] = {1.07, 1.09, 1.09, 1.09};
Double_t tWScale[] = {1.07, 1.09, 1.09, 1.09};
Double_t WWScale[] = {1.00, 1.00, 1.00, 1.00};
Double_t ZjScale[] = {1.00, 1.00, 1.00, 1.00};//{1.00, 4.20, 1.80, 4.00};


//------------------------------------------------------------------------------
// drawDistributions
//------------------------------------------------------------------------------
void drawDistributions(Int_t    njet       = 0,
		       TString  channel    = "OF",
		       Double_t luminosity = 19365,
		       TString  format     = "png",
		       Bool_t   drawRatio  = true,
		       Bool_t   dataDriven = false,
		       Bool_t   savePlots  = true,
		       Bool_t   setLogy    = false,
		       TString  scaleXS    = "1")
{

  //loading macros
  gInterpreter->LoadMacro ("draw/draw.C+");
  gInterpreter->ExecuteMacro("draw/ChargeRatioStyle.C");
  gStyle ->SetOptStat (0);
  gStyle ->SetPalette (1);

  _channel    = channel;
  _njet       = njet;
  _luminosity = luminosity;
  _format     = format;
  _drawRatio  = drawRatio;
  _dataDriven = dataDriven;
  _savePlots  = savePlots;
  _setLogy    = setLogy;

  gSystem->mkdir(_format, kTRUE);

  gStyle->SetHatchesLineWidth(1.00);
  gStyle->SetHatchesSpacing  (0.55);



  // Read input files
  //----------------------------------------------------------------------------
  //TString path = Form("rootfiles/%djet/%s/", _njet, _channel.Data());

  TString path = Form("rootFiles/AllJet/%s/", _channel.Data());
  path = path + scaleXS + "/";

  for (UInt_t ip=0; ip<nProcesses; ip++)
    input[ip] = new TFile(path + process[ip] + ".root", "read");



  // Differential distributions
  //----------------------------------------------------------------------------
  /*
  if (1) {

    DrawHistogram("hPtLepton1WWLevel_Diff",  "p_{T}^{max}"   ,  1, 0, "GeV","");
    DrawHistogram("hDileptonWWLevel_Diff",  "p_{T} (ll) "    ,  1, 0, "GeV","");
    DrawHistogram("hDeltaPhiWWLevel_Diff",  "#Delta#phi_{ll}",  1, 0, "rad","");
    DrawHistogram("hMinvWWLevel_Diff",  "m_{#font[12]{ll}}"  ,  1, 0, "GeV","");

  }
  */
  // DY distributions
  //----------------------------------------------------------------------------
  if (0) {
    DrawHistogram("hMassInZevents45.0",  "m_{#font[12]{ll}}", 1, 0, "GeV","", 76, 106);
    DrawHistogram("hMassOutZevents45.0", "m_{#font[12]{ll}}", 5, 0, "GeV","");
    
    DrawHistogram("hMassInZevents30.0",  "m_{#font[12]{ll}}", 1, 0, "GeV","", 76, 106);
    DrawHistogram("hMassOutZevents30.0", "m_{#font[12]{ll}}", 5, 0, "GeV","");
  }


  // Top distributions
  //----------------------------------------------------------------------------
  if (0) {
    DrawHistogram("hbTagDisNTopTaggedTopControlRegion", "2^{nd} jet TCHE","", 5, 1, "NULL","", -999, -999, false);
    DrawHistogram("hbTagDisNTopControlRegion",          "2^{nd} jet TCHE","", 5, 1, "NULL","", -999, -999, false);
    DrawHistogram("hbTagDisTopTaggedEvents",            "2^{nd} jet TCHE","", 5, 1, "NULL","", -999, -999, false);
  }


  // PAS distributions
  //----------------------------------------------------------------------------
  if (1) {
    for(Int_t _nje = 0; _nje < 4; ++_nje){

      TString nj;
      nj.Form("%d Jet",_nje);
      if (_nje == 2) nj.Form("%d+ Jet",_nje);
      if (_nje == 3) nj = "Inclusive";

      char name[80];

      sprintf(name,"hPtLepton1WWLevel%d",_nje);
      DrawHistogram(name,  "p_{T}^{max}",           10, 0, "GeV",nj.Data(), 20, 200);

      sprintf(name,"hpfMetWWLevel%d",_nje);
      DrawHistogram(name, "pfMET", 5, 0, "GeV",nj.Data());

      sprintf(name,"hPtDiLeptonWWLevel%d",_nje);
      DrawHistogram(name, "p_{T}^{#font[12]{ll}}", 5, 0, "GeV",nj.Data());

      sprintf(name,"hMinvWWLevel%d",_nje);
      DrawHistogram(name,       "m_{#font[12]{ll}}",     5, 0, "GeV",nj.Data());

      sprintf(name,"hDeltaPhiLeptonsWWLevel%d",_nje);
      DrawHistogram(name, "#Delta#phi_{ll}",  1, 0, "rad",nj.Data());

      
      sprintf(name,"hPtLepton1WWLevelNoHt%d",_nje);
      DrawHistogram(name,  "p_{T}^{max}",           10, 0, "GeV",nj.Data(), 20, 200);
      
      sprintf(name,"hpfMetWWLevelNoHt%d",_nje);
      DrawHistogram(name, "pfMET", 5, 0, "GeV",nj.Data());

      sprintf(name,"hPtDiLeptonWWLevelNoHt%d",_nje);
      DrawHistogram(name, "p_{T}^{#font[12]{ll}}", 5, 0, "GeV",nj.Data());

      sprintf(name,"hMinvWWLevelNoHt%d",_nje);
      DrawHistogram(name,       "m_{#font[12]{ll}}",     5, 0, "GeV",nj.Data());

      sprintf(name,"hDeltaPhiLeptonsWWLevelNoHt%d",_nje);
      DrawHistogram(name, "#Delta#phi_{ll}",  1, 0, "rad",nj.Data());

      sprintf(name,"hHt%d",_nje);
      DrawHistogram(name, "Ht"             , 10, 1, "GeV",nj.Data(),0.,500.);

      sprintf(name,"hHtAfter%d",_nje);
      DrawHistogram(name, "Ht"             , 10, 1, "GeV",nj.Data(),0.,500.);
      /*
      sprintf(name,"hMtLepton1WWLevel%d",_nje);
      DrawHistogram(name, "m_{T}"           , 10, 1, "GeV",nj.Data(),0.,200.);

      sprintf(name,"hMtLepton2WWLevel%d",_nje);
      DrawHistogram(name, "m_{T}"           , 10, 1, "GeV",nj.Data(),0.,200.);

      sprintf(name,"hMtSumWWLevel%d",_nje);
      DrawHistogram(name, "m_{T}"           , 10, 1, "GeV",nj.Data(),0.,350.);
      */
    }
    DrawHistogram("hWnJets", "number of jets",  1, 1,"","",0.,8.);
    DrawHistogram("hWnBtaggedJets", "number of jets",  1, 1,"","",0.,8.);


  }

  if(0){
    _dataDriven = false;
    systError[iWW]     = 0.1; // 10%
    systError[iHWW]    = 0.1; // 10%
    //systError[iQCD]    = 0.1; // 7% 
    systError[iTop]    = 0.1; // 8%;
    systError[iTTJets] = 0.1; // 8%;
    systError[iVBF]    = 0.1; // 32%
    //    systError[iWJets]  = 0.1; // 32%

    //    systError[ittFakes]   = 0.; // ?

    _setLogy = false;
    
    DrawHistogram("hPtLepton1TwoLeptonsLevel",  "p_{T}^{max}",           4, 0, "GeV","", 20, 200);
    DrawHistogram("hPtLepton2TwoLeptonsLevel",  "p_{T}^{min}",           4, 0, "GeV","", 20, 100);
    DrawHistogram("hPtDiLeptonTwoLeptonsLevel", "p_{T}^{#font[12]{ll}}", 4, 0, "GeV","", 20, 200);
    DrawHistogram("hMinvTwoLeptonsLevel",       "m_{#font[12]{ll}}",     4, 0, "GeV","");
    DrawHistogram("hNJetsPF30TwoLeptonsLevel",  "nJets",                 1, 0, "","");
    DrawHistogram("hDeltaRLeptonsTwoLeptonsLevel",   "#DeltaR_{ll}",     4, 0, "rad","");
    DrawHistogram("hDeltaPhiLeptonsTwoLeptonsLevel", "#Delta#phi_{ll}",  4, 0, "rad","");
    DrawHistogram("hpminMetTwoLeptonsLevel",         "min.proj.E_{T}",   10, 0, "GeV","", 0, 200);
    DrawHistogram("hppfMetTwoLeptonsLevel",          "proj.PFE_{T}",     10, 0, "GeV","", 0, 200);
    
  }
}


//------------------------------------------------------------------------------
// DrawHistogram
//------------------------------------------------------------------------------
void DrawHistogram(TString  hname,
		   TString  xtitle,
		   Int_t    ngroup       = -1,
		   Int_t    precision    = 1,
		   TString  units        = "NULL",
		   TString  jetChannel   = "0-jet",
		   Double_t xmin         = -999,
		   Double_t xmax         =  999,
		   Bool_t   moveOverflow = true)

{
  TCanvas* canvas;
 
  TPad* pad1;
  TPad* pad2;

  if (_drawRatio) {
    canvas = new TCanvas(hname, hname, 550, 1.2*600);
    canvas->Update();

    pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetTopMargin   (0.05);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3); 
    pad2->SetTopMargin   (0.08);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
  }
  else {
    canvas = new TCanvas(hname, hname, 550, 600);
    canvas -> Update();
    pad1 = (TPad*)canvas->GetPad(0);
  }


  //----------------------------------------------------------------------------
  // pad1
  //----------------------------------------------------------------------------
  pad1->cd();

  pad1->SetLogy(_setLogy);

  THStack* hstack = new THStack(hname, hname);

  TH1F* hist[nProcesses];

  for (UInt_t ip=0; ip<nProcesses; ip++) {

    hist[ip] = (TH1F*)input[ip]->Get(hname);
    hist[ip]->SetName(hname + process[ip]);

    if (ngroup > 0) hist[ip]->Rebin(ngroup);

    if (moveOverflow) MoveOverflowBins  (hist[ip], xmin, xmax);
    else              ZeroOutOfRangeBins(hist[ip], xmin, xmax);
    
    if (ip == iData) {
      hist[ip]->SetMarkerStyle(kFullCircle);
    }
    else {
      hist[ip]->SetFillColor(color[ip]);
      hist[ip]->SetFillStyle(1001);
      hist[ip]->SetLineColor(color[ip]);

      //if(ip != iData)
      //hist[ip]->Scale(1.25);

      //      if (ip == iH125) hist[ip]->SetFillStyle(3354);

      if (_dataDriven && ip == itt)    hist[ip]->Scale(ttScale[_njet]);
      if (_dataDriven && ip == itW)    hist[ip]->Scale(tWScale[_njet]);
      if (_dataDriven && ip == iWW)    hist[ip]->Scale(WWScale[_njet]);
      if (_dataDriven && ip == iqqWW)    hist[ip]->Scale(WWScale[_njet]);
      if (_dataDriven && ip == iDY)    hist[ip]->Scale(ZjScale[_njet]);
      if (_dataDriven && ip == iDYtau) hist[ip]->Scale(ZjScale[_njet]);

      hstack->Add(hist[ip]);
    }
  }


  // Include relative systematics 
  //----------------------------------------------------------------------------

  Int_t defineChannel = -1; 
  
  if ( _njet == 0 && _channel == "SF")  defineChannel  = jet0_SF;
  if ( _njet == 0 && _channel == "OF")  defineChannel  = jet0_OF;
  if ( _njet == 1 && _channel == "SF")  defineChannel  = jet1_SF;
  if ( _njet == 1 && _channel == "OF")  defineChannel  = jet1_OF;

  
  systError[iWW]     = 0.100; // 10%
  systError[iHWW]    = 0.100; // 10%
  //  systError[iQCD]    = 0.00; // 10%
  systError[iTop]    = 0.100; // 10%
  systError[iTTJets] = 0.100; // 10%
  systError[iVBF]    = 0.100; // 10%
  //  systError[iWJets]  = 0.100; // 10%

  
  //[*] WW / ggWW/  Vg/ Wjets/ VVV/ DY/ VV/ ggH/ VgS
  /*
  systError[iDY]     = Systematics [DYtt+defineChannel]/ 1e2;     
  systError[iWW]     = Systematics [ggWW+defineChannel]/ 1e2;
  systError[iqqWW]   = Systematics [WW+defineChannel]  / 1e2;
  systError[iWZ]     = Systematics [VV+defineChannel]  / 1e2;
  systError[iZZ]     = Systematics [VV+defineChannel]  / 1e2;
  systError[iWg]     = Systematics [Vg+defineChannel]  / 1e2;
  systError[iWgS]    = Systematics [VgS+defineChannel] / 1e2;
  systError[iWj]     = Systematics [Wjets+defineChannel] / 1e2;
  systError[iDYtau]  = Systematics [DYtt+defineChannel]/ 1e2;     
  systError[iZg]     = Systematics [Vg+defineChannel] / 1e2;
  systError[iVVV]    = Systematics [VVV+defineChannel]/ 1e2;
  systError[iH125]   = Systematics [ggH+defineChannel]/ 1e2;
  systError[iVBF]    = Systematics [ggH+defineChannel]/ 1e2;
  systError[ittFakes]    = Systematics [ggH+defineChannel]/ 1e2;
  */

  // All MC
  //----------------------------------------------------------------------------
  TH1F* allmc = (TH1F*)hist[iData]->Clone("allmc");

  allmc->SetFillColor  (kGray+2);
  allmc->SetFillStyle  (   3345);
  allmc->SetLineColor  (kGray+2);
  allmc->SetMarkerColor(kGray+2);
  allmc->SetMarkerSize (      0);

  for (UInt_t ibin=1; ibin<=allmc->GetNbinsX(); ibin++) {

    Double_t binValue = 0;
    Double_t binError = 0;

    for (UInt_t ip=0; ip<nProcesses; ip++) {

      if (ip == iData) continue;

      Double_t binContent = hist[ip]->GetBinContent(ibin);
      
      binValue += binContent;
      binError += (hist[ip]->GetBinError(ibin) * hist[ip]->GetBinError(ibin));

    
      if (_dataDriven) { 
	//binError += (systError[ip]*binContent * systError[ip]*binContent);
	//}
    
	if ( nProcesses != itW &&  nProcesses != itt) { 
	  binError += (systError[ip]*binContent * systError[ip]*binContent);
	  // luminosity syst
	  binError += (0.026*binContent * 0.026*binContent);
	} else {
	  binError += (systError[ip]*binContent*ttScale[_njet] * systError[ip]*binContent*ttScale[_njet]);
	  // luminosity syst
	  binError += (0.026*binContent *ttScale[_njet]* 0.026*binContent*tScale[_njet]);
	}
      }

      binError = sqrt(binError);

      allmc->SetBinContent(ibin, binValue);
      allmc->SetBinError  (ibin, binError);
    }

  }


  // Axis labels
  //----------------------------------------------------------------------------
  TAxis* xaxis = hist[iData]->GetXaxis();
  TAxis* yaxis = hist[iData]->GetYaxis();

  TString ytitle = Form("entries / %s.%df", "%", precision);

  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(Form(ytitle.Data(), hist[iData]->GetBinWidth(0)));
  yaxis->SetTitleOffset(1.6);

  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.025);
  xaxis->SetLabelSize  (   0.04);
  xaxis->SetNdivisions (   505);
  //xaxis->SetTitle      (xtitle);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.35);
  xaxis->SetTitleSize  (  0.045);
  
  yaxis->CenterTitle   (      );
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset(  0.02);
  yaxis->SetLabelSize  (   0.04);
  yaxis->SetNdivisions (   505);
  //yaxis->SetTitle      (ytitle);
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  1.80);
  yaxis->SetTitleSize  (  0.045);



  if (!units.Contains("NULL")) {
    
    xaxis->SetTitle(Form("%s [%s]", xaxis->GetTitle(), units.Data()));
    yaxis->SetTitle(Form("%s %s",   yaxis->GetTitle(), units.Data()));
  }


  // Draw
  //----------------------------------------------------------------------------
  xaxis->SetRangeUser(xmin, xmax);

  //  Float_t yrange = max(hist[iData]->GetBinContent(hist[iData]->GetMaximumBin()),);
  //  yaxis->SetRangeUser(0., 1.6*yrange);

  //  cout<<hist[iData]->GetBinContent(hist[iData]->GetMaximumBin())<<endl;

  hist[iData]->Draw("ep");

  hstack     ->Draw("hist,same");
  allmc      ->Draw("e2,same");
  hist[iData]->Draw("ep,same");


  // Adjust scale
  //----------------------------------------------------------------------------
  Float_t theMax   = GetMaximumIncludingErrors(hist[iData]);
  Float_t theMaxMC = GetMaximumIncludingErrors(allmc);

  //  if (theMaxMC > theMax) theMax = theMaxMC;

  if (pad1->GetLogy()) {

    theMax = TMath::Power(10, TMath::Log10(theMax) + 4.0);//2.7

    hist[iData]->SetMinimum(0.05);
  }
  if (theMax > theMaxMC) theMax *= 1.6;
  else theMax = theMaxMC * 1.6;

  hist[iData]->SetMaximum(theMax);
  
  // Legend
  //----------------------------------------------------------------------------
  Double_t yoffset = (_drawRatio) ? 0.054 : 0.048;
  Double_t x0      = (_drawRatio) ? 0.370 : 0.400; 

  TString allmcTitle = (_dataDriven) ? " stat #oplus syst" : " #sigma_{stat}";

  DrawLegend(0.23, 0.74 + 2.*(yoffset+0.001), hist[iData],   " data",    "lp", 0.035, 0.2, yoffset);
  DrawLegend(0.23, 0.74 + 1.*(yoffset+0.001), hist[iWW],     " WW",      "f",  0.035, 0.2, yoffset);
  DrawLegend(0.23, 0.74                     , hist[iHWW],    " HWW",     "f",  0.035, 0.2, yoffset);
  //DrawLegend(0.23, 0.74 - 1.*(yoffset+0.001), hist[iQCD],    " QCD",     "f",  0.035, 0.2, yoffset);
  DrawLegend(0.23, 0.74 - 1.*(yoffset+0.001), hist[iTop],    " Top",     "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74 + 1.*(yoffset+0.001), hist[iTTJets], " TTJets",  "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74                     , hist[iVBF],    " VBF",     "f",  0.035, 0.2, yoffset);  
  //  DrawLegend(  x0, 0.74 - 1.*(yoffset+0.001), hist[iWJets],  " WJets",   "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74 + 2.*(yoffset+0.001), allmc,       allmcTitle,   "f",  0.035, 0.2, yoffset);
  
  // Additional titles
  //----------------------------------------------------------------------------
  Double_t deltaY = (_drawRatio) ? 0.02 : 0.0;

  DrawTLatex(0.9, 0.860 + deltaY, 0.04, "CMS preliminary");
  DrawTLatex(0.9, 0.815 + deltaY, 0.03, Form("L = %.3f fb^{-1}", _luminosity/1e3));

  TString channelLabel; 

  //if (_njet == 0) channelLabel = "0-jet";
  channelLabel = jetChannel;

  //CJ
  DrawTLatex(0.9, 0.725 + deltaY, 0.06, _channel);
  DrawTLatex(0.9, 0.695 - deltaY, 0.06, channelLabel);
  
  if (_dataDriven) {
    DrawTLatex(0.9, 0.770 + deltaY, 0.06, "data-driven normalization");
  }

  //----------------------------------------------------------------------------
  // pad2
  //----------------------------------------------------------------------------
  if (_drawRatio) {

    pad2->cd();
    
    TH1F* ratio       = hist[iData]->Clone("ratio");
    TH1F* uncertainty = allmc->Clone("uncertainty");

    for (UInt_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {

      Double_t mcValue = allmc->GetBinContent(ibin);
      Double_t mcError = allmc->GetBinError  (ibin);

      Double_t dtValue = ratio->GetBinContent(ibin);
      Double_t dtError = ratio->GetBinError  (ibin);

      Double_t ratioValue       = (mcValue > 0) ? dtValue/mcValue : 0.0;
      Double_t ratioError       = (mcValue > 0) ? dtError/mcValue : 0.0;
      Double_t uncertaintyError = (mcValue > 0) ? mcError/mcValue : 0.0;

      ratio->SetBinContent(ibin, ratioValue);
      ratio->SetBinError  (ibin, ratioError);

      uncertainty->SetBinContent(ibin, 1.0);
      uncertainty->SetBinError  (ibin, uncertaintyError);
    }


    TAxis* uaxis = (TAxis*)uncertainty->GetXaxis();
    
    uaxis->SetRangeUser(xmin, xmax);
    
    
    uncertainty->Draw("e2");
    ratio      ->Draw("ep,same");

    uncertainty->GetYaxis()->SetRangeUser(0.4, 1.6);

    Pad2TAxis(uncertainty, hist[iData]->GetXaxis()->GetTitle(), "data / prediction");
  }


  // Save
  //----------------------------------------------------------------------------
  pad1->cd();
  pad1->GetFrame()->DrawClone();
  pad1->RedrawAxis();
  pad1->Update();

  if (_drawRatio) {
    pad2->cd();
    pad2->GetFrame()->DrawClone();
    pad2->RedrawAxis();
    pad2->Update();
  }

  if (_savePlots) {
    canvas->cd();
    pad1->Modified();
    if (_drawRatio) 
      pad2->Modified();
    canvas -> cd();
    canvas->Update();
    canvas->SaveAs(Form("%s/%sJet%s.%s",//("%s/%s_%djet_s.%s",
			_format.Data(),
			hname.Data(),
			//_njet,
			_channel.Data(),
			_format.Data())
		   );
  }
}


//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F*    h)

{
  Float_t maxWithErrors = 0;  
  
  Float_t binHeight = h->GetBinContent(h->GetMaximumBin());// + h->GetBinError(h->GetMaximumBin());
  if (binHeight > maxWithErrors) maxWithErrors = binHeight;

  return maxWithErrors;
}


//------------------------------------------------------------------------------
// MoveOverflowBins
//------------------------------------------------------------------------------
void MoveOverflowBins(TH1* h,
		      Double_t xmin,
		      Double_t xmax) const
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  Int_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax != -999) ? axis->FindBin(xmax) : nbins;

  Double_t firstVal = 0;
  Double_t firstErr = 0;

  Double_t lastVal = 0;
  Double_t lastErr = 0;

  for (UInt_t i=0; i<=nbins+1; i++) {

    if (i <= firstBin) {
      firstVal += h->GetBinContent(i);
      firstErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i >= lastBin) {
      lastVal += h->GetBinContent(i);
      lastErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i < firstBin || i > lastBin) {
      h->SetBinContent(i, 0);
      h->SetBinError  (i, 0);
    }
  }

  firstErr = sqrt(firstErr);
  lastErr  = sqrt(lastErr);

  h->SetBinContent(firstBin, firstVal);
  h->SetBinError  (firstBin, firstErr);

  h->SetBinContent(lastBin, lastVal);
  h->SetBinError  (lastBin, lastErr);
}


//------------------------------------------------------------------------------
// ZeroOutOfRangeBins
//------------------------------------------------------------------------------
void ZeroOutOfRangeBins(TH1* h, Double_t xmin, Double_t xmax) const
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  Int_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax != -999) ? axis->FindBin(xmax) : nbins;

  for (UInt_t i=0; i<=nbins+1; i++) {

    if (i < firstBin || i > lastBin) {
      h->SetBinContent(i, 0);
      h->SetBinError  (i, 0);
    }
  }
}


//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t x, Double_t y, Double_t tsize, const char* text)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC();
  tl->SetTextAlign(   32);
  tl->SetTextFont (   42);
  tl->SetTextSize (tsize);

  tl->Draw("same");
}


//------------------------------------------------------------------------------
// Pad2TAxis
//------------------------------------------------------------------------------
void Pad2TAxis(TH1* hist, TString xtitle, TString ytitle)
{
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.025);
  xaxis->SetLabelSize  (   0.1);
  xaxis->SetNdivisions (   505);
  xaxis->SetTitle      (xtitle);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.35);
  xaxis->SetTitleSize  (  0.11);
  
  yaxis->CenterTitle   (      );
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset(  0.02);
  yaxis->SetLabelSize  (   0.1);
  yaxis->SetNdivisions (   505);
  yaxis->SetTitle      (ytitle);
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  0.75);
  yaxis->SetTitleSize  (  0.11);
}

//------------------------------------------------------------------------------                                                               
// Draw Legend
//------------------------------------------------------------------------------                                                                
TLegend* DrawLegend(Float_t x1,
                    Float_t y1,
                    TH1*    hist,
                    TString label,
                    TString option,
                    Float_t tsize,
                    Float_t xoffset,
                    Float_t yoffset)
{
  TLegend* legend = new TLegend(x1,
                                y1,
                                x1 + xoffset,
                                y1 + yoffset);

  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (tsize);

  legend->AddEntry(hist, label.Data(), option.Data());
  legend->Draw();

  return legend;
}
