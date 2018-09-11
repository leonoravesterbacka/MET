#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <list>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include <math.h>
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
using namespace std;

TGraphErrors *h1= new TGraphErrors;
TGraphErrors *h2= new TGraphErrors;
TGraphErrors *h1up= new TGraphErrors;
TGraphErrors *h1down = new TGraphErrors;
TGraphErrors *h1EE = new TGraphErrors;
TGraphErrors *h2EE = new TGraphErrors;
TGraphErrors *hMCEEup = new TGraphErrors;
TGraphErrors *hMCEEdown = new TGraphErrors;
TGraphErrors *hMCEEjerup = new TGraphErrors;
TGraphErrors *hMCEEjerdown = new TGraphErrors;
TGraphErrors *hMCEEupUncl = new TGraphErrors;
TGraphErrors *hMCEEdownUncl = new TGraphErrors;
TGraphErrors *h1MM = new TGraphErrors;
TGraphErrors *h2MM = new TGraphErrors;
TGraphErrors *hMCMMup = new TGraphErrors;
TGraphErrors *hMCMMdown = new TGraphErrors;
TGraphErrors *hMCMMjerup = new TGraphErrors;
TGraphErrors *hMCMMjerdown = new TGraphErrors;
TGraphErrors *hMCMMupUncl = new TGraphErrors;
TGraphErrors *hMCMMdownUncl = new TGraphErrors;
TGraphErrors *h1GG = new TGraphErrors;
TGraphErrors *h2GG = new TGraphErrors;
TGraphErrors *hMCGGup = new TGraphErrors;
TGraphErrors *hMCGGdown = new TGraphErrors;
TGraphErrors *hMCGGjerup = new TGraphErrors;
TGraphErrors *hMCGGjerdown = new TGraphErrors;
TGraphErrors *hMCGGupUncl = new TGraphErrors;
TGraphErrors *hMCGGdownUncl = new TGraphErrors;


void testRatioEEMuMuGJetsComparison(TString, TString, TString);                                            
//TMultiGraph mg;
Double_t chi2;
Double_t mm_da_0;
Double_t mm_da_0_e;
Double_t ee_da_0;
Double_t ee_da_0_e;
Double_t gg_da_0;
Double_t gg_da_0_e;
Double_t mm_da_s;
Double_t mm_da_s_e;
Double_t ee_da_s;
Double_t ee_da_s_e;
Double_t gg_da_s;
Double_t gg_da_s_e;
Double_t mm_mc_0;
Double_t mm_mc_0_e;
Double_t ee_mc_0;
Double_t ee_mc_0_e;
Double_t gg_mc_0;
Double_t gg_mc_0_e;
Double_t mm_mc_s;
Double_t mm_mc_s_e;
Double_t ee_mc_s;
Double_t ee_mc_s_e;
Double_t gg_mc_s;
Double_t gg_mc_s_e;
Double_t mm_da_p0;
Double_t mm_da_p0_e;
Double_t ee_da_p0;
Double_t ee_da_p0_e;
Double_t gg_da_p0;
Double_t gg_da_p0_e;
Double_t mm_da_p;
Double_t mm_da_p_e;
Double_t ee_da_p;
Double_t ee_da_p_e;
Double_t gg_da_p;
Double_t gg_da_p_e;
Double_t mm_mc_p0;
Double_t mm_mc_p0_e;
Double_t ee_mc_p0;
Double_t ee_mc_p0_e;
Double_t gg_mc_p0;
Double_t gg_mc_p0_e;
Double_t mm_mc_p;
Double_t mm_mc_p_e;
Double_t ee_mc_p;
Double_t ee_mc_p_e;
Double_t gg_mc_p;
Double_t gg_mc_p_e;
Double_t mm_damc_0_e;
Double_t mm_damc_s_e;
Double_t ee_damc_0_e;
Double_t ee_damc_s_e;
Double_t gg_damc_0_e;
Double_t gg_damc_s_e;
Double_t mm_damc_p0_e;
Double_t mm_damc_p_e;
Double_t ee_damc_p0_e;
Double_t ee_damc_p_e;
Double_t gg_damc_p0_e;
Double_t gg_damc_p_e;

TGraphAsymmErrors staterror(TH1F *, TH1F * );
TGraphAsymmErrors JESerror(TH1F *, TH1F *, TH1F *, TH1F * );
TGraphAsymmErrors JERerror(TH1F *, TH1F *, TH1F *, TH1F * , TH1F *, TH1F *);
TGraphAsymmErrors TOTTOTerror(TH1F *, TH1F *, TH1F *, TH1F *, TH1F *, TH1F * , TH1F *, TH1F * );
TGraphAsymmErrors TOTNEWerror(TH1F *, TH1F *, TH1F *, TH1F *, TH1F *, TH1F * , TH1F *);

void PlotWithRatioEEMuMuGJetsComparison(TCanvas *c1, TGraphErrors *mcdyee,  TGraphErrors *datadyee,  TGraphErrors *mcdymm,  TGraphErrors *datadymm,TGraphErrors *mcdygg,  TGraphErrors *datadygg,  TGraphErrors *mcdyeeup, TGraphErrors *mcdyeedown,  TGraphErrors *mcdyeejerup, TGraphErrors *mcdyeejerdown, TGraphErrors *mcdyeeupUncl, TGraphErrors *mcdyeedownUncl,  TGraphErrors *mcdymmup, TGraphErrors *mcdymmdown,  TGraphErrors *mcdymmjerup, TGraphErrors *mcdymmjerdown, TGraphErrors *mcdymmupUncl, TGraphErrors *mcdymmdownUncl, TGraphErrors *mcdyggup, TGraphErrors *mcdyggdown,  TGraphErrors *mcdyggjerup, TGraphErrors *mcdyggjerdown, TGraphErrors *mcdyggupUncl, TGraphErrors *mcdyggdownUncl, TString EjeX, TString histogramaname, TString Type) {

TH1::SetDefaultSumw2();
  c1->cd(); 
//Float_t xbins[] = { 4, 6,8, 10, 12,14, 16, 18, 20,22, 24, 26, 28, 32,  40};
//     double nbins=14.;
                 //this is the OG one Float_t xbins[] = {400,500,600, 700,800, 900, 1000, 1100, 1300, 1600, 2500};
//Float_t xbins[] = {300, 450,550, 650, 750, 820, 900, 990, 1090, 1210, 1400, 1700, 2400}; 
Float_t xbins[] = {350,500, 600, 700, 800, 900, 1000, 1100, 1200, 1400, 1700, 2400}; 
   //this is the rightFloat_t xbins[] = {200, 360, 410, 460, 510, 570, 630, 690,  750, 820, 900, 990, 1090, 1210, 1400, 1700, 2400}; 
double nbins=11.;
 
//double nbins=12.;
// Float_t xbins[] = {50, 60, 70, 80, 100,120, 140,170, 200,250, 330};
// double nbins=9.;
 //     cout << "plotting resolution as a function of qt" << endl;
 //   if (histogramaname.Contains("over")){
//  Float_t xbins[] = {18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 210, 235, 260, 290, 350};
//    double nbins=19.;
//  Float_t xbins[] = {0, 10, 18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 200, 225, 250, 275, 305, 335, 365, 400, 440};
//  Float_t xbins[] = {0, 10, 18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 200, 225, 250, 275, 305, 335, 365, 400, 440, 500};
//    double nbins=25.;
//    double nbins=26.;
 //     cout << "plotting scale " << endl;
 //   }



   TH1F * hdatadyee= new TH1F("","", nbins,xbins);
   TH1F * hmcdyee=        new TH1F("hmcdyee"+histogramaname+Type,        "hmcdyee"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyeeup=      new TH1F("hmcdyeeup"+histogramaname+Type,      "hmcdyeeup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyeedown=    new TH1F("hmcdyeedown"+histogramaname+Type,    "hmcdyeedown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyeejerup=   new TH1F("hmcdyeejerup"+histogramaname+Type,   "hmcdyeejerup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyeejerdown= new TH1F("hmcdyeejerdown"+histogramaname+Type, "hmcdyeejerdown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyeeupUncl=  new TH1F("hmcdyeeupUncl"+histogramaname+Type,  "hmcdyeeupUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyeedownUncl=new TH1F("hmcdyeedownUncl"+histogramaname+Type,"hmcdyeedownUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hdadyee=        new TH1F("hdadyee"+histogramaname+Type,        "hdadyee"+histogramaname+Type, nbins,xbins);          
   TH1F * hmcdymm=        new TH1F("hmcdymm"+histogramaname+Type,        "hmcdymm"+histogramaname+Type, nbins,xbins);           
   TH1F * hmcdymmup=      new TH1F("hmcdymmup"+histogramaname+Type,      "hmcdymmup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdymmdown=    new TH1F("hmcdymmdown"+histogramaname+Type,    "hmcdymmdown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdymmjerup=   new TH1F("hmcdymmjerup"+histogramaname+Type,   "hmcdymmjerup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdymmjerdown= new TH1F("hmcdymmjerdown"+histogramaname+Type, "hmcdymmjerdown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdymmupUncl=  new TH1F("hmcdymmupUncl"+histogramaname+Type,  "hmcdymmupUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdymmdownUncl=new TH1F("hmcdymmdownUncl"+histogramaname+Type,"hmcdymmdownUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hdadymm=        new TH1F("hdadymm"+histogramaname+Type,        "hdadymm"+histogramaname+Type, nbins,xbins);          
   TH1F * hmcdygg=        new TH1F("hmcdygg"+histogramaname+Type,        "hmcdygg"+histogramaname+Type, nbins,xbins);          
   TH1F * hmcdyggup=      new TH1F("hmcdyggup"+histogramaname+Type,      "hmcdyggup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyggdown=    new TH1F("hmcdyggdown"+histogramaname+Type,    "hmcdyggdown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyggjerup=   new TH1F("hmcdyggjerup"+histogramaname+Type,   "hmcdyggjerup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyggjerdown= new TH1F("hmcdyggjerdown"+histogramaname+Type, "hmcdyggjerdown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyggupUncl=  new TH1F("hmcdyggupUncl"+histogramaname+Type,  "hmcdyggupUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyggdownUncl=new TH1F("hmcdyggdownUncl"+histogramaname+Type,"hmcdyggdownUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hdadygg=        new TH1F("hdadygg"+histogramaname+Type,        "hdadygg"+histogramaname+Type, nbins,xbins);          
   TH1F * hdatadymm=      new TH1F("hdatadymm"+histogramaname+Type,"hdatadymm"+histogramaname+Type, nbins,xbins);
   TH1F * hdatadygg=      new TH1F("hdatadygg"+histogramaname+Type,"hdatadygg"+histogramaname+Type,  nbins,xbins);
   //TH1F * hmcdygg= new TH1F("hmcdygg"+histogramaname+Type,"hmcdygg"+histogramaname+Type, nbins,xbins);
   
   for (int i = 0, n = nbins; i < n; ++i) {
       hdatadyee->SetBinContent(i+1,datadyee->GetY()[i]);
       hmcdyee->SetBinContent(i+1,mcdyee->GetY()[i]);
 	   hmcdyeedown->    SetBinContent(i+1,mcdyeedown->GetY()[i]);
       hmcdyeeup->      SetBinContent(i+1,mcdyeeup->GetY()[i]);     
       hmcdyeejerdown-> SetBinContent(i+1,mcdyeejerdown->GetY()[i]);
       hmcdyeejerup->   SetBinContent(i+1,mcdyeejerup->GetY()[i]);     
       hmcdyeedownUncl->SetBinContent(i+1,mcdyeedownUncl->GetY()[i]);
       hmcdyeeupUncl->  SetBinContent(i+1,mcdyeeupUncl->GetY()[i]);      
       hmcdymmdown->    SetBinContent(i+1,mcdymmdown->GetY()[i]);
       hmcdymmup->      SetBinContent(i+1,mcdymmup->GetY()[i]);     
       hmcdymmjerdown-> SetBinContent(i+1,mcdymmjerdown->GetY()[i]);
       hmcdymmjerup->   SetBinContent(i+1,mcdymmjerup->GetY()[i]);     
       hmcdymmdownUncl->SetBinContent(i+1,mcdymmdownUncl->GetY()[i]);
       hmcdymmupUncl->  SetBinContent(i+1,mcdymmupUncl->GetY()[i]);      
       hmcdyggdown->    SetBinContent(i+1,mcdyggdown->GetY()[i]);
       hmcdyggup->      SetBinContent(i+1,mcdyggup->GetY()[i]);     
       hmcdyggjerdown-> SetBinContent(i+1,mcdyggjerdown->GetY()[i]);
       hmcdyggjerup->   SetBinContent(i+1,mcdyggjerup->GetY()[i]);     
       hmcdyggdownUncl->SetBinContent(i+1,mcdyggdownUncl->GetY()[i]);
       hmcdyggupUncl->  SetBinContent(i+1,mcdyggupUncl->GetY()[i]);      
       hdatadyee->SetBinError(i+1,datadyee->GetErrorY(i));
       hmcdyee->SetBinError(i+1,mcdyee->GetErrorY(i));          
       hdatadymm->SetBinContent(i+1,datadymm->GetY()[i]);
       hmcdymm->SetBinContent(i+1,mcdymm->GetY()[i]);
       hdatadymm->SetBinError(i+1,datadymm->GetErrorY(i));
       hmcdymm->SetBinError(i+1,mcdymm->GetErrorY(i));                  
       hdatadygg->SetBinContent(i+1,datadygg->GetY()[i]);
       hmcdygg->SetBinContent(i+1,mcdygg->GetY()[i]);
       hdatadygg->SetBinError(i+1,datadygg->GetErrorY(i));
       hmcdygg->SetBinError(i+1,mcdygg->GetErrorY(i));                  
    }

    TGraphAsymmErrors toterrorEEMC= TOTNEWerror(hmcdyee,hmcdyeedown,hmcdyeeup, hmcdyeejerdown,hmcdyeejerup,hmcdyeedownUncl,hmcdyeeupUncl);
    TGraphAsymmErrors toterrorMMMC= TOTNEWerror(hmcdymm,hmcdymmdown,hmcdymmup, hmcdymmjerdown,hmcdymmjerup,hmcdymmdownUncl,hmcdymmupUncl);
    TGraphAsymmErrors toterrorGGMC= TOTNEWerror(hmcdygg,hmcdyggdown,hmcdyggup, hmcdyggjerdown,hmcdyggjerup,hmcdyggdownUncl,hmcdyggupUncl);
    //TGraphAsymmErrors toterrorMMMC= TOTNEWerror(hmcdymm,hmcdymmdown,hmcdymmup, hmcdymmjerdown,hmcdymmjerup,hmcdymm,hmcdymm);
    //TGraphAsymmErrors toterrorEEMC= TOTNEWerror(hmcdyee,hmcdyeedown,hmcdyeeup, hmcdyeejerdown,hmcdyeejerup,hmcdyee,hmcdyee);
    //TGraphAsymmErrors toterrorGGMC= TOTNEWerror(hmcdygg,hmcdyggdown,hmcdyggup, hmcdyggjerdown,hmcdyggjerup,hmcdygg,hmcdygg);

   for (int i = 0, n = nbins; i < n; ++i) {
       hmcdyee->SetBinError(i+1,toterrorEEMC.GetErrorYhigh(i));
       hmcdymm->SetBinError(i+1,toterrorMMMC.GetErrorYhigh(i));
       hmcdygg->SetBinError(i+1,toterrorGGMC.GetErrorYhigh(i));
    }

    hdatadyee->SetMarkerStyle(22);
    hdatadyee->SetMarkerColor(46);
    hdatadyee->SetLineColor(46);
    hmcdyee->SetMarkerColor(46);
    hmcdyee->SetMarkerStyle(22); 
    hdatadymm->SetMarkerStyle(23);
    hdatadymm->SetMarkerColor(4);
    hdatadymm->SetLineColor(4);
    hmcdymm->SetMarkerColor(22);
    hmcdymm->SetMarkerStyle(23);
    hdatadygg->SetMarkerStyle(20);
    hdatadygg->SetMarkerColor(8);
    hdatadygg->SetLineColor(8);
    hmcdygg->SetMarkerColor(8);
    hmcdygg->SetMarkerStyle(20);    
    TH1F *ratiodyee=(TH1F*)hdatadyee->Clone(); 
    TH1F *ratiodymm=(TH1F*)hdatadymm->Clone(); 
    TH1F *ratiodygg=(TH1F*)hdatadygg->Clone(); 
   
    ratiodyee->Divide(hmcdyee);
    ratiodymm->Divide(hmcdymm);
    ratiodygg->Divide(hmcdygg);
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.014);
    pad1->SetLeftMargin(0.1);
    pad1->Draw();
    pad1->cd();
    hdatadyee->GetYaxis()->SetRangeUser(0., 40.);
    hdatadyee->GetXaxis()->SetLabelFont(43); //font in pixels
    hdatadyee->GetXaxis()->SetLabelSize(1); //in pixels
    hdatadyee->GetYaxis()->SetLabelFont(43);
    hdatadyee->GetYaxis()->SetLabelSize(16);
    hdatadyee->GetYaxis()->SetTitleFont(43); //font in pixels
    hdatadyee->GetYaxis()->SetTitleSize(20); //in pixels
    hdatadyee->GetXaxis()->SetTitleOffset(3.3);
    hdatadyee->SetFillColor(54);
    hdatadyee->SetFillStyle(3345);                                                                                                                                    
    hdatadyee->Draw("e1");
    pad1->Update();
    TString ylabel="";

if (histogramaname=="met_uPara_over_qt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_uParaRaw_over_qt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_uPara_over_zll_pt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_uParaPuppi_over_zll_pt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_uParaRaw_over_zll_pt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_uParaPuppi_over_zll_pt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_uPara__vs_zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerp__vs_zll_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_met_sumEt-zll_pt") ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_met_sumEt-zll_pt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPara__vs"))ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname.Contains("uPara_vs"))ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname.Contains("uParaPuppi__vs"))ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname.Contains("uPerp__vs"))ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPerp_vs"))ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPerpPuppi__vs"))ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.5,1.1);  


    hdatadymm->Draw("e1, X0 same");
    hdatadyee->Draw("AXIS");
    hdatadygg->Draw("e1, X0 same");
    hdatadymm->Draw("e1, X0 same");
    if (Type == "PF_Type1FIT"   or Type == "Puppi_Type1FIT" or Type == "PF_Type1FITSumEt" or Type == "PF_Type1FITnVert"){
        hdatadyee->Draw("e1, X0 same");
    }
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hdatadyee->GetXaxis ()->GetXmin (), 1, hdatadyee->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw();
    }                                                                                                              
    
    TLegend* leg(0);
    leg = new TLegend(0.25,0.05,0.77,0.26);
    if (histogramaname.Contains("over")){
    leg = new TLegend(0.25,0.05,0.77,0.26);
    }
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillColor(0);
    if (Type == "PF_Type1" or Type == "PF_Type1FIT"  or Type == "Puppi_Type1FIT" or Type == "PF_Type1MC" or Type == "PF_Type1FITnVert" or Type == "PF_Type1FITSumEt"){
    leg->AddEntry(hdatadymm, "Type 1 PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "Type 1 PF p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadygg, "Type 1 PF p_{T}^{miss} #gamma, q_{T} > 50 GeV", "lp");
    }                                                                                     
    if (Type == "Puppi_Type1"){
    leg->AddEntry(hdatadymm, "Type 1 Puppi p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadygg, "Type 1 Puppi p_{T}^{miss} Z#rightarrow ee", "lp");
    }                                                                                     
    if (histogramaname.Contains("sum")){
    cout.precision(3);
    //TF1 *f2 = new TF1("f2","[0]+[1]*sqrt(x)");
//    if(histogramaname == "met_uPara__vs_met_sumEt-zll_pt"){
//        cout << "uPara @@@@@@@@@@@@" << endl;
//           TF1 *f2 = new TF1("f2","[0]+[1]*sqrt(x)",400, 2500);
//           //TF1 *f2 = new TF1("f2","[0]+[1]*sqrt(x)",200, 2500);
//    //f2->SetParLimits(0,0,10);
//    //Double_t chi2 = f2->GetChisquare();
//    //f2->SetParLimits(0,0,2);
//    }
//    else{TF1 *f2 = new TF1("f2","[0]+[1]*sqrt(x)", 400,2500);}
//gStyle->SetOptFit(1111);
    TF1 *f1 = new TF1("f1","[0]+[1]*sqrt(x)", 350,2400);       
    TF1 *f2 = new TF1("f2","[0]+[1]*sqrt(x)", 600,2400);       
    TF1 *f3 = new TF1("f3","[0]+[1]*sqrt(x)", 350,2400);       
    TF1 *f4 = new TF1("f4","[0]+[1]*sqrt(x)", 350,2400);       
    TF1 *f5 = new TF1("f5","[0]+[1]*sqrt(x)", 350,2400);       
    TF1 *f6 = new TF1("f6","[0]+[1]*sqrt(x)", 350,2400);
    f1->SetLineColor(2);f1->SetLineStyle(2); 
    f2->SetLineColor(2); 
    f3->SetLineColor(4);f3->SetLineStyle(2); 
    f4->SetLineColor(4); 
    f5->SetLineColor(3);f5->SetLineStyle(2); 
    f6->SetLineColor(3);                                      
//    f1->SetParLimits(0,0,100);
//    f2->SetParLimits(0,0,100);
//    f3->SetParLimits(0,0,100);
//    f4->SetParLimits(0,0,100);
//    f5->SetParLimits(0,0,100);
//    f6->SetParLimits(0,0,100);

    hdatadyee->Fit("f1", "R, SAME");                                     
    ee_da_s   = TVirtualFitter::GetFitter()->GetParameter(1); 
    ee_da_s_e = TVirtualFitter::GetFitter()->GetParError(1);  
    ee_da_0   = TVirtualFitter::GetFitter()->GetParameter(0); 
    ee_da_0_e = TVirtualFitter::GetFitter()->GetParError(0);  
    hmcdyee->Fit("f2", "R", "SAME");                                       
    ee_mc_s   = TVirtualFitter::GetFitter()->GetParameter(1); 
    ee_mc_s_e = TVirtualFitter::GetFitter()->GetParError(1);  
    ee_mc_0   = TVirtualFitter::GetFitter()->GetParameter(0); 
    ee_mc_0_e = TVirtualFitter::GetFitter()->GetParError(0);  
    hdatadymm->Fit("f3", "R", "SAME");
    mm_da_s   = TVirtualFitter::GetFitter()->GetParameter(1);  
    mm_da_s_e = TVirtualFitter::GetFitter()->GetParError(1);  
    mm_da_0   = TVirtualFitter::GetFitter()->GetParameter(0);  
    mm_da_0_e = TVirtualFitter::GetFitter()->GetParError(0);  
    hmcdymm->Fit("f4", "R", "SAME");                                       
    mm_mc_s   = TVirtualFitter::GetFitter()->GetParameter(1); 
    mm_mc_s_e = TVirtualFitter::GetFitter()->GetParError(1);  
    mm_mc_0   = TVirtualFitter::GetFitter()->GetParameter(0); 
    mm_mc_0_e = TVirtualFitter::GetFitter()->GetParError(0);  
    hdatadygg->Fit("f5", "R", "SAME");                                      
    gg_da_s   = TVirtualFitter::GetFitter()->GetParameter(1); 
    gg_da_s_e = TVirtualFitter::GetFitter()->GetParError(1);  
    gg_da_0   = TVirtualFitter::GetFitter()->GetParameter(0); 
    gg_da_0_e = TVirtualFitter::GetFitter()->GetParError(0);  
    hmcdygg->Fit("f6", "R", "SAME");                                       
    gg_mc_s   = TVirtualFitter::GetFitter()->GetParameter(1); 
    gg_mc_s_e = TVirtualFitter::GetFitter()->GetParError(1);  
    gg_mc_0   = TVirtualFitter::GetFitter()->GetParameter(0); 
    gg_mc_0_e = TVirtualFitter::GetFitter()->GetParError(0); 
    mm_damc_0_e = (sqrt((mm_da_0_e/mm_da_0)*(mm_da_0_e/mm_da_0) +(mm_mc_0_e/mm_mc_0)*(mm_mc_0_e/mm_mc_0)))*(mm_da_0/mm_mc_0);  
    mm_damc_s_e = (sqrt((mm_da_s_e/mm_da_s)*(mm_da_s_e/mm_da_s) +(mm_mc_s_e/mm_mc_s)*(mm_mc_s_e/mm_mc_s)))*(mm_da_s/mm_mc_s);  
    ee_damc_0_e = (sqrt((ee_da_0_e/ee_da_0)*(ee_da_0_e/ee_da_0) +(ee_mc_0_e/ee_mc_0)*(ee_mc_0_e/ee_mc_0)))*(ee_da_0/ee_mc_0);  
    ee_damc_s_e = (sqrt((ee_da_s_e/ee_da_s)*(ee_da_s_e/ee_da_s) +(ee_mc_s_e/ee_mc_s)*(ee_mc_s_e/ee_mc_s)))*(ee_da_s/ee_mc_s);  
    gg_damc_0_e = (sqrt((gg_da_0_e/gg_da_0)*(gg_da_0_e/gg_da_0) +(gg_mc_0_e/gg_mc_0)*(gg_mc_0_e/gg_mc_0)))*(gg_da_0/gg_mc_0);  
    gg_damc_s_e = (sqrt((gg_da_s_e/gg_da_s)*(gg_da_s_e/gg_da_s) +(gg_mc_s_e/gg_mc_s)*(gg_mc_s_e/gg_mc_s)))*(gg_da_s/gg_mc_s); 

    if(histogramaname=="met_uPara__vs_met_sumEt-zll_pt"){cout << " uPara " << endl;}
    else{cout << " uPerp " << endl;}                                              
    
    //cout<<"\n Z$\\rightarrow \\mu\\mu$   & "<< mm_da_0<<" $\\pm$ "<<mm_da_0_e<<" & " <<mm_da_0/mm_mc_0<<" $\\pm$ "<< mm_damc_0_e <<" & "<< mm_da_s<<" $\\pm$ "<<mm_da_s_e<<" & " <<mm_da_s/mm_mc_s<<" $\\pm$ "<< mm_damc_s_e<<"\\\\"<<endl; 
    //cout<<"\n Z$\\rightarrow$ ee       & "<< ee_da_0<<" $\\pm$ "<<ee_da_0_e<<" & " <<ee_da_0/ee_mc_0<<" $\\pm$ "<< ee_damc_0_e <<" & "<< ee_da_s<<" $\\pm$ "<<ee_da_s_e<<" & " <<ee_da_s/ee_mc_s<<" $\\pm$ "<< ee_damc_s_e<<"\\\\"<<endl; 
    //cout<<"\n $\\gamma$+jets           & "<< gg_da_0<<" $\\pm$ "<<gg_da_0_e<<" & " <<gg_da_0/gg_mc_0<<" $\\pm$ "<< gg_damc_0_e <<" & "<< gg_da_s<<" $\\pm$ "<<gg_da_s_e<<" & " <<gg_da_s/gg_mc_s<<" $\\pm$ "<< gg_damc_s_e<<"\\\\"<<endl; 
    //cout<<"\n"<< endl; 
    cout<<"\n Z$\\rightarrow \\mu\\mu$   & "<< mm_da_0<<" $\\pm$ "<<mm_da_0_e<<" & " <<mm_mc_0<<" $\\pm$ "<< mm_mc_0_e <<" & "<< mm_da_s<<" $\\pm$ "<<mm_da_s_e<<" & " <<mm_da_s/mm_mc_s<<" $\\pm$ "<< mm_damc_s_e<<"\\\\"<<endl; 
    cout<<"\n Z$\\rightarrow$ ee         & "<< ee_da_0<<" $\\pm$ "<<ee_da_0_e<<" & " <<ee_mc_0<<" $\\pm$ "<< ee_mc_0_e <<" & "<< ee_da_s<<" $\\pm$ "<<ee_da_s_e<<" & " <<ee_da_s/ee_mc_s<<" $\\pm$ "<< ee_damc_s_e<<"\\\\"<<endl; 
    cout<<"\n $\\gamma$+jets             & "<< gg_da_0<<" $\\pm$ "<<gg_da_0_e<<" & " <<gg_mc_0<<" $\\pm$ "<< gg_mc_0_e <<" & "<< gg_da_s<<" $\\pm$ "<<gg_da_s_e<<" & " <<gg_da_s/gg_mc_s<<" $\\pm$ "<< gg_damc_s_e<<"\\\\"<<endl; 
    cout<<"\n"<< endl; 
   }                                                                                                                                                             
    
          if (histogramaname.Contains("nVert")){
          cout.precision(3);
          if(histogramaname == "met_uParaPuppi__vs_nVert"){
              cout << " uPara " << endl;
              TF1 *f1 = new TF1("f1","sqrt([0]**2+ (x*([1]**2))/0.7)",8, 40);
          }
          else{
          TF1 *f1 = new TF1("f1","sqrt([0]**2+ (x*([1]**2))/0.7)",0, 40);
          }
          hdatadymm->Fit("f1", "R");
          mm_da_p   = TVirtualFitter::GetFitter()->GetParameter(1);  
          mm_da_p_e = TVirtualFitter::GetFitter()->GetParError(1);  
          mm_da_p0   = TVirtualFitter::GetFitter()->GetParameter(0);  
          mm_da_p0_e = TVirtualFitter::GetFitter()->GetParError(0);  
          hmcdymm->Fit("f1", "R");                                       
          mm_mc_p   = TVirtualFitter::GetFitter()->GetParameter(1); 
          mm_mc_p_e = TVirtualFitter::GetFitter()->GetParError(1);  
          mm_mc_p0   = TVirtualFitter::GetFitter()->GetParameter(0); 
          mm_mc_p0_e = TVirtualFitter::GetFitter()->GetParError(0);  
          hdatadyee->Fit("f1", "R");                                      
          ee_da_p   = TVirtualFitter::GetFitter()->GetParameter(1); 
          ee_da_p_e = TVirtualFitter::GetFitter()->GetParError(1);  
          ee_da_p0   = TVirtualFitter::GetFitter()->GetParameter(0); 
          ee_da_p0_e = TVirtualFitter::GetFitter()->GetParError(0);  
          hmcdyee->Fit("f1", "R");                                       
          ee_mc_p   = TVirtualFitter::GetFitter()->GetParameter(1); 
          ee_mc_p_e = TVirtualFitter::GetFitter()->GetParError(1);  
          ee_mc_p0   = TVirtualFitter::GetFitter()->GetParameter(0); 
          ee_mc_p0_e = TVirtualFitter::GetFitter()->GetParError(0);  
          hdatadygg->Fit("f1", "R");                                      
          gg_da_p   = TVirtualFitter::GetFitter()->GetParameter(1); 
          gg_da_p_e = TVirtualFitter::GetFitter()->GetParError(1);  
          gg_da_p0   = TVirtualFitter::GetFitter()->GetParameter(0); 
          gg_da_p0_e = TVirtualFitter::GetFitter()->GetParError(0);  
          hmcdygg->Fit("f1", "R");                                       
          gg_mc_p   = TVirtualFitter::GetFitter()->GetParameter(1); 
          gg_mc_p_e = TVirtualFitter::GetFitter()->GetParError(1);  
          gg_mc_p0   = TVirtualFitter::GetFitter()->GetParameter(0); 
          gg_mc_p0_e = TVirtualFitter::GetFitter()->GetParError(0); 
          mm_damc_p0_e = (sqrt((mm_da_p0_e/mm_da_p0)*(mm_da_p0_e/mm_da_p0) +(mm_mc_p0_e/mm_mc_p0)*(mm_mc_p0_e/mm_mc_p0)))*(mm_da_p0/mm_mc_p0); 
          mm_damc_p_e = (sqrt((mm_da_p_e/mm_da_p)*(mm_da_p_e/mm_da_p) +(mm_mc_p_e/mm_mc_p)*(mm_mc_p_e/mm_mc_p)))*(mm_da_p/mm_mc_p); 
          ee_damc_p0_e = (sqrt((ee_da_p_e/ee_da_p0)*(ee_da_p0_e/ee_da_p0) +(ee_mc_p0_e/ee_mc_p0)*(ee_mc_p0_e/ee_mc_p0)))*(ee_da_p0/ee_mc_p0); 
          ee_damc_p_e = (sqrt((ee_da_p_e/ee_da_p)*(ee_da_p_e/ee_da_p) +(ee_mc_p_e/ee_mc_p)*(ee_mc_p_e/ee_mc_p)))*(ee_da_p/ee_mc_p); 
          gg_damc_p0_e = (sqrt((gg_da_p0_e/gg_da_p0)*(gg_da_p0_e/gg_da_p0) +(gg_mc_p0_e/gg_mc_p0)*(gg_mc_p0_e/gg_mc_p0)))*(gg_da_p0/gg_mc_p0); 
          gg_damc_p_e = (sqrt((gg_da_p_e/gg_da_p)*(gg_da_p_e/gg_da_p) +(gg_mc_p_e/gg_mc_p)*(gg_mc_p_e/gg_mc_p)))*(gg_da_p/gg_mc_p);                 
          
          if(histogramaname=="met_uPara_vs_nVert" or histogramaname == "met_uParaPuppi__vs_nVert"){cout << " uPara " << endl;}
          else{cout << " uPerp " << endl;}                                              
          cout<<"\n Z$\\rightarrow \\mu\\mu$ & "<< mm_da_p0<<" $\\pm$ "<<mm_da_p0_e<<" & " <<mm_mc_p0<<" $\\pm$ "<<mm_mc_p0_e <<" & "<< mm_da_p<<" $\\pm$ "<<mm_da_p_e<<" & " <<mm_da_p/mm_mc_p<<" $\\pm$ "<< mm_damc_p_e<<"\\\\"<<endl; 
          cout<<"\n Z$\\rightarrow$ ee       & "<< ee_da_p0<<" $\\pm$ "<<ee_da_p0_e<<" & " <<ee_mc_p0<<" $\\pm$ "<<ee_mc_p0_e <<" & "<< ee_da_p<<" $\\pm$ "<<ee_da_p_e<<" & " <<ee_da_p/ee_mc_p<<" $\\pm$ "<< ee_damc_p_e<<"\\\\"<<endl; 
          if(histogramaname=="met_uPara__vs_nVert" or histogramaname == "met_uPerp__vs_nVert"){
          cout<<"\n $\\gamma$+jets           & "<< gg_da_p0<<" $\\pm$ "<<gg_da_p0_e<<" & " <<gg_mc_p0<<" $\\pm$ "<<gg_mc_p0_e <<" & "<< gg_da_p<<" $\\pm$ "<<gg_da_p_e<<" & " <<gg_da_p/gg_mc_p<<" $\\pm$ "<< gg_damc_p_e<<"\\\\"<<endl; 
          }
          cout<<"\n"<< endl; 
   }
    leg->Draw();
    c1->Modified();
    c1->cd();
 
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.04);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();
    
    ratiodyee->GetYaxis()->SetTitle("Data / MC");
    ratiodyee->GetXaxis()->SetTitle(EjeX);
    ratiodyee->GetYaxis()->SetRangeUser(0.6, 1.4);
    if (histogramaname.Contains("over")) 
        ratiodyee->GetYaxis()->SetRangeUser(0.9,1.1);
    if (Type == "tk")ratiodyee->GetYaxis()->SetRangeUser(0.8,1.2);
    ratiodyee->SetMarkerStyle(22);
    ratiodyee->SetMarkerColor(46);
    ratiodyee->SetLineColor(46);
    ratiodymm->SetMarkerStyle(23);
    ratiodymm->SetMarkerColor(4);
    ratiodymm->SetLineColor(4); 
    ratiodygg->SetMarkerStyle(20);
    ratiodygg->SetMarkerColor(8);
    ratiodygg->SetLineColor(8);     
    ratiodyee->GetXaxis()->SetLabelFont(43); 
    ratiodyee->GetXaxis()->SetLabelSize(16); 
    ratiodyee->GetYaxis()->SetLabelFont(43);
    ratiodyee->GetYaxis()->SetLabelSize(16);
    ratiodyee->GetXaxis()->SetTitleFont(43); 
    ratiodyee->GetXaxis()->SetTitleSize(25); 
    ratiodyee->GetYaxis()->SetTitleFont(43);
    ratiodyee->GetYaxis()->SetTitleSize(25);
    ratiodyee->GetXaxis()->SetTitleOffset(3);
    ratiodyee->GetYaxis()->SetTitleOffset(1.2);
    ratiodyee->SetFillColor(54);
    ratiodyee->SetFillStyle(3345);                                                                
  
    //==========  JES ============================================                                                                                                            
    
    //============================================================
    TGraphAsymmErrors tottoterree= TOTTOTerror(hmcdyee,hmcdyeedown,hmcdyeeup, hmcdyeejerdown,hmcdyeejerup,hmcdyeedownUncl,hmcdyeeupUncl,hdatadyee);
    TGraphAsymmErrors jererree= JERerror(hmcdyee,hmcdyeedown,hmcdyeeup, hmcdyeejerdown,hmcdyeejerup,hdatadyee);
    TGraphAsymmErrors erree= JESerror(hmcdyee,hmcdyeedown,hmcdyeeup,hdatadyee);
    TGraphAsymmErrors staterree= staterror(hmcdyee,hdatadyee);

   TAxis *erree_yaxis = erree.GetYaxis ();
   erree_yaxis->SetRangeUser (0, 2);
   erree.SetTitle (0);
   erree.SetFillColor (kGreen);
   erree.SetFillStyle (3002);                                   
   
   TAxis *staterree_yaxis = staterree.GetYaxis ();
   staterree_yaxis->SetRangeUser (0, 2);
   staterree.SetTitle (0);
   staterree.SetFillColor (kRed+1);
   staterree.SetFillStyle (3002);  
   
   TAxis *jererree_yaxis = jererree.GetYaxis ();
   jererree_yaxis->SetRangeUser (0, 2);
   jererree.SetTitle (0);
   jererree.SetFillColor (kBlue+2);
   jererree.SetFillStyle (3002);                                   
    
   TAxis *tottoterree_yaxis = tottoterree.GetYaxis ();       
   tottoterree_yaxis->SetRangeUser (0, 2);
   tottoterree.SetTitle (0);
   tottoterree.SetFillColor (kGray+1);
   tottoterree.SetFillStyle (3002);                      


   ratiodyee->Draw("e1,X0 "); //"e2"
   ratiodyee->Draw("AXIS"); //"e2"
   ratiodygg->Draw("e1,X0  same "); //"e2"
   ratiodyee->Draw("e1,X0  same"); //"e2"
   ratiodymm->Draw("e1,X0  same"); //"e2"
   tottoterree.Draw ("2 same");
   //jererree.Draw ("2 same");
   //erree.Draw ("2 same");
   //staterree.Draw("2 same");
   //ratiodygg->Draw("e1,X0  same "); //"e2"
   ratiodymm->Draw("e1,X0  same"); //"e2"
   TLine *lineR =  new TLine ( ratiodyee->GetXaxis ()->GetXmin (), 1, ratiodyee->GetXaxis ()->GetXmax (), 1);
   lineR->SetLineColor (kBlack ); lineR->SetLineWidth (2); lineR->SetLineStyle (2); lineR->Draw();     
        
    TLegend* legratio(0);
   
    if (Type == "tk" or Type == "tkRes"){
        legratio = new TLegend(0.77,0.78,0.775,0.785);
    }
    else{
        legratio = new TLegend(0.27,0.78,0.87,0.95);
    }
    legratio->SetFillColor(0);
    legratio->SetBorderSize(0);
    legratio-> SetNColumns(2);
    if (Type == "tk" or Type == "tkRes"){
        legratio->Draw();
    }
    else{
    legratio->AddEntry(&tottoterree, "Uncl p_{T}^{miss} + JER+ JES + Stat","f");
    legratio->AddEntry(&jererree, "JER + JES + Stat","f");
    legratio->AddEntry(&erree, "JES + Stat","f");
    legratio->AddEntry(&staterree, "Stat","f");
    //legratio->Draw();
    }
    //    ratiodygg->Draw("e1,  same "); //"e2"
//    ratiodymm->Draw("e1,  same"); //"e2"
//    ratiodyee->Draw("e1,  same"); //"e2"

    pad1->cd();
    TLatex latex;                  
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.08);
    latex.DrawLatex(0.26, 0.80, "#bf{CMS}");

//    if ((histogramaname.Contains("sumEt"))  or (histogramaname.Contains("nVert")) ) {
        TLatex latexz;
        latexz.SetNDC();
        latexz.SetTextAngle(0);
        latexz.SetTextFont(42);
        latexz.SetTextAlign(31);
        latexz.SetTextSize(0.05);
   if  (histogramaname.Contains("vs")){ 
        latexz.DrawLatex(0.6, 0.28, "Response corrected");
    }
    
                                                    
    TLatex latexc;
    latexc.SetNDC();
    latexc.SetTextAngle(0);
    latexc.SetTextFont(42);
    latexc.SetTextAlign(31);
    latexc.SetTextSize(0.06);
    latexc.DrawLatex(0.90, 0.91, "35.9 fb^{-1} (13 TeV)");              
    TLatex latexd;
    latexd.SetNDC();
    latexd.SetTextAngle(0);
    latexd.SetTextFont(42);
    latexd.SetTextAlign(31);
    latexd.SetTextSize(0.06);
    latexd.SetTextAngle(90); 
    latexd.DrawLatex(0.038, 0.87, ylabel);      



  
    pad1->Modified();
    pad2->Modified();
    pad1->Update();
    pad2->Update();
    pad1->cd();
    c1->Update();
    c1->cd();
    c1->Print("~/www/met/paper/comparisons/June19/"+histogramaname+Type+".png");
    c1->Print("~/www/met/paper/comparisons/June19/"+histogramaname+Type+".root");
    c1->Print("~/www/met/paper/comparisons/June19/"+histogramaname+Type+".pdf");  


}


void Comparison(TString histograma, TString histograma2, TString EjeX, TString Type) {                                               

TCanvas *c1 = new TCanvas("c1","example",600,700);//Feb29
TFile file1EE ("../tgraphs/DYEZllResolutionMay25sumEtRMS.root"); //Mar27_qT50
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25sumEtRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25sumEtRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay25sumEtRMS.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25sumEtRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay25sumEtRMS.root");
 TFile file1MCEEup        ("../tgraphs/DY_up_jes_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCEEdown      ("../tgraphs/DY_down_jes_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCEEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCEEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCEEupUncl    ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25sumEtRMS.root");
 TFile file1MCEEdownUncl  ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25sumEtRMS.root");
 TFile file1MCMMup        ("../tgraphs/DY_up_jes_DYMZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCMMdown      ("../tgraphs/DY_down_jes_DYMZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCMMjerup     ("../tgraphs/DY_up_jer_DYMZllResolutionMay25sumEtRMS.root"); 
 TFile file1MCMMjerdown   ("../tgraphs/DY_down_jer_DYMZllResolutionMay25sumEtRMS.root");
 TFile file1MCMMupUncl    ("../tgraphs/DY_up_uncl_DYMZllResolutionMay25sumEtRMS.root");
 TFile file1MCMMdownUncl  ("../tgraphs/DY_down_uncl_DYMZllResolutionMay25sumEtRMS.root");
 TFile file1MCGGup        ("../tgraphs/GJets_up_jes_DYgjetsResolutionMay25sumEtRMS.root"); 
 TFile file1MCGGdown      ("../tgraphs/GJets_down_jes_DYgjetsResolutionMay25sumEtRMS.root"); 
 TFile file1MCGGjerup     ("../tgraphs/GJets_up_jer_DYgjetsResolutionMay25sumEtRMS.root"); 
 TFile file1MCGGjerdown   ("../tgraphs/GJets_down_jer_DYgjetsResolutionMay25sumEtRMS.root");
 TFile file1MCGGupUncl    ("../tgraphs/GJets_up_uncl_DYgjetsResolutionMay25sumEtRMS.root");
 TFile file1MCGGdownUncl  ("../tgraphs/GJets_down_uncl_DYgjetsResolutionMay25sumEtRMS.root"); 

// TFile file1EE ("../tgraphs/DYEZllResolutionMay25nVertRMS.root");
// TFile file2EE ("../tgraphs/DataEZllResolutionMay25nVertRMS.root");
// TFile file1MM ("../tgraphs/DYMZllResolutionMay25nVertRMS.root"); 
// TFile file2MM ("../tgraphs/DataMZllResolutionMay25nVertRMS.root");
// TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25nVertRMS.root"); 
// TFile file2GG ("../tgraphs/DatagjetsResolutionMay25nVertRMS.root");
// TFile file1MCEEup        ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
// TFile file1MCEEdown      ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
// TFile file1MCEEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
// TFile file1MCEEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
// TFile file1MCEEupUncl    ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
// TFile file1MCEEdownUncl  ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
// TFile file1MCMMup        ("../tgraphs/DY_up_jes_DYMZllResolutionMay25nVertRMS.root"); 
// TFile file1MCMMdown      ("../tgraphs/DY_down_jes_DYMZllResolutionMay25nVertRMS.root"); 
// TFile file1MCMMjerup     ("../tgraphs/DY_up_jer_DYMZllResolutionMay25nVertRMS.root"); 
// TFile file1MCMMjerdown   ("../tgraphs/DY_down_jer_DYMZllResolutionMay25nVertRMS.root");
// TFile file1MCMMupUncl    ("../tgraphs/DY_up_uncl_DYMZllResolutionMay25nVertRMS.root");
// TFile file1MCMMdownUncl  ("../tgraphs/DY_down_uncl_DYMZllResolutionMay25nVertRMS.root");
// TFile file1MCGGup        ("../tgraphs/GJets_up_jes_DYgjetsResolutionMay25nVertRMS.root"); 
// TFile file1MCGGdown      ("../tgraphs/GJets_down_jes_DYgjetsResolutionMay25nVertRMS.root"); 
// TFile file1MCGGjerup     ("../tgraphs/GJets_up_jer_DYgjetsResolutionMay25nVertRMS.root"); 
// TFile file1MCGGjerdown   ("../tgraphs/GJets_down_jer_DYgjetsResolutionMay25nVertRMS.root");
// TFile file1MCGGupUncl    ("../tgraphs/GJets_up_uncl_DYgjetsResolutionMay25nVertRMS.root");
// TFile file1MCGGdownUncl  ("../tgraphs/GJets_down_uncl_DYgjetsResolutionMay25nVertRMS.root"); 


 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 hMCEEup = (TGraphErrors*)   file1MCEEup.Get("E_"+ histograma); 
 hMCEEdown = (TGraphErrors*) file1MCEEdown.Get("E_"+ histograma); 
 hMCEEjerup = (TGraphErrors*)   file1MCEEjerup.Get("E_"+ histograma); 
 hMCEEjerdown = (TGraphErrors*) file1MCEEjerdown.Get("E_"+ histograma); 
 hMCEEupUncl = (TGraphErrors*)   file1MCEEupUncl.Get("E_"+ histograma); 
 hMCEEdownUncl = (TGraphErrors*)   file1MCEEdownUncl.Get("E_"+ histograma);                 
 hMCMMup = (TGraphErrors*)   file1MCMMup.Get("M_"+ histograma); 
 hMCMMdown = (TGraphErrors*) file1MCMMdown.Get("M_"+ histograma); 
 hMCMMjerup = (TGraphErrors*)   file1MCMMjerup.Get("M_"+ histograma); 
 hMCMMjerdown = (TGraphErrors*) file1MCMMjerdown.Get("M_"+ histograma); 
 hMCMMupUncl = (TGraphErrors*)   file1MCMMupUncl.Get("M_"+ histograma); 
 hMCMMdownUncl = (TGraphErrors*)   file1MCMMdownUncl.Get("M_"+ histograma); 
 hMCGGup = (TGraphErrors*)   file1MCGGup.Get(histograma2); 
 hMCGGdown = (TGraphErrors*) file1MCGGdown.Get(histograma2); 
 hMCGGjerup = (TGraphErrors*)   file1MCGGjerup.Get(histograma2); 
 hMCGGjerdown = (TGraphErrors*) file1MCGGjerdown.Get(histograma2); 
 hMCGGupUncl = (TGraphErrors*)   file1MCGGupUncl.Get(histograma2); 
 hMCGGdownUncl = (TGraphErrors*)   file1MCGGdownUncl.Get(histograma2);   




    cout << " h2EE " << h2EE << endl ;            
    cout << " h1MM " << h1MM << endl ;               
    cout << " h2MM " << h2MM << endl ;               
    cout << " h1GG " << h1GG << endl ;               
    cout << " h2GG " << h2GG << endl ;               
    cout << " h1EEup " << hMCEEup << endl ;
    cout << " h1EEdown " << hMCEEdown << endl ;
    cout << " hMCGGup "        << hMCGGup << endl ;
    cout << " hMCGGdown "      << hMCGGdown << endl ;
    cout << " hMCGGjerup "     << hMCGGjerup << endl ;
    cout << " hMCGGjerdown "   << hMCGGjerdown << endl ;
    cout << " hMCGGupUncl "    << hMCGGupUncl << endl ;
    cout << " hMCGGdownUncl "  << hMCGGdownUncl << endl ;

   gStyle->SetOptStat(0);
//   gStyle->SetOptFit(1111);
   gStyle->SetErrorX(0.5);     

   PlotWithRatioEEMuMuGJetsComparison(c1, h1EE, h2EE, h1MM, h2MM, h1GG,h2GG, hMCEEup, hMCEEdown,  hMCEEjerup, hMCEEjerdown, hMCEEupUncl, hMCEEdownUncl,hMCMMup, hMCMMdown,  hMCMMjerup, hMCMMjerdown, hMCMMupUncl, hMCMMdownUncl, hMCGGup, hMCGGdown,  hMCGGjerup, hMCGGjerdown, hMCGGupUncl, hMCGGdownUncl, EjeX, histograma, Type);
}


TGraphAsymmErrors JESerror(TH1F *Background, TH1F *Backgrounddown, TH1F *Backgroundup, TH1F * hdata){
    
    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");
	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	
    const Int_t nvar = 27;
	Double_t x[nvar];
	Double_t y[nvar];
	Double_t exl[nvar];
	Double_t eyl[nvar];
	Double_t exh[nvar];
	Double_t eyh[nvar];
	Double_t x1[nvar];
	Double_t y1[nvar];
	Double_t exl1[nvar];
	Double_t eyl1[nvar];
	Double_t exh1[nvar];
	Double_t eyh1[nvar];

    TH1F *ratiop = (TH1F *) Background->Clone ("backgroundratiop");
    TH1F *ratiom = (TH1F *) Background->Clone ("backgroundratiom");

	double  ymax = 2.;

	for (int km = 0; km <= Background->GetNbinsX (); km++){
        double conte1 =
            sqrt (Background->GetBinError (km) *
	        Background->GetBinError (km) +
	        (Background->GetBinContent (km) -
	        Backgroundup->GetBinContent (km)) *
	        (Background->GetBinContent (km) -
	        Backgroundup->GetBinContent (km)));
        double conte2 =
            sqrt (Background->GetBinError (km) *
	        Background->GetBinError (km) +
	        (Background->GetBinContent (km) -
	        Backgrounddown->GetBinContent (km)) *
	        (Background->GetBinContent (km) -
	        Backgrounddown->GetBinContent (km)));


	    den1->SetBinContent (km,
	  		   Background->GetBinContent (km) + conte1);
	    den2->SetBinContent (km,
	  		   Background->GetBinContent (km) - conte2);
	          ymax = Background->GetBinContent(km) + conte1;
	    x1[km] = Background->GetBinCenter (km);
	    y1[km] = Background->GetBinContent (km);
	    exl1[km] = Background->GetBinWidth (km) / 2;
	    exh1[km] = Background->GetBinWidth (km) / 2;
	    eyl1[km] = conte2;
	    eyh1[km] = conte1;
	  }

	  ratiop->Divide (den1);
	  ratiom->Divide (den2);
	  
	  TH1F *ratio = (TH1F *) hdata->Clone ("ratiodata");
	  ratio->Divide (Background);

	  for (int km = 0; km <= ratio->GetNbinsX (); km++){
	      if (ratio->GetBinContent (km) > ymax)
            ymax = ratio->GetBinContent (km) + ratio->GetBinError (km);
	      x[km] = ratio->GetBinCenter (km);
	      y[km] = 1;	
	      exl[km] = ratio->GetBinWidth (km) / 2;
	      exh[km] = ratio->GetBinWidth (km) / 2;


	      if (ratiop->GetBinContent (km) != 0)
		eyh[km] = (1. / ratiop->GetBinContent (km) - 1)*ratio->GetBinContent (km);
	      else
		eyh[km] = 0;

	      if (ratiom->GetBinContent (km) != 0)
		eyl[km] = (1 - 1. / ratiom->GetBinContent (km))*ratio->GetBinContent (km);
	      else
		eyl[km] = 0;
          

	    }
    
TGraphAsymmErrors *err = new TGraphAsymmErrors (nvar, x, y, exl, exh, eyl, eyh);


 return *err;


}

TGraphAsymmErrors staterror(TH1F *Background, TH1F * hdata){

   

    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");

	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	const Int_t nvar=27;// (const Int_t)Background->GetNbinsX ();
	Double_t x[nvar];
	Double_t y[nvar];
	Double_t exl[nvar];
	Double_t eyl[nvar];
	Double_t exh[nvar];
	Double_t eyh[nvar];
	Double_t x1[nvar];
	Double_t y1[nvar];
	Double_t exl1[nvar];
	Double_t eyl1[nvar];
	Double_t exh1[nvar];
	Double_t eyh1[nvar];



    TH1F *ratiop = (TH1F *) Background->Clone ("backgroundratiop");
    TH1F *ratiom = (TH1F *) Background->Clone ("backgroundratiom");

	  double  ymax = 2.;


	  for (int km = 0; km <= Background->GetNbinsX (); km++)
	    {

	      double conte1 = Background->GetBinError (km) ;
		      
	      double conte2 = Background->GetBinError (km);
		      


	      den1->SetBinContent (km,
				   Background->GetBinContent (km) + conte1);
	      den2->SetBinContent (km,
				   Background->GetBinContent (km) - conte2);


	            ymax = Background->GetBinContent(km) + conte1;
	      x1[km] = Background->GetBinCenter (km);
	      y1[km] = Background->GetBinContent (km);
	      exl1[km] = Background->GetBinWidth (km) / 2;
	      exh1[km] = Background->GetBinWidth (km) / 2;
	      eyl1[km] = conte2;
	      eyh1[km] = conte1;
	    }


	  ratiop->Divide (den1);
	  ratiom->Divide (den2);
	  
	  TH1F *ratio = (TH1F *) hdata->Clone ("ratiodata");
	  ratio->Divide (Background);

	  for (int km = 0; km <= ratio->GetNbinsX (); km++)
	    {
	      if (ratio->GetBinContent (km) > ymax)
     		ymax = ratio->GetBinContent (km) + ratio->GetBinError (km);
	      x[km] = ratio->GetBinCenter (km);
	      y[km] = 1;	
	      exl[km] = ratio->GetBinWidth (km) / 2;
	      exh[km] = ratio->GetBinWidth (km) / 2;


	      if (ratiop->GetBinContent (km) != 0)
		eyh[km] =( 1. / ratiop->GetBinContent (km) - 1)*ratio->GetBinContent (km);
	      else
		eyh[km] = 0;

	      if (ratiom->GetBinContent (km) != 0)
		eyl[km] = (1 - 1. / ratiom->GetBinContent (km))*ratio->GetBinContent (km);
	      else
		eyl[km] = 0;


	    }


                                                                                                                                          
TGraphAsymmErrors *err = new TGraphAsymmErrors (nvar, x, y, exl, exh, eyl, eyh);


 return *err;


}


TGraphAsymmErrors JERerror(TH1F *Background, TH1F *Backgrounddown, TH1F *Backgroundup,  TH1F *BackgroundUncldown, TH1F *BackgroundUnclup, TH1F * hdata){
    
    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");
	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	
    const Int_t nvar = 27;
	Double_t x[nvar];
	Double_t y[nvar];
	Double_t exl[nvar];
	Double_t eyl[nvar];
	Double_t exh[nvar];
	Double_t eyh[nvar];
	Double_t x1[nvar];
	Double_t y1[nvar];
	Double_t exl1[nvar];
	Double_t eyl1[nvar];
	Double_t exh1[nvar];
	Double_t eyh1[nvar];

    TH1F *ratiop = (TH1F *) Background->Clone ("backgroundratiop");
    TH1F *ratiom = (TH1F *) Background->Clone ("backgroundratiom");

	double  ymax = 2.;

	for (int km = 0; km <= Background->GetNbinsX (); km++){
        double conte1 =
            sqrt (Background->GetBinError (km) *Background->GetBinError (km) + (Background->GetBinContent (km) -BackgroundUnclup->GetBinContent (km)) *(Background->GetBinContent (km) -BackgroundUnclup->GetBinContent (km))+ (Background->GetBinContent (km) -Backgroundup->GetBinContent (km)) *(Background->GetBinContent (km) -Backgroundup->GetBinContent (km)));
        double conte2 =
            sqrt (Background->GetBinError (km) *Background->GetBinError (km) + (Background->GetBinContent (km) -BackgroundUncldown->GetBinContent (km)) *(Background->GetBinContent (km) -BackgroundUncldown->GetBinContent (km)) +(Background->GetBinContent (km) -Backgrounddown->GetBinContent (km)) *(Background->GetBinContent (km) -Backgrounddown->GetBinContent (km)));

	    if (conte1> conte2){
        den1->SetBinContent (km,
	  		   Background->GetBinContent (km) + conte1);
	    den2->SetBinContent (km,
	  		   Background->GetBinContent (km) - conte1);
        }
        else{
        den1->SetBinContent (km,
        	   Background->GetBinContent (km) + conte2);
        den2->SetBinContent (km,
        	   Background->GetBinContent (km) - conte2);
        }
	    ymax = Background->GetBinContent(km) + conte1;
	    x1[km] = Background->GetBinCenter (km);
	    y1[km] = Background->GetBinContent (km);
	    exl1[km] = Background->GetBinWidth (km) / 2;
	    exh1[km] = Background->GetBinWidth (km) / 2;
	    eyl1[km] = conte2;
	    eyh1[km] = conte1;
	  }

	  ratiop->Divide (den1);
	  ratiom->Divide (den2);
	  
	  TH1F *ratio = (TH1F *) hdata->Clone ("ratiodata");
	  ratio->Divide (Background);

	  for (int km = 0; km <= ratio->GetNbinsX (); km++){
	      if (ratio->GetBinContent (km) > ymax)
            ymax = ratio->GetBinContent (km) + ratio->GetBinError (km);
	      x[km] = ratio->GetBinCenter (km);
	      y[km] = 1;	
	      exl[km] = ratio->GetBinWidth (km) / 2;
	      exh[km] = ratio->GetBinWidth (km) / 2;


	      if (ratiop->GetBinContent (km) != 0)
		eyh[km] = (1. / ratiop->GetBinContent (km) - 1)*ratio->GetBinContent (km);
	      else
		eyh[km] = 0;

	      if (ratiom->GetBinContent (km) != 0)
		eyl[km] = (1 - 1. / ratiom->GetBinContent (km))*ratio->GetBinContent (km);
	      else
		eyl[km] = 0;
          

	    }                                                                                                                                                                                           
    
TGraphAsymmErrors *err = new TGraphAsymmErrors (nvar, x, y, exl, exh, eyl, eyh);


 return *err;

}                                                                                                                                         
                                                                                                                                                                                                                                                                                 


TGraphAsymmErrors TOTTOTerror(TH1F *Background, TH1F *Backgrounddown, TH1F *Backgroundup, TH1F *Backgroundjerdown, TH1F *Backgroundjerup, TH1F *BackgroundUncldown, TH1F *BackgroundUnclup, TH1F * hdata){
    
    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");
	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	
    const Int_t nvar = 27;
	Double_t x[nvar];
	Double_t y[nvar];
	Double_t exl[nvar];
	Double_t eyl[nvar];
	Double_t exh[nvar];
	Double_t eyh[nvar];
	Double_t x1[nvar];
	Double_t y1[nvar];
	Double_t exl1[nvar];
	Double_t eyl1[nvar];
	Double_t exh1[nvar];
	Double_t eyh1[nvar];

    TH1F *ratiop = (TH1F *) Background->Clone ("backgroundratiop");
    TH1F *ratiom = (TH1F *) Background->Clone ("backgroundratiom");

	double  ymax = 2.;

	for (int km = 0; km <= Background->GetNbinsX (); km++){
            
        double conte1 =
            sqrt (Background->GetBinError (km) *Background->GetBinError (km) + (Background->GetBinContent (km) -BackgroundUnclup->GetBinContent (km)) *(Background->GetBinContent (km) -BackgroundUnclup->GetBinContent (km))+ (Background->GetBinContent (km) -Backgroundup->GetBinContent (km)) *(Background->GetBinContent (km) -Backgroundup->GetBinContent (km))+ (Background->GetBinContent (km) -Backgroundjerup->GetBinContent (km)) *(Background->GetBinContent (km) -Backgroundjerup->GetBinContent (km)));
        double conte2 =
            sqrt (Background->GetBinError (km) *Background->GetBinError (km) + (Background->GetBinContent (km) -BackgroundUncldown->GetBinContent (km)) *(Background->GetBinContent (km) -BackgroundUncldown->GetBinContent (km)) +(Background->GetBinContent (km) -Backgrounddown->GetBinContent (km)) *(Background->GetBinContent (km) -Backgrounddown->GetBinContent (km))+(Background->GetBinContent (km) -Backgroundjerdown->GetBinContent (km)) *(Background->GetBinContent (km) -Backgroundjerdown->GetBinContent (km))  );
        if  (abs(conte1-conte2)/(conte1+conte2) <1){
	        if (conte1> conte2){
                den1->SetBinContent (km,Background->GetBinContent (km) + conte1);
	            den2->SetBinContent (km,Background->GetBinContent (km) - conte1);
            }
            if (conte2> conte1){
                den1->SetBinContent (km,Background->GetBinContent (km) + conte2);
	            den2->SetBinContent (km,Background->GetBinContent (km) - conte2);                     
            }
        }
        else{
            if (conte1> conte2){
                den1->SetBinContent (km,Background->GetBinContent (km) + conte2);
                den2->SetBinContent (km,Background->GetBinContent (km) - conte2);                     
            }                                                                                                 
            else{
                den1->SetBinContent (km,Background->GetBinContent (km) + conte1);
                den2->SetBinContent (km,Background->GetBinContent (km) - conte1);  
            } 
        }
        
        ymax = Background->GetBinContent(km) + conte1;
	    x1[km] = Background->GetBinCenter (km);
	    y1[km] = Background->GetBinContent (km);
	    exl1[km] = Background->GetBinWidth (km) / 2;
	    exh1[km] = Background->GetBinWidth (km) / 2;
	    eyl1[km] = conte2;
	    eyh1[km] = conte1;
	  }

	  ratiop->Divide (den1);
	  ratiom->Divide (den2);
	  
	  TH1F *ratio = (TH1F *) hdata->Clone ("ratiodata");
	  ratio->Divide (Background);

	  for (int km = 0; km <= ratio->GetNbinsX (); km++){
	      if (ratio->GetBinContent (km) > ymax)
            ymax = ratio->GetBinContent (km) + ratio->GetBinError (km);
	      x[km] = ratio->GetBinCenter (km);
	      y[km] = 1;	
	      exl[km] = ratio->GetBinWidth (km) / 2;
	      exh[km] = ratio->GetBinWidth (km) / 2;


	      if (ratiop->GetBinContent (km) != 0)
		eyh[km] = (1. / ratiop->GetBinContent (km) - 1)*ratio->GetBinContent (km);
	      else
		eyh[km] = 0;

	      if (ratiom->GetBinContent (km) != 0)
		eyl[km] = (1 - 1. / ratiom->GetBinContent (km))*ratio->GetBinContent (km);
	      else
		eyl[km] = 0;
          

	    }                                                                                                                                                                                           
    
TGraphAsymmErrors *err = new TGraphAsymmErrors (nvar, x, y, exl, exh, eyl, eyh);


 return *err;


}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   



TGraphAsymmErrors TOTNEWerror(TH1F *h, TH1F *hdown, TH1F *hup, TH1F *hjerdown, TH1F *hjerup, TH1F *hUncldown, TH1F *hUnclup){
    
	
    const Int_t nvar = 27;
	Double_t x[nvar];
	Double_t y[nvar];
	Double_t exl[nvar];
	Double_t eyl[nvar];
	Double_t exh[nvar];
	Double_t eyh[nvar];
	Double_t x1[nvar];
	Double_t y1[nvar];

	double  ymax = 2.;

	for (int km = 0; km <= h->GetNbinsX (); km++){
        double conte1 =
            sqrt ( (h->GetBinContent (km) -hUnclup->GetBinContent (km)) *(h->GetBinContent (km) -hUnclup->GetBinContent (km))+ (h->GetBinContent (km) -hup->GetBinContent (km)) *(h->GetBinContent (km) -hup->GetBinContent (km))+ (h->GetBinContent (km) -hjerup->GetBinContent (km)) *(h->GetBinContent (km) -hjerup->GetBinContent (km)));
        //cout << "bin    " << km << "***********************************" << endl;
        //cout << "jes Up " << (h->GetBinContent (km) -hup->GetBinContent (km)) << endl;
        //cout << "jes Dn " << (h->GetBinContent (km) -hdown->GetBinContent (km)) << endl;
        //cout << "jer Up " << (h->GetBinContent (km) -hjerup->GetBinContent (km)) << endl;
        //cout << "jer Dn " << (h->GetBinContent (km) -hjerdown->GetBinContent (km)) << endl;
        //cout << "unc Up " << (h->GetBinContent (km) -hUnclup->GetBinContent (km)) << endl;
        //cout << "unc Dn " << (h->GetBinContent (km) -hUncldown->GetBinContent (km)) << endl;
        
        
        
        
        
        
        
        ymax = h->GetBinContent(km) + conte1;
	    x1[km] = h->GetBinCenter (km);
	    y1[km] = h->GetBinContent (km);
	    exl[km] = h->GetBinWidth (km) / 2;
	    exh[km] = h->GetBinWidth (km) / 2;
	    eyl[km] = conte1;
	    eyh[km] = conte1;
	  }

      //ratiop->Divide (den1);
	  //ratiom->Divide (den2);
	  
	  //
//      TH1F *ratio = (TH1F *) h->Clone ("ratiodata");
//	  //ratio->Divide (h);
//
//	  for (int km = 0; km <= ratio->GetNbinsX (); km++){
//	      if (ratio->GetBinContent (km) > ymax)
//            ymax = ratio->GetBinContent (km) + ratio->GetBinError (km);
//	      x[km] = ratio->GetBinCenter (km);
//	      y[km] = 1;	
//	      exl[km] = ratio->GetBinWidth (km) / 2;
//	      exh[km] = ratio->GetBinWidth (km) / 2;
//
//
//	      if (ratiop->GetBinContent (km) != 0)
//		eyh[km] = (ratiop->GetBinContent (km) - ratio->GetBinContent (km));
//	      else
//		eyh[km] = 0;
//
//	      if (ratiom->GetBinContent (km) != 0)
//		eyl[km] = (ratio->GetBinContent (km) - ratiom->GetBinContent (km));
//	      else
//		eyl[km] = 0;
//          
//
//	    }                                                                                                                                                                                           
//    
TGraphAsymmErrors *err = new TGraphAsymmErrors (nvar, x, y, exl, exh, eyl, eyh);


 return *err;


}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
