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
TGraphErrors *h1EEup = new TGraphErrors;
TGraphErrors *h1EEdown = new TGraphErrors;
TGraphErrors *h1EEjerup = new TGraphErrors;
TGraphErrors *h1EEjerdown = new TGraphErrors;
TGraphErrors *h1EEupUncl = new TGraphErrors;
TGraphErrors *h1EEdownUncl = new TGraphErrors;
TGraphErrors *h1MM = new TGraphErrors;
TGraphErrors *h2MM = new TGraphErrors;
TGraphErrors *h1MMup = new TGraphErrors;
TGraphErrors *h1MMdown = new TGraphErrors;
TGraphErrors *h1GG = new TGraphErrors;
TGraphErrors *h2GG = new TGraphErrors;
TGraphErrors *h1GGup = new TGraphErrors;
TGraphErrors *h1GGdown = new TGraphErrors;
void testRatioEEMuMuGJetsComparison(TString, TString, TString);                                            
//TMultiGraph mg;
Int_t doErrorLegend;
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

void PlotWithRatioEEMuMuGJetsComparison(TCanvas *c1, TGraphErrors *mcdyee,  TGraphErrors *datadyee,  TGraphErrors *mcdymm,  TGraphErrors *datadymm,TGraphErrors *mcdygg,  TGraphErrors *datadygg,  TGraphErrors *mcdyeeup, TGraphErrors *mcdyeedown,  TGraphErrors *mcdyeejerup, TGraphErrors *mcdyeejerdown, TGraphErrors *mcdyeeupUncl, TGraphErrors *mcdyeedownUncl,  TString EjeX, TString histogramaname, TString Type) {

TH1::SetDefaultSumw2();
  c1->cd(); 
//  Float_t xbins[] = { 4, 6, 8, 10, 12,14, 16, 18, 20,22, 24, 26, 28, 32,  40};
//double nbins=14.;
//this is the OG one Float_t xbins[] = {400,500,600, 700,800, 900, 1000, 1100, 1300, 1600, 2500};
//[350, 420],[420, 500],  [500, 550], [550, 650],[650, 750],/Float_t xbins[] = {350,500, 600, 700, 800, 900, 1000, 1100, 1200, 1400, 1700, 2400}; 
Float_t xbins[] = {350,500, 600, 700, 800, 900, 1000, 1100, 1200, 1400, 1700, 2400}; 
double nbins=11.;
// Float_t xbins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100,120, 140,170, 200,250, 330};
// [[50, 60],[70, 80],[80,100], [100, 120], [120, 140], [140,170],[170, 200],[200, 250] ,[250, 330]]
 
//  Float_t xbins[] = {50, 60, 70, 80, 100,120, 140,170, 200,250, 330};
// double nbins=14.;
// double nbins=10.;
 //     cout << "plotting resolution as a function of qt" << endl;
 //   if (histogramaname.Contains("over")){
//  Float_t xbins[] = {18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 210, 235, 260, 290, 350};
//    double nbins=19.;
//  Float_t xbins[] = {0, 10, 18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 200, 225, 250, 275, 305, 335, 365, 400, 440, 500};
//    double nbins=27.;
 //     cout << "plotting scale " << endl;
 //     this is the new scale binning
  //Float_t xbins[] = {0,10,  18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 200, 225, 250, 275, 305, 335, 365, 400};
//  Float_t xbins[] = {0, 18,24, 30, 38, 46, 52, 60, 68, 76,84,92,100, 115,130,  150, 175, 200, 225, 250, 275, 305, 335, 365, 400, 440, 500};
    //double nbins=25.;
//    double nbins=26.;
//  Float_t xbins[] = {0,20, 30, 40,50,  65,80,95, 110, 130, 160,190, 230, 280, 350, 420, 500};
//    double nbins=16.;
 //   }

   TH1F * hdatadyee= new TH1F("","", nbins,xbins);
   TH1F * hmcdyee= new TH1F("hmcdyee"+histogramaname+Type,"hmcdyee"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyeeup= new TH1F("hmcdyeeup"+histogramaname+Type,"hmcdyeeup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyeedown= new TH1F("hmcdyeedown"+histogramaname+Type,"hmcdyeedown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyeejerup= new TH1F("hmcdyeejerup"+histogramaname+Type,"hmcdyeejerup"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyeejerdown= new TH1F("hmcdyeejerdown"+histogramaname+Type,"hmcdyeejerdown"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdyeeupUncl= new TH1F("hmcdyeeupUncl"+histogramaname+Type,"hmcdyeeupUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hmcdyeedownUncl= new TH1F("hmcdyeedownUncl"+histogramaname+Type,"hmcdyeedownUncl"+histogramaname+Type, nbins, xbins);
   TH1F * hdatadymm= new TH1F("hdatadymm"+histogramaname+Type,"hdatadymm"+histogramaname+Type, nbins,xbins);
   TH1F * hmcdymm= new TH1F("hmcdymm"+histogramaname+Type,"hmcdymm"+histogramaname+Type, nbins,xbins);
   TH1F * hdatadygg= new TH1F("hdatadygg"+histogramaname+Type,"hdatadygg"+histogramaname+Type,  nbins,xbins);
   TH1F * hmcdygg= new TH1F(" "," ", nbins,xbins);
   //TH1F * hmcdygg= new TH1F("hmcdygg"+histogramaname+Type,"hmcdygg"+histogramaname+Type, nbins,xbins);
   
   for (int i = 0, n = nbins; i < n; ++i) {
       hdatadyee->SetBinContent(i+1,datadyee->GetY()[i]);
       hmcdyee->SetBinContent(i+1,mcdyee->GetY()[i]);
       hmcdyeedown->SetBinContent(i+1,mcdyeedown->GetY()[i]);
       hmcdyeeup->SetBinContent(i+1,mcdyeeup->GetY()[i]);     
       hmcdyeejerdown->SetBinContent(i+1,mcdyeejerdown->GetY()[i]);
       hmcdyeejerup->SetBinContent(i+1,mcdyeejerup->GetY()[i]);     
       hmcdyeedownUncl->SetBinContent(i+1,mcdyeedownUncl->GetY()[i]);
       hmcdyeeupUncl->SetBinContent(i+1,mcdyeeupUncl->GetY()[i]);
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
  
    if (Type.Contains("MC")){
        ratiodyee->Divide(hmcdyee);
        ratiodymm->Divide(hmcdymm);
        ratiodygg->Divide(hmcdygg);
        cout << "doing this"<< endl;
    }
    else{
        ratiodyee->Divide(hdatadyee, hmcdyee, 1., 1., "B");
        ratiodymm->Divide(hdatadymm, hmcdymm, 1., 1., "B");
        ratiodygg->Divide(hdatadygg, hmcdygg, 1., 1., "B");
    }
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.014);
    pad1->SetLeftMargin(0.1);
    pad1->SetTicks();
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

if (histogramaname=="met_uPara_over_met_sumEt-zll_pt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_uPara_over_qt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_uParaRaw_over_qt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_uPara_over_zll_pt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_uParaPuppi_over_zll_pt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_uParaRaw_over_zll_pt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_uParaPuppi_over_zll_pt") ylabel="-<u_{||}> / <q_{T}>";
if (histogramaname=="met_metx_over_nVert") ylabel="   <p_{X}^{miss}> ";
if (histogramaname=="met_mety_over_nVert") ylabel="   <p_{Y}^{miss}> ";
if (histogramaname=="met_uPara__vs_zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerp__vs_zll_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_met_sumEt-zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerp_vs_met_sumEt-zll_pt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uParaRaw_vs_zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerpRaw_vs_zll_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uParaRaw_vs_met_sumEt-zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerpRaw_vs_met_sumEt-zll_pt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uParaRaw_vs_nVert") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerpRaw_vs_nVert")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara__vs_met_sumEPara")ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerp__vs_met_sumEPerp")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPara__vs"))ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname.Contains("uPara_vs"))ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname.Contains("uParaPuppi__vs"))ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname.Contains("uPerp__vs"))ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPerp_vs"))ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPerpPuppi__vs"))ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_metx_vs_nVert")ylabel="#sigma ( E^{miss} X) [GeV]";
if (histogramaname=="met_mety_vs_nVert")ylabel="#sigma ( E^{miss} Y) [GeV]";
if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.45,1.15);  
//if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.9,1.1);  

if (Type == "PF_Raw" or Type == "Raw_PFvsPuppiMM" or Type == "Raw_PFvsPuppiEE" or Type == "Raw_PFvsPuppiGJets" or Type == "PF_Raw_0jet" or Type == "PF_Raw_1jet" ){
    hdatadyee->GetYaxis()->SetRangeUser(0.5,1.1);
}
if (Type == "tk"){
    hdatadyee->GetYaxis()->SetRangeUser(0.,1.2);
}
if (Type == "tkRes"){
    hdatadyee->GetYaxis()->SetRangeUser(0.,60.);
}
if (Type == "PF_Type1_1jet" or Type == "PF_Type1_0jet"){
    hdatadyee->GetYaxis()->SetRangeUser(0.8,1.2);
}                                                   

    hdatadymm->Draw("e, X0 same");
    hdatadyee->Draw("AXIS");
    if (Type == "PF_Type1qTAlpha" or Type == "PF_Type1nVertAlpha" or Type == "PF_Type1SumEtAlpha" or Type == "PF_Type1nVertZeta" or Type == "PF_Type1nVertRaw" or Type == "PF_Type1nVertbkgSub" or Type == "PF_Type1nVertnJet50" or Type ==  "PF_Type1nVertnJet60" or Type == "PF_Type1nVertnJet70" or Type == "PF_Type1nVertZetanJet50" or Type == "PF_Type1IdTest" or Type == "PF_Type1qTNoSC" or Type == "PF_Type1nVertNoSC" or  Type == "PF_Type1SumEtNoSC" or Type == "PF_Type1qTMean" or Type == "PF_Type1NewBinning" or Type == "PF_Type1nVertRMS" or Type == "PF_Type1qTRMS" or Type == "PF_Type1SumEtRMS" or Type == "PF_Type1" or Type == "PF_Type1Mean" or Type == "PF_Type1_Scale_sumEt" or Type == "PF_Type1_Scale_nVert" or Type == "PF_Type1MC" or Type == "PF_Raw" or Type=="gammaqt" or Type == "PF_Raw_Res" or Type == "PF_Type1qT" or Type == "PF_Type1nVert"   or Type == "PF_Type1SumEt" or Type == "PF_Type1qTMC" or Type == "PF_Type1nVertMC"   or Type == "PF_Type1SumEtMC"){
        hdatadygg->Draw("e, X0 same");
    }
    hdatadyee->Draw("e, X0 same");
    hdatadymm->Draw("e, X0 same");
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hdatadyee->GetXaxis ()->GetXmin (), 1, hdatadyee->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw();
    }                                                                                                              
    
    TLegend* leg(0);
    leg = new TLegend(0.14,0.04,0.79,0.315);
    if (Type == "Puppi_Type1Mean" or Type == "Puppi_Type1" or Type == "Type1_PFvsPuppiMM" or Type == "Type1_PFvsPuppiEE"){
    leg = new TLegend(0.14,0.04,0.79,0.28);
    }
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.051);
    leg->SetFillColor(0);
    if (Type == "PF_Type1qTNoSC" or Type == "PF_Type1Mean" or Type == "PF_Type1" or  Type == "PF_Type1qTRMS" or  Type == "PF_Type1Mean" or Type == "PF_Type1_Scale_sumEt" or Type == "PF_Type1_Scale_nVert" or Type == "Puppi_Type1FIT" or Type == "PF_Type1MC" or Type == "PF_Type1qT" or Type == "PF_Type1qTMC"   or Type == "PF_Type1qTlowZpt"){
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "PF p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadygg, "PF p_{T}^{miss} #gamma+jets", "lp");
    }                                                                                     
    if (Type == "PF_Type1SumEtNoSC" or Type == "PF_Type1nVert"  or Type == "PF_Type1SumEt"  or Type == "PF_Type1nVertRMS" or Type == "PF_Type1nVertMC"   or Type == "PF_Type1SumEtMC" or Type == "PF_Type1SumEtRMS"){
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "PF p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadygg, "PF p_{T}^{miss} #gamma+jets", "lp");
    }
    if ( Type == "Puppi_Type1MC" or Type == "Puppi_Type1" or  Type == "PF_Type1qTCentral"  or Type == "PF_Type1nVertCentral" or Type == "PF_Type1SumEtCentral"  or Type == "Puppi_Type1_Scale_nVert"  or Type == "Puppi_Type1RMS" or Type == "Puppi_Type1Mean" or Type == "Puppi_Type1Mean"){
    leg->AddEntry(hdatadymm, "PUPPI p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "PUPPI p_{T}^{miss} Z#rightarrow ee", "lp");
    }                                                                                     
    if (Type == "PF_Type1nVertRaw" or Type == "PF_Raw"){
    leg->AddEntry(hdatadymm, "Raw PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "Raw PF p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadygg, "Raw PF p_{T}^{miss} #gamma, q_{T} > 50 GeV", "lp");
    }                                                                                   
    if (Type == "PF_Type1nVertAlpha" or Type == "PF_Type1qTAlpha" or Type == "PF_Type1SumEtAlpha"){                                                        
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss} Z#rightarrow #mu#mu, #alpha < 0.3", "lp");
    leg->AddEntry(hdatadyee, "PF p_{T}^{miss} Z#rightarrow ee, #alpha < 0.3", "lp");
    leg->AddEntry(hdatadygg, "PF p_{T}^{miss} #gamma, q_{T} > 50 GeV, #alpha < 0.3", "lp");
    }                                                                                  
    if (Type == "PF_Type1IdTest"){
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, Tight  Id", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee, Medium Id", "lp");
    leg->AddEntry(hdatadygg, "Z#rightarrow ee, Loose  Id", "lp");
    }                                                                                   
    if (Type == "PF_Type1JetEF"){
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, both electron emEM > 0.90", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee, both electron emEM < 0.90", "lp");
    }                                                                      
    if (Type == "PF_Type1SumEPara" or Type == "PF_Type1SumEPerp"){
    leg->AddEntry(hdatadyee, "Z#rightarrow ee", "lp");
    }

    if (Type == "PF_Type1MatchEle"){                                                                 
    leg->AddEntry(hdatadymm, "PF,  not gen matched electron", "lp");
    leg->AddEntry(hdatadyee, "PF,  gen matched electrons ", "lp");
    //leg->AddEntry(hdatadygg, "p_{T}^{miss} no mathcing", "lp");
    }                                                                                   
    if (Type == "PF_Type1nVertNoJetEta"){                                           
    leg->AddEntry(hdatadymm, "#gamma + jets", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow #mu#mu ", "lp");
    //leg->AddEntry(hdatadygg, "p_{T}^{miss} no mathcing", "lp");
    }                                                                          
    if (Type == "PF_Type1nVertJetEta"){         
    leg->AddEntry(hdatadymm, "#gamma + jets, jet |#eta| < 1.3", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow #mu#mu, jet |#eta| < 1.3", "lp");
    }                                                                          
    if (Type == "PF_Type1nVertRunBCD"){                                           
    leg->AddEntry(hdatadymm, "#gamma + jets, Run BCD", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow #mu#mu, Run BCD", "lp");
    }                                                                          
    if (Type == "PF_Type1FitVsMean"){                                        
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss}, Z#rightarrow ee, mean of fit", "lp");
    leg->AddEntry(hdatadyee, "PF p_{T}^{miss}, Z#rightarrow ee, mean", "lp");
    }                                                                                           
    if (Type == "PF_Type1qTFitVsRMS" or Type == "PF_Type1nVertFitVsRMS"){                                        
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss}, Z#rightarrow ee, FWHM", "lp");
    leg->AddEntry(hdatadyee, "PF p_{T}^{miss}, Z#rightarrow ee, RMS", "lp");
    }                                                                                          

    if (Type == "PF_Type1nVertRunEF"){                                   
    leg->AddEntry(hdatadymm, "#gamma + jets, Run EF", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow #mu#mu, Run EF", "lp");
    }                                                                   
    if (Type == "PF_Type1nVertRunGH"){                                   
    leg->AddEntry(hdatadymm, "#gamma + jets, Run GH", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow #mu#mu, Run GH", "lp");
    }                                                                   

    if (Type == "PF_Type1nVertZeta"){                                                 
    leg->AddEntry(hdatadygg, "#gamma+jets", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, Z #eta < 1.4", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                                                                                
    if (Type == "PF_Type1nVertnJet50"){                                  
    leg->AddEntry(hdatadygg, "#gamma+jets", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, nJet > 0, jet p_{T} > 50", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                                                                                 
    if (Type == "PF_Type1nVertnJet60"){                                  
    leg->AddEntry(hdatadygg, "#gamma+jets", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, nJet > 0, jet p_{T} > 60", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                                                                                
    if (Type == "PF_Type1nVertnJet70"){                                  
    leg->AddEntry(hdatadygg, "#gamma+jets", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, nJet > 0, jet p_{T} > 70", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                              
    if (Type == "PF_Type1nVertZetanJet50"){                                  
    leg->AddEntry(hdatadygg, "#gamma+jets", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, jet p_{T} > 50, Z #eta < 1.4", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                                                                                
    if (Type == "PF_Type1nVertNoSC"){                                  
    leg->AddEntry(hdatadygg, "#gamma+jets, no scale correction", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, no scale correction", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                                                                                
    if (Type == "PF_Type1nVertbkgSub"){                                  
    leg->AddEntry(hdatadygg, "#gamma+jets, subtracted QCD", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee, subtracted tt", "lp");
    leg->AddEntry(hdatadymm, "Z#rightarrow ee", "lp");
    }                                                                            

if (Type == "PF_Raw_Res"){
    leg->AddEntry(hdatadymm, "Raw PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "Raw PF p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadygg, "Raw PF p_{T}^{miss} #gamma, q_{T} > 50 GeV", "lp");
    }                                                                                   
    if (Type=="gammaqt"){
    leg->AddEntry(hdatadymm, "Z#rightarrow #mu#mu PF p_{T}^{miss}", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee     PF p_{T}^{miss}", "lp");
    leg->AddEntry(hdatadygg, "#gamma+jets PF p_{T}^{miss}, q_{T} > 50 GeV", "lp");
    }                                                                                        
    if (Type == "Type1_PFvsPuppiMM"){
    leg->AddEntry(hdatadyee, "PF       p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadymm, "PUPPI p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    }                                                                                    
    if (Type == "Type1_PFvsPuppiEE"){
    leg->AddEntry(hdatadyee, "PF       p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadymm, "PUPPI p_{T}^{miss} Z#rightarrow ee", "lp");
    }                                                                                    
    if (Type == "PF_Raw_0jet" or Type == "PF_Raw_1jet"){                                                                  
    leg->AddEntry(hdatadymm, "Raw PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadygg, "Raw PF p_{T}^{miss} #gamma, q_{T} > 50 GeV", "lp");
    }                                                                                   
    if (Type == "PF_Type1_0jet"){                                                         
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss} Z#rightarrow #mu#mu, 0 jet", "lp");
    leg->AddEntry(hdatadygg, "PF p_{T}^{miss} Z#rightarrow ee, 0 jet", "lp");
    }                                                                                               
    if (Type == "PF_Type1_1jet"){                                                         
    leg->AddEntry(hdatadymm, "PF p_{T}^{miss} Z#rightarrow #mu#mu, > 0 jet", "lp");
    leg->AddEntry(hdatadygg, "PF p_{T}^{miss} Z#rightarrow ee, > 0 jet", "lp");
    }                                                                                              
    if (Type == "Raw_PFvsPuppiMM"){
    leg->AddEntry(hdatadymm, "Raw PF      p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadygg, "Raw PUPPI   p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    }                                                                                       
    if (Type == "Raw_PFvsPuppiEE"){                                                          
    leg->AddEntry(hdatadymm, "Raw PUPPI p_{T}^{miss} Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadyee, "Raw PF    p_{T}^{miss} Z#rightarrow ee", "lp");
    }                                                                                     
    if (Type == "Puppi_RawvsType1"){
    leg->AddEntry(hdatadymm, "Type1 PUPPI p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadygg, "Raw   PUPPI p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    }                                                                                    
    if (Type == "PF_RawvsType1"){
    leg->AddEntry(hdatadymm, "Type1 PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadygg, "Raw   PF p_{T}^{miss} Z#rightarrow #mu#mu", "lp");
    }                                                                                        
    if (Type == "Type1_PFvsPuppiGJets" or Type == "Raw_PFvsPuppiGJets"){
    leg->AddEntry(hdatadymm, "Puppi p_{T}^{miss}", "lp");
    leg->AddEntry(hdatadygg, "Puppi w/ Fix p_{T}^{miss} ", "lp");               
    } 
    
    leg->Draw();
    c1->Modified();
    c1->cd();
 
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.04);
    pad2->SetBottomMargin(0.3);
    pad2->SetTicks();    
    pad2->Draw();
    pad2->cd();
    
    ratiodyee->GetYaxis()->SetTitle("Data / MC");
    ratiodyee->GetXaxis()->SetTitle(EjeX);
    ratiodyee->GetYaxis()->SetRangeUser(0.6, 1.4);
    if (histogramaname.Contains("Puppi_over")) 
        ratiodyee->GetYaxis()->SetRangeUser(0.85,1.15);
    if (histogramaname.Contains("Para_over")) 
        ratiodyee->GetYaxis()->SetRangeUser(0.85,1.15);
    if (histogramaname.Contains("ParaRaw_over")) 
        ratiodyee->GetYaxis()->SetRangeUser(0.8,1.2);    
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
    ratiodyee->GetYaxis()->SetTitleOffset(1.3);
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
   tottoterree.SetFillColor (kGray+2);
   tottoterree.SetFillStyle (3002);                      


   ratiodyee->Draw("e,X0 "); //"e2"
   ratiodyee->Draw("AXIS"); //"e2"
   //ratiodygg->Draw("e1,X0  same "); //"e2"
   if (Type =="PF_Type1qTNoSC" or Type == "PF_Type1qTAlpha" or Type == "PF_Type1SumEtAlpha" or Type == "PF_Type1nVertAlpha" or Type == "PF_Type1nVertRaw" or Type == "PF_Type1SumEtNoSC" or Type == "PF_Type1nVertNoSC" or Type == "PF_Type1MatchEle" or Type == "PF_Type1Mean" or  Type == "PF_Type1Mean" or Type == "PF_Type1" or Type == "PF_Type1nVertRMS"  or Type == "PF_Type1SumEtRMS" or Type == "PF_Type1qTRMS" or Type == "PF_Type1_Scale_sumEt" or Type == "PF_Type1_Scale_nVert" or Type == "PF_Type1MC" or Type == "PF_Raw"  or Type == "PF_Raw_Res" or Type=="gammaqt" or Type == "PF_Type1nVert" or Type == "PF_Type1SumEt" or Type == "PF_Type1qT" or Type == "PF_Type1qTMC" or Type == "PF_Type1nVertMC"   or Type == "PF_Type1SumEtMC"){
       ratiodygg->Draw("e,X0  same"); //"e2"
   }
   ratiodymm->Draw("e,X0  same"); //"e2"
   doErrorLegend = 0;
   if (doErrorLegend){
       tottoterree.Draw ("2 same");
       jererree.Draw ("2 same");
       erree.Draw ("2 same");
       staterree.Draw("2 same");
   }
   else{
       tottoterree.Draw ("2 same");
   }
   //ratiodygg->Draw("e1,X0  same "); //"e2"
   ratiodyee->Draw("e,X0 same"); //"e2"
   ratiodymm->Draw("e,X0  same"); //"e2"
   TLine *lineR =  new TLine ( ratiodyee->GetXaxis ()->GetXmin (), 1, ratiodyee->GetXaxis ()->GetXmax (), 1);
   lineR->SetLineColor (kBlack ); lineR->SetLineWidth (2); lineR->SetLineStyle (2); lineR->Draw();     
        
    TLegend* legratio(0);
   
    legratio = new TLegend(0.25,0.3,0.85,0.5);
    legratio->SetFillColor(0);
    legratio->SetBorderSize(0);
    if (doErrorLegend){
        legratio-> SetNColumns(2);
        legratio->AddEntry(&tottoterree, "Uncl p_{T}^{miss} + JER+ JES + Stat","f");
        legratio->AddEntry(&jererree, "JER + JES + Stat","f");
        legratio->AddEntry(&erree, "JES + Stat","f");
        legratio->AddEntry(&staterree, "Stat","f");
        legratio->Draw();
    }
    else{
        leg->AddEntry(&tottoterree, "Uncertainty","f");
    }
    //    ratiodygg->Draw("e1,  same "); //"e2"
//    ratiodymm->Draw("e1,  same"); //"e2"
//    ratiodyee->Draw("e1,  same"); //"e2"

    pad1->cd();
    //legratio->Draw();
    TLatex latex;                  
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.08);
    latex.DrawLatex(0.26, 0.80, "#bf{CMS}");
    TLatex latexP;                  
    latexP.SetNDC();
    latexP.SetTextAngle(0);
    latexP.SetTextFont(42);
    latexP.SetTextAlign(31);
    latexP.SetTextSize(0.055);
    if ((Type.Contains("FitVs"))) {
        latexP.DrawLatex(0.405, 0.74, "#it{Supplementary}");
    }
    else{
    //latexP.DrawLatex(0.34, 0.74, "#it{Preliminary}");
    latexP.DrawLatex(0.34, 0.74, "");
    }


//    if ((histogramaname.Contains("sumEt"))  or (histogramaname.Contains("nVert")) ) {
        TLatex latexz;
        latexz.SetNDC();
        latexz.SetTextAngle(0);
        latexz.SetTextFont(42);
        latexz.SetTextAlign(31);
        latexz.SetTextSize(0.051);
        if  (histogramaname.Contains("vs")){ 
            if (Type == "PF_Type1qTMC" or Type == "PF_Type1nVertMC"   or Type == "PF_Type1SumEtMC" ){    
                latexz.DrawLatex(0.77, 0.365, "MC only, Response corrected");
            } 
            else if (Type == "Type1_PFvsPuppiMM" or Type == "Type1_PFvsPuppiEE" or Type == "Puppi_Type1") {
                latexz.DrawLatex(0.63, 0.3, "Response corrected");
            }
            else{
                latexz.DrawLatex(0.63, 0.335, "Response corrected");
            }
        }
        else{ 
            if (Type == "PF_Type1MC" or Type == "Puppi_Type1MC"){   
                latexz.DrawLatex(0.4, 0.28, "MC only");
            }
            else{
                latexz.DrawLatex(0.75, 0.28, "");
            }
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
    c1->Print("~/www/met/paper/comparisons/Aug31/"+histogramaname+Type+".png");
    c1->Print("~/www/met/paper/comparisons/Aug31/"+histogramaname+Type+".root");
    c1->Print("~/www/met/paper/comparisons/Aug31/"+histogramaname+Type+".pdf");  


}


void Comparison(TString histograma, TString histograma2, TString EjeX, TString Type) {                                               


TCanvas *c1 = new TCanvas("c1","example",600,700);
if (histograma.Contains("over")){
if (Type == "Type1_PFvsPuppiMM"){
 TFile file1EE ("../tgraphs/DYMZllScaleNov30.root");                              
 TFile file2EE ("../tgraphs/DataMZllScaleNov30.root");                            
 TFile file1MM ("../tgraphs/DYMZllScaleNov30.root");                              
 TFile file2MM ("../tgraphs/DataMZllScaleNov30.root");                            
 TFile file1GG ("../tgraphs/DYMZllScaleNov10Puppi.root");                         
 TFile file2GG ("../tgraphs/DataMZllScaleNov10Puppi.root");                       
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllScaleNov30.root");              
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllScaleNov30.root");            
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYMZllScaleNov30.root");           
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllScaleNov30.root");         
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllScaleNov30.root");             
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllScaleNov30.root");         
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma);                              
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma);                              
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma);                               
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma);                               
 h1GG = (TGraphErrors*) file1GG.Get(histograma2);                                   
 h2GG = (TGraphErrors*) file2GG.Get(histograma2);                                   
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma);                        
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma);                      
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma);                  
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma);                
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma);                
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma);            
}                                                                                   

if (Type == "Raw_PFvsPuppiMM" ){
 TFile file1EE ("../tgraphs/DYMZllScaleNov30.root");                             
 TFile file2EE ("../tgraphs/DataMZllScaleNov30.root");                           
 TFile file1MM ("../tgraphs/DYMZllScaleNov30.root");                             
 TFile file2MM ("../tgraphs/DataMZllScaleNov30.root");                           
 TFile file1GG ("../tgraphs/DYMZllScaleMar13PuppiRaw.root");                        
 TFile file2GG ("../tgraphs/DataMZllScaleMar13PuppiRaw.root");                      
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllScaleNov30.root");             
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllScaleNov30.root");           
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYMZllScaleNov30.root");          
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllScaleNov30.root");        
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllScaleNov30.root");            
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllScaleNov30.root");        
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma);                           
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma);                           
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma);                            
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma);                            
 h1GG = (TGraphErrors*) file1GG.Get(histograma2);                                
 h2GG = (TGraphErrors*) file2GG.Get(histograma2);                                
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma);                     
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma);                   
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma);               
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma);             
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma);             
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma);         
}                                                                                



if (Type == "Type1_PFvsPuppiEE"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay25Mean.root");    //Nov30                         
 TFile file2EE ("../tgraphs/DataEZllScaleMay25Mean.root");                           
 TFile file1MM ("../tgraphs/DYEZllScaleMay25PuppiMean.root");                             
 TFile file2MM ("../tgraphs/DataEZllScaleMay25PuppiMean.root");                           
 TFile file1GG ("../tgraphs/DYEZllScaleMay25PuppiMean.root");      //Nov22                  
 TFile file2GG ("../tgraphs/DataEZllScaleMay25PuppiMean.root");                      
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay25Mean.root");             
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay25Mean.root");           
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay25Mean.root");          
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay25Mean.root");        
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay25Mean.root");            
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay25Mean.root");        
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma);                             
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma);                             
 h1MM = (TGraphErrors*) file1MM.Get(histograma2);                              
 h2MM = (TGraphErrors*) file2MM.Get(histograma2);                              
 h1GG = (TGraphErrors*) file1GG.Get(histograma2);                                  
 h2GG = (TGraphErrors*) file2GG.Get(histograma2);                                  
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma);                       
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma);                     
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma);                 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma);               
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma);               
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);           
}                                                                                  

if (Type == "Raw_PFvsPuppiEE" ){
 TFile file1EE ("../tgraphs/DYEZllScaleNov07.root");                             
 TFile file2EE ("../tgraphs/DataEZllScaleNov07.root");                           
 TFile file1MM ("../tgraphs/DYEZllScaleMar13PuppiRaw.root");                             
 TFile file2MM ("../tgraphs/DataEZllScaleMar13PuppiRaw.root");                           
 TFile file1GG ("../tgraphs/DYEZllScaleMar13PuppiRaw.root");                        
 TFile file2GG ("../tgraphs/DataEZllScaleMar13PuppiRaw.root");                      
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleNov07.root");             
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleNov07.root");           
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleNov07.root");          
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleNov07.root");        
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleNov07.root");            
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleNov07.root");        
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma);                             
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma);                             
 h1MM = (TGraphErrors*) file1MM.Get(histograma2);                              
 h2MM = (TGraphErrors*) file2MM.Get(histograma2);                              
 h1GG = (TGraphErrors*) file1GG.Get(histograma2);                                  
 h2GG = (TGraphErrors*) file2GG.Get(histograma2);                                  
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma);                       
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma);                     
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma);                 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma);               
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma);               
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);           
}                                                                                  


if (Type == "Type1_PFvsPuppiGJets" or Type == "Raw_PFvsPuppiGJets" ){
 TFile file1EE ("../tgraphs/DatagjetsScaleJuly21PuppiNoFix.root");
 TFile file2EE ("../tgraphs/DatagjetsScaleJuly21PuppiNoFix.root");     
 TFile file1MM ("../tgraphs/DatagjetsScaleJuly21PuppiNoFix.root");    
 TFile file2MM ("../tgraphs/DatagjetsScaleJuly21PuppiNoFix.root");     
 TFile file1GG ("../tgraphs/DatagjetsScaleJuly21PuppiFix.root");                        
 TFile file2GG ("../tgraphs/DatagjetsScaleJuly21PuppiFix.root");                      
 TFile file1EEup     ("../tgraphs/GJets_up_jes_DYgjetsScaleApril17.root"); 
 TFile file1EEdown   ("../tgraphs/GJets_down_jes_DYgjetsScaleApril17.root"); 
 TFile file1EEjerup     ("../tgraphs/GJets_up_jer_DYgjetsScaleApril17.root"); 
 TFile file1EEjerdown   ("../tgraphs/GJets_down_jer_DYgjetsScaleApril17.root");
 TFile file1EEupUncl ("../tgraphs/GJets_up_uncl_DYgjetsScaleApril17.root");
 TFile file1EEdownUncl ("../tgraphs/GJets_down_uncl_DYgjetsScaleApril17.root");
 h1EE = (TGraphErrors*) file1EE.Get(histograma2);                             
 h2EE = (TGraphErrors*) file2EE.Get(histograma2);                             
 h1MM = (TGraphErrors*) file1MM.Get(histograma2);                              
 h2MM = (TGraphErrors*) file2MM.Get(histograma2);                              
 h1GG = (TGraphErrors*) file1GG.Get(histograma2);                                  
 h2GG = (TGraphErrors*) file2GG.Get(histograma2);                                  
 h1EEup = (TGraphErrors*)   file1EEup.Get(histograma);                       
 h1EEdown = (TGraphErrors*) file1EEdown.Get(histograma);                     
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get(histograma);                 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get(histograma);               
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get(histograma);               
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get(histograma);           
}                                                                                            

if (Type == "PF_Raw" ){
 TFile file1EE ("../tgraphs/DYEZllScaleNov07.root");
 TFile file2EE ("../tgraphs/DataEZllScaleNov07.root");
 TFile file1MM ("../tgraphs/DYMZllScaleNov07.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleNov07.root");
 TFile file1GG ("../tgraphs/GJetsgjetsScaleNov07.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleNov07.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleNov07.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleNov07.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleNov07.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleNov07.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleNov07.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleNov07.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                   

if (Type == "PF_Type1"){
 TFile file1EE ("../tgraphs/DYEZllScaleApril30.root");//Nov30
 TFile file2EE ("../tgraphs/DataEZllScaleNov30.root"); 
 TFile file1MM ("../tgraphs/DYMZllScaleApril30.root"); //Feb21
 TFile file2MM ("../tgraphs/DataMZllScaleFeb21.root");//Feb21
 TFile file1GG ("../tgraphs/GJetsgjetsScaleApril30.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleNov30.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril30.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril30.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril30.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril30.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril30.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril30.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                         
}                                                                                             

if (Type == "PF_Type1IdTest"){
 TFile file1EE ("../tgraphs/DataEZllScaleMay02FirstBinFix.root");//Nov30
 TFile file2EE ("../tgraphs/DataEZllScaleMay02FirstBinFix.root"); 
 TFile file1MM ("../tgraphs/DataEZllScaleMay08MediumId.root"); //Feb21
 TFile file2MM ("../tgraphs/DataEZllScaleMay08MediumId.root");//Feb21
 TFile file1GG ("../tgraphs/DataEZllScaleMay08LooseId.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleMay08LooseId.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril30.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril30.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril30.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril30.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril30.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril30.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                   


if (Type == "PF_Type1JetEF"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay09matchDoubleEle.root");//Nov30
 TFile file2EE ("../tgraphs/DYEZllScaleMay09matchDoubleEle.root"); 
 TFile file1MM ("../tgraphs/DYEZllScaleMay09matchNoEle.root"); //Feb21
 TFile file2MM ("../tgraphs/DYEZllScaleMay09matchNoEle.root");//Feb21
 TFile file1GG ("../tgraphs/DYEZllScaleMay09matchNoEle.root"); 
 TFile file2GG ("../tgraphs/DYEZllScaleMay09matchNoEle.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril30.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril30.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril30.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril30.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril30.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril30.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                   








if (Type == "PF_Type1NewBinning"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay02FirstBinFix.root");//Nov30
 TFile file2EE ("../tgraphs/DataEZllScaleMay02FirstBinFix.root"); 
 TFile file1MM ("../tgraphs/DYMZllScaleMay02FirstBinFix.root"); //Feb21
 TFile file2MM ("../tgraphs/DataMZllScaleMay02FirstBinFix.root");//Feb21
 TFile file1GG ("../tgraphs/GJetsgjetsScaleMay02FirstBinFix.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleMay02FirstBinFix.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay02FirstBinFix.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay02FirstBinFix.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                      
}                                                                                             

if (Type == "PF_Type1MC"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay02FirstBinFix.root");//Nov30
 TFile file2EE ("../tgraphs/DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1MM ("../tgraphs/DYMZllScaleMay02FirstBinFix.root"); //Feb21
 TFile file2MM ("../tgraphs/DYMZllScaleMay02FirstBinFix.root");//Feb21
 TFile file1GG ("../tgraphs/GJetsgjetsScaleMay02FirstBinFix.root"); 
 TFile file2GG ("../tgraphs/GJetsgjetsScaleMay02FirstBinFix.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay02FirstBinFix.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay02FirstBinFix.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay02FirstBinFix.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                      
}                                                                                             

if (Type == "PF_Type1Mean"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay25Mean.root");//Nov30
 TFile file2EE ("../tgraphs/DataEZllScaleMay25Mean.root"); 
 TFile file1MM ("../tgraphs/DYMZllScaleMay25Mean.root"); //Feb21
 TFile file2MM ("../tgraphs/DataMZllScaleMay25Mean.root");//Feb21
 TFile file1GG ("../tgraphs/GJetsgjetsScaleMay25Mean.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleMay25Mean.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay25Mean.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay25Mean.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay25Mean.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay25Mean.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay25Mean.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay25Mean.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                        
}                                                                                         

if (Type == "PF_Type1FitVsMean"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay25Mean.root");//Nov30
 TFile file2EE ("../tgraphs/DataEZllScaleMay25Mean.root"); 
 TFile file1MM ("../tgraphs/DYEZllScaleMay02FirstBinFix.root"); //Feb21
 TFile file2MM ("../tgraphs/DataEZllScaleMay02FirstBinFix.root");//Feb21
 TFile file1GG ("../tgraphs/DYEZllScaleMay02FirstBinFix.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleMay02FirstBinFix.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay25Mean.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay25Mean.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay25Mean.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay25Mean.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay25Mean.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay25Mean.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                  
}                                                                                         





if (Type == "PF_Type1_Scale_nVert"){
 TFile file1EE ("../tgraphs/DYEZllScaleApril18nVertZ50.root");
 TFile file2EE ("../tgraphs/DataEZllScaleApril18nVertZ50.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril18nVertZ50.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril18nVertZ50.root");
 TFile file1GG ("../tgraphs/GJetsgjetsScaleApril18nVert.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleApril18nVert.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril18nVertZ50.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril18nVertZ50.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril18nVertZ50.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril18nVertZ50.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril18nVertZ50.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril18nVertZ50.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                                      
if (Type == "PF_Type1_Scale_sumEt"){
 TFile file1EE ("../tgraphs/DYEZllScaleApril18SumEtZ50.root");
 TFile file2EE ("../tgraphs/DataEZllScaleApril18SumEtZ50.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril18SumEtZ50.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril18SumEtZ50.root");
 TFile file1GG ("../tgraphs/GJetsgjetsScaleApril18SumEt.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleApril18SumEt.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril18SumEtZ50.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril18SumEtZ50.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril18SumEtZ50.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril18SumEtZ50.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril18SumEtZ50.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril18SumEtZ50.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                            

if (Type == "Puppi_Type1_Scale_nVert"){
 TFile file1EE ("../tgraphs/DYEZllScaleApril23PuppinVertZ100.root");
 TFile file2EE ("../tgraphs/DataEZllScaleApril23PuppinVertZ100.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril23PuppinVertZ100.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril23PuppinVertZ100.root");
 TFile file1GG ("../tgraphs/DYEZllScaleApril23PuppinVertZ100.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleApril23PuppinVertZ100.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril23PuppinVertZ100.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril23PuppinVertZ100.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril23PuppinVertZ100.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril23PuppinVertZ100.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril23PuppinVertZ100.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril23PuppinVertZ100.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                                  

if (Type == "PF_Type1MatchEle"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay07matchEle0.root");
 TFile file2EE ("../tgraphs/DYEZllScaleMay07matchEle0.root");
 TFile file1MM ("../tgraphs/DYEZllScaleMay07matchNoEle0.root"); 
 TFile file2MM ("../tgraphs/DYEZllScaleMay07matchNoEle0.root");
 TFile file1GG ("../tgraphs/DYEZllScaleMay07matchEle.root"); 
 TFile file2GG ("../tgraphs/DYEZllScaleMay07matchEle.root");
 TFile file1EEup     ("../tgraphs/DYEZllScaleMay07matchEle.root"); 
 TFile file1EEdown   ("../tgraphs/DYEZllScaleMay07matchEle.root"); 
 TFile file1EEjerup     ("../tgraphs/DYEZllScaleMay07matchEle.root"); 
 TFile file1EEjerdown   ("../tgraphs/DYEZllScaleMay07matchEle.root"); 
 TFile file1EEupUncl ("../tgraphs/DYEZllScaleMay07matchEle.root");
 TFile file1EEdownUncl ("../tgraphs/DYEZllScaleMay07matchEle.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                            


//if (Type == "PF_Type1MC"){
// TFile file1EE ("../tgraphs/DYEZllScaleApril30.root");
// TFile file2EE ("../tgraphs/DYEZllScaleApril30.root");
// TFile file1MM ("../tgraphs/DYMZllScaleApril30.root"); 
// TFile file2MM ("../tgraphs/DYMZllScaleApril30.root");
// TFile file1GG ("../tgraphs/GJetsgjetsScaleApril30.root"); 
// TFile file2GG ("../tgraphs/GJetsgjetsScaleApril30.root");
// TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril30.root"); 
// TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril30.root"); 
// TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril30.root"); 
// TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril30.root"); 
// TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril30.root");
// TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril30.root");
// h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
// h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
// h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
// h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
// h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
// h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
// h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
// h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
// h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
// h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
// h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
// h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
//}                                                                            

if (Type == "Puppi_Type1NewBinning"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay03PuppiFirstBinFix.root");//April30
 TFile file2EE ("../tgraphs/DataEZllScaleMay02PuppiFirstBinFix.root");//April17
 TFile file1MM ("../tgraphs/DYMZllScaleMay02PuppiFirstBinFix.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleMay02PuppiFirstBinFix.root");
 TFile file1GG ("../tgraphs/DYEZllScaleMay02PuppiFirstBinFix.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleMay02PuppiFirstBinFix.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay03PuppiFirstBinFix.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay03PuppiFirstBinFix.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 //h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 //h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma);//Nov22 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                                     

if (Type == "Puppi_Type1MC"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay03PuppiFirstBinFix.root");//April30
 TFile file2EE ("../tgraphs/DYEZllScaleMay02PuppiFirstBinFix.root");//April17
 TFile file1MM ("../tgraphs/DYMZllScaleMay02PuppiFirstBinFix.root"); 
 TFile file2MM ("../tgraphs/DYMZllScaleMay02PuppiFirstBinFix.root");
 TFile file1GG ("../tgraphs/DYEZllScaleMay02PuppiFirstBinFix.root"); 
 TFile file2GG ("../tgraphs/DYEZllScaleMay02PuppiFirstBinFix.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay03PuppiFirstBinFix.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay03PuppiFirstBinFix.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay03PuppiFirstBinFix.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 //h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 //h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma);//Nov22 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                                    





if (Type == "Puppi_Type1"){
 TFile file1EE ("../tgraphs/DYEZllScaleApril30Puppi.root");//April30
 TFile file2EE ("../tgraphs/DataEZllScaleApril17Puppi.root");//April17
 TFile file1MM ("../tgraphs/DYMZllScaleApril30Puppi.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril17Puppi.root");
 TFile file1GG ("../tgraphs/DYEZllScaleApril30Puppi.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleApril17Puppi.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleApril30Puppi.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleApril30Puppi.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleApril30Puppi.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleApril30Puppi.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleApril30Puppi.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleApril30Puppi.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 //h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 //h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma);//Nov22 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                    

if (Type == "Puppi_Type1Mean"){
 TFile file1EE ("../tgraphs/DYEZllScaleMay25PuppiMean.root");//April30
 TFile file2EE ("../tgraphs/DataEZllScaleMay25PuppiMean.root");//April17
 TFile file1MM ("../tgraphs/DYMZllScaleMay25PuppiMean.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleMay25PuppiMean.root");
 TFile file1GG ("../tgraphs/DYEZllScaleMay25PuppiMean.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleMay25PuppiMean.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleMay25PuppiMean.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleMay25PuppiMean.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleMay25PuppiMean.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleMay25PuppiMean.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleMay25PuppiMean.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleMay25PuppiMean.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 //h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 //h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma);//Nov22 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                    


if (Type == "PF_Raw_0jet" or Type == "PF_Type1_0jet"){
 TFile file1EE ("../tgraphs/DYEZllScale0jetNov07.root");
 TFile file2EE ("../tgraphs/DataEZllScale0jetNov07.root");
 TFile file1MM ("../tgraphs/DYEZllScale0jetNov07.root"); 
 TFile file2MM ("../tgraphs/DataMZllScale0jetNov07.root");
 TFile file1GG ("../tgraphs/DYEZllScale0jetNov07.root"); 
 TFile file2GG ("../tgraphs/DataEZllScale0jetNov07.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScale0jetNov07.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScale0jetNov07.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScale0jetNov07.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jes_DYEZllScale0jetNov07.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScale0jetNov07.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScale0jetNov07.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                  

if (Type == "PF_Raw_1jet" or Type == "PF_Type1_1jet"){
 TFile file1EE ("../tgraphs/DYEZllScaleNov071jet.root");
 TFile file2EE ("../tgraphs/DataEZllScaleNov071jet.root");
 TFile file1MM ("../tgraphs/DYMZllScaleNov071jet.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleNov071jet.root");
 TFile file1GG ("../tgraphs/DYEZllScaleNov071jet.root"); 
 TFile file2GG ("../tgraphs/DataEZllScaleNov071jet.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllScaleNov071jet.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllScaleNov071jet.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllScaleNov071jet.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllScaleNov071jet.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllScaleNov071jet.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllScaleNov071jet.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 
}                                                                                  

if (Type == "PF_RawvsType1"){
 TFile file1EE ("../tgraphs/DYMZllScaleApril17.root");
 TFile file2EE ("../tgraphs/DataMZllScaleApril17.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril17.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril17.root");
 TFile file1GG ("../tgraphs/DYMZllScaleApril17.root"); 
 TFile file2GG ("../tgraphs/DataMZllScaleApril17.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllScaleApril17.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllScaleApril17.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYMZllScaleApril17.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllScaleApril17.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllScaleApril17.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllScaleApril17.root");
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("M_"+histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get("M_"+histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma); 
}                                                                                 

if (Type == "Puppi_RawvsType1"){
 TFile file1EE ("../tgraphs/DYMZllScaleApril17Puppi.root");
 TFile file2EE ("../tgraphs/DataMZllScaleApril17Puppi.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril17Puppi.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril17Puppi.root");
 TFile file1GG ("../tgraphs/DYMZllScaleApril17Puppi.root"); 
 TFile file2GG ("../tgraphs/DataMZllScaleApril17Puppi.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_down_jer_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllScaleApril17Puppi.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllScaleApril17Puppi.root");
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("M_"+histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get("M_"+histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma); 
}                                                                                       

if (Type == "Type1Puppi_comp" ){
 TFile file1EE ("../tgraphs/DYEZllScaleMay22Puppi.root");
 TFile file2EE ("../tgraphs/DataEZllScaleMay22Puppi.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril17Puppi.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril17Puppi.root");
 TFile file1GG ("../tgraphs/GJetsgjetsScaleMay23.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleMay23.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_down_jer_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllScaleApril17Puppi.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllScaleApril17Puppi.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma); 
}                                                                                      

if (Type == "Type1_PFvsPuppi_" ){
 TFile file1EE ("../tgraphs/DYEZllScaleMay29.root");
 TFile file2EE ("../tgraphs/DataEZllScaleMay29.root");
 TFile file1MM ("../tgraphs/DYMZllScaleApril17.root"); 
 TFile file2MM ("../tgraphs/DataMZllScaleApril17.root");
 TFile file1GG ("../tgraphs/DYMZllScaleApril17.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleMay23.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_down_jer_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllScaleApril17Puppi.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllScaleApril17Puppi.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllScaleApril17Puppi.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma); 
}                                                                                      





if (Type == "gammaType1"){
 TFile file1EE ("../tgraphs/GJetsgjetsScaleApril17.root");
 TFile file2EE ("../tgraphs/DatagjetsScaleApril17.root");
 TFile file1MM ("../tgraphs/GJetsgjetsScaleApril17.root"); 
 TFile file2MM ("../tgraphs/DatagjetsScaleApril17.root");  
 TFile file1GG ("../tgraphs/GJetsgjetsScaleFeb15_ReReco.root"); 
 TFile file2GG ("../tgraphs/DatagjetsScaleFeb15_ReReco.root");
 TFile file1EEup     ("../tgraphs/GJets_up_jes_DYgjetsScaleApril17.root"); 
 TFile file1EEdown   ("../tgraphs/GJets_down_jes_DYgjetsScaleApril17.root"); 
 TFile file1EEjerup     ("../tgraphs/GJets_up_jer_DYgjetsScaleApril17.root"); 
 TFile file1EEjerdown   ("../tgraphs/GJets_down_jer_DYgjetsScaleApril17.root"); 
 TFile file1EEupUncl ("../tgraphs/GJets_up_uncl_DYgjetsScaleApril17.root");
 TFile file1EEdownUncl ("../tgraphs/GJets_down_uncl_DYgjetsScaleApril17.root");
 h1EE = (TGraphErrors*) file1EE.Get( histograma); 
 h2EE = (TGraphErrors*) file2EE.Get( histograma); 
 h1MM = (TGraphErrors*) file1MM.Get(histograma); 
 h2MM = (TGraphErrors*) file2MM.Get(histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get(histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get(histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get(histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get(histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get(histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get(histograma); 
}                                                                                         
}
else{

if (Type == "PF_Type1qT"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25qTRMS.root");   //Nov16 is OG one//done w/o the < 1 error requirement
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25qTRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25qTRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay25qTRMS.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25qTRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay25qTRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25qTRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25qTRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+ histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                 
}                                                                                                                                                                    

if (Type == "PF_Type1qTFitVsRMS"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25qTRMS.root");   //Nov16 is OG one//done w/o the < 1 error requirement
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25qTRMS.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay01qT.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay01qT.root");
 TFile file1GG ("../tgraphs/DYEZllResolutionMay01qT.root"); 
 TFile file2GG ("../tgraphs/DataEZllResolutionMay01qT.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25qTRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25qTRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                 
}                                                                                                                                  








if (Type == "PF_Type1qTAlpha"){
 TFile file1EE ("../tgraphs/DYEZllResolutionJune06qTRMS_alpha.root");   //Nov16 is OG one//done w/o the < 1 error requirement
 TFile file2EE ("../tgraphs/DataEZllResolutionJune06qTRMS_alpha.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionJune06qTRMS_alpha.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionJune06qTRMS_alpha.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25qTRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay25qTRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionJune06qTRMS_alpha.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionJune06qTRMS_alpha.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionJune06qTRMS_alpha.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionJune06qTRMS_alpha.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionJune06qTRMS_alpha.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionJune06qTRMS_alpha.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+ histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                 
}                                                                                                                                        






if (Type == "PF_Type1qTMC"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25qTRMS.root");   //Nov16 is OG one//done w/o the < 1 error requirement
 TFile file2EE ("../tgraphs/DYEZllResolutionMay25qTRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25qTRMS.root"); 
 TFile file2MM ("../tgraphs/DYMZllResolutionMay25qTRMS.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25qTRMS.root"); 
 TFile file2GG ("../tgraphs/GJetsgjetsResolutionMay25qTRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25qTRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25qTRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25qTRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+ histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                 
}                                                                                                                                     


if (Type == "PF_Type1qTRMS"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay01qTMean.root");   //Nov16 is OG one//done w/o the < 1 error requirement
 TFile file2EE ("../tgraphs/DataEZllResolutionMay01qTMean.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay01qTMean.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay01qTMean.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay01qTMean.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay01qTMean.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay01qTMean.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay01qTMean.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay01qTMean.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay01qTMean.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay01qTMean.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay01qTMean.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+ histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                 
}                                                                                                                                      


if (Type == "PF_Type1qTNoSC"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay04qTNoSC.root");   //Nov16 is OG one
 TFile file2EE ("../tgraphs/DataEZllResolutionMay04qTNoSC.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay04qTNoSC.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay04qTNoSC.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay04qTNoSC.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay04qTNoSC.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay04qTNoSC.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay04qTNoSC.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay04qTNoSC.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay04qTNoSC.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay04qTNoSC.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay04qTNoSC.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);              
}                                                                                     

if (Type == "PF_Type1SumEt"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25sumEtRMS.root");  //Mar27_qT50
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25sumEtRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25sumEtRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay25sumEtRMS.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25sumEtRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay25sumEtRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25sumEtRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25sumEtRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                                       


if (Type == "PF_Type1SumEtAlpha"){
 TFile file1EE ("../tgraphs/DYEZllResolutionJune06sumEtRMS_alpha.root");  //Mar27_qT50
 TFile file2EE ("../tgraphs/DataEZllResolutionJune06sumEtRMS_alpha.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionJune06sumEtRMS_alpha.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionJune06sumEtRMS_alpha.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionJune06sumEtRMS_alpha.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionJune06sumEtRMS_alpha.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionJune06sumEtRMS_alpha.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionJune06sumEtRMS_alpha.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionJune06sumEtRMS_alpha.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionJune06sumEtRMS_alpha.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionJune06sumEtRMS_alpha.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionJune06sumEtRMS_alpha.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                                







if (Type == "PF_Type1SumEPara"){
 TFile file1EE ("../tgraphs/DataEZllResolutionMay25sumEParaRMS.root");  //Mar27_qT50
 TFile file2EE ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root");
 TFile file1MM ("../tgraphs/DataEZllResolutionMay25sumEParaRMS.root"); 
 TFile file2MM ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root");
 TFile file1GG ("../tgraphs/DataEZllResolutionMay25sumEParaRMS.root"); 
 TFile file2GG ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root");
 TFile file1EEup     ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DYEZllResolutionMay25sumEParaRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                                      

if (Type == "PF_Type1SumEPerp"){
 TFile file1EE ("../tgraphs/DataEZllResolutionMay25sumEPerpRMS.root");  //Mar27_qT50
 TFile file2EE ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root");
 TFile file1MM ("../tgraphs/DataEZllResolutionMay25sumEPerpRMS.root"); 
 TFile file2MM ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root");
 TFile file1GG ("../tgraphs/DataEZllResolutionMay25sumEPerpRMS.root"); 
 TFile file2GG ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root");
 TFile file1EEup     ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DYEZllResolutionMay25sumEPerpRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                              


if (Type == "PF_Type1SumEtMC"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25sumEtRMS.root");  //Mar27_qT50
 TFile file2EE ("../tgraphs/DYEZllResolutionMay25sumEtRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25sumEtRMS.root"); 
 TFile file2MM ("../tgraphs/DYMZllResolutionMay25sumEtRMS.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25sumEtRMS.root"); 
 TFile file2GG ("../tgraphs/GJetsgjetsResolutionMay25sumEtRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25sumEtRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25sumEtRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25sumEtRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                       



if (Type == "PF_Type1SumEtNoSC"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay04sumEtNoSC.root");  //Mar27_qT50
 TFile file2EE ("../tgraphs/DataEZllResolutionMay04sumEtNoSC.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay04sumEtNoSC.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay04sumEtNoSC.root");
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay04sumEtNoSC.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay04sumEtNoSC.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay04sumEtNoSC.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay04sumEtNoSC.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay04sumEtNoSC.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay04sumEtNoSC.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay04sumEtNoSC.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay04sumEtNoSC.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);            
}                                                                                       



if (Type == "PF_Type1nVert"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25nVertRMS.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25nVertRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay25nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay25nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                               
}                                                                                           


if (Type == "PF_Type1nVertFitVsRMS"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25nVertRMS.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25nVertRMS.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay03nVert.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay01nVert.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/DYEZllResolutionMay03nVert.root"); 
 TFile file2GG ("../tgraphs/DataEZllResolutionMay01nVert.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           





if (Type == "PF_Type1nVertRunBCD"){
 TFile file1EE  ("../tgraphs/DataEZllResolutionJune08nVert_runBCD.root"); //Mar20_nVert
 TFile file2EE  ("../tgraphs/DataEZllResolutionJune08nVert_runBCD.root");
 TFile file1MM ("../tgraphs/DatagjetsResolutionJune08nVert_runBCD.root"); 
 TFile file2MM ("../tgraphs/DatagjetsResolutionJune08nVert_runBCD.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/DatagjetsResolutionJune08nVert_runBCD.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionJune08nVert_runBCD.root");
 TFile file1EEup     ("../tgraphs/DataEZllResolutionJune08nVert_runB.root"); 
 TFile file1EEdown   ("../tgraphs/DataEZllResolutionJune08nVert_runB.root"); 
 TFile file1EEjerup     ("../tgraphs/DataEZllResolutionJune08nVert_runB.root"); 
 TFile file1EEjerdown   ("../tgraphs/DataEZllResolutionJune08nVert_runB.root"); 
 TFile file1EEupUncl ("../tgraphs/DataEZllResolutionJune08nVert_runB.root");
 TFile file1EEdownUncl ("../tgraphs/DataEZllResolutionJune08nVert_runB.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get(histograma2); 
 h2MM = (TGraphErrors*) file2MM.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           




if (Type == "PF_Type1nVertJetEta"){
 TFile file1EE ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root");
 TFile file1MM ("../tgraphs/GJetsgjetsResolutionJune06nVert_JetEta.root"); 
 TFile file2MM ("../tgraphs/GJetsgjetsResolutionJune06nVert_JetEta.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionJune06nVert_JetEta.root"); 
 TFile file2GG ("../tgraphs/GJetsgjetsResolutionJune06nVert_JetEta.root");
 TFile file1EEup     ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root"); 
 TFile file1EEdown   ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root"); 
 TFile file1EEjerup     ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root"); 
 TFile file1EEjerdown   ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root"); 
 TFile file1EEupUncl ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root");
 TFile file1EEdownUncl ("../tgraphs/DYMZllResolutionJune07nVertJetEta.root");
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get(histograma2); 
 h2MM = (TGraphErrors*) file2MM.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma);                    
}                                                                                           



if (Type == "PF_Type1nVertNoJetEta"){
 TFile file1EE ("../tgraphs/DYMZllResolutionJune07nVert.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DYMZllResolutionJune07nVert.root");
 TFile file1MM ("../tgraphs/GJetsgjetsResolutionJune06nVert.root"); 
 TFile file2MM ("../tgraphs/GJetsgjetsResolutionJune06nVert.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionJune06nVert.root"); 
 TFile file2GG ("../tgraphs/GJetsgjetsResolutionJune06nVert.root");
 TFile file1EEup     ("../tgraphs/DYMZllResolutionJune07nVert.root"); 
 TFile file1EEdown   ("../tgraphs/DYMZllResolutionJune07nVert.root"); 
 TFile file1EEjerup     ("../tgraphs/DYMZllResolutionJune07nVert.root"); 
 TFile file1EEjerdown   ("../tgraphs/DYMZllResolutionJune07nVert.root"); 
 TFile file1EEupUncl ("../tgraphs/DYMZllResolutionJune07nVert.root");
 TFile file1EEdownUncl ("../tgraphs/DYMZllResolutionJune07nVert.root");
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get(histograma2); 
 h2MM = (TGraphErrors*) file2MM.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma);                    
}                                                                                           







if (Type == "PF_Type1nVertAlpha"){
 TFile file1EE ("../tgraphs/DYEZllResolutionJune06nVertRMS_alpha.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionJune06nVertRMS_alpha.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionJune06nVertRMS_alpha.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionJune06nVertRMS_alpha.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionJune06nVertRMS_alpha.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionJune06nVertRMS_alpha.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionJune06nVertRMS_alpha.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionJune06nVertRMS_alpha.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionJune06nVertRMS_alpha.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionJune06nVertRMS_alpha.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionJune06nVertRMS_alpha.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionJune06nVertRMS_alpha.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           






if (Type == "PF_Type1nVertRaw"){
 TFile file1EE ("../tgraphs/DYEZllResolutionJune05nVertRMS.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionJune05nVertRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionJune05nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionJune05nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionJune05nVert.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionJune05nVert.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionJune05nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionJune05nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionJune05nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionJune05nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_jes_DYEZllResolutionJune05nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_up_jes_DYEZllResolutionJune05nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           




if (Type == "PF_Type1nVertZeta"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS_Zeta.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_Zeta.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           

if (Type == "PF_Type1nVertnJet50"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS_nJetPt50.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_nJetPt50.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           


if (Type == "PF_Type1nVertnJet60"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS_nJetPt60.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_nJetPt60.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           

if (Type == "PF_Type1nVertnJet70"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS_nJetPt70.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_nJetPt70.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           


if (Type == "PF_Type1nVertZetanJet50"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS_nJetPt50_Zeta.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_nJetPt50_Zeta.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                                        

if (Type == "PF_Type1nVertNoSC"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS_NoSC.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_NoSC.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS_NoSC.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS_NoSC.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                                       


if (Type == "PF_Type1nVertbkgSub"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DataEZllResolutionMay30nVertRMS_bkgSub.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay30nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay30nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay30nVertRMS.root"); 
 TFile file2GG ("../tgraphs/DatagjetsResolutionMay30nVertRMS_bkgSub.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("E_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("E_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                                





if (Type == "PF_Type1nVertMC"){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25nVertRMS.root"); //Mar20_nVert
 TFile file2EE ("../tgraphs/DYEZllResolutionMay25nVertRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25nVertRMS.root"); 
 TFile file2MM ("../tgraphs/DYMZllResolutionMay25nVertRMS.root");//Mar09_qT100
 TFile file1GG ("../tgraphs/GJetsgjetsResolutionMay25nVertRMS.root"); 
 TFile file2GG ("../tgraphs/GJetsgjetsResolutionMay25nVertRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);                    
}                                                                                           



if (Type == "Puppi_Type1" ){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25PuppiRMS.root"); //Nov22
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25PuppiRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25PuppiRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay25PuppiRMS.root");
 TFile file1GG ("../tgraphs/DYEZllResolutionMay25PuppiRMS.root"); 
 TFile file2GG ("../tgraphs/DataEZllResolutionMay25PuppiRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25PuppiRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25PuppiRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25PuppiRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25PuppiRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25PuppiRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25PuppiRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get("E_"+ histograma); 
 h2GG = (TGraphErrors*) file2GG.Get("E_"+ histograma); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);               
}                                                                                      

if (Type == "Type1_PFvsPuppiMM" or Type == "Raw_PFvsPuppiMM" ){
 TFile file1EE ("../tgraphs/DYMZllResolutionMay025nVertRMS.root");
 TFile file2EE ("../tgraphs/DataMZllResolutionMay025nVertRMS.root");
 TFile file1MM ("../tgraphs/DYMZllResolutionMay25PuppiRMS.root"); 
 TFile file2MM ("../tgraphs/DataMZllResolutionMay25PuppiRMS.root");
 TFile file1GG ("../tgraphs/DYMZllResolutionMay25PuppiRMS.root");
 TFile file2GG ("../tgraphs/DataMZllResolutionMay25PuppiRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYMZllResolutionMay025nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYMZllResolutionMay025nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYMZllResolutionMay025nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYMZllResolutionMay025nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYMZllResolutionMay025nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYMZllResolutionMay025nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("M_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("M_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get(histograma2); 
 h2MM = (TGraphErrors*) file2MM.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("M_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("M_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("M_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("M_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("M_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("M_"+ histograma);            
}                                                                                           

if (Type == "Type1_PFvsPuppiEE" or Type == "Raw_PFvsPuppiEE" ){
 TFile file1EE ("../tgraphs/DYEZllResolutionMay25nVertRMS.root"); //Nov16
 TFile file2EE ("../tgraphs/DataEZllResolutionMay25nVertRMS.root");
 TFile file1MM ("../tgraphs/DYEZllResolutionMay25PuppiRMS.root"); 
 TFile file2MM ("../tgraphs/DataEZllResolutionMay25PuppiRMS.root");
 TFile file1GG ("../tgraphs/DYEZllResolutionMay25PuppiRMS.root");//Nov22
 TFile file2GG ("../tgraphs/DataEZllResolutionMay25PuppiRMS.root");
 TFile file1EEup     ("../tgraphs/DY_up_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEdown   ("../tgraphs/DY_down_jes_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerup     ("../tgraphs/DY_up_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEjerdown   ("../tgraphs/DY_down_jer_DYEZllResolutionMay25nVertRMS.root"); 
 TFile file1EEupUncl ("../tgraphs/DY_up_uncl_DYEZllResolutionMay25nVertRMS.root");
 TFile file1EEdownUncl ("../tgraphs/DY_down_uncl_DYEZllResolutionMay25nVertRMS.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get(histograma2); 
 h2MM = (TGraphErrors*) file2MM.Get(histograma2); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEjerup = (TGraphErrors*)   file1EEjerup.Get("E_"+ histograma); 
 h1EEjerdown = (TGraphErrors*) file1EEjerdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma);           
}                                                                                  

}

    cout << " h1EE " << h1EE << endl ;
    cout << " h2EE " << h2EE << endl ;            
    cout << " h1MM " << h1MM << endl ;               
    cout << " h2MM " << h2MM << endl ;               
    cout << " h1GG " << h1GG << endl ;               
    cout << " h2GG " << h2GG << endl ;               
    cout << " h1EEup " << h1EEup << endl ;
    cout << " h1EEdown " << h1EEdown << endl ;
    cout << " h1EEunclUp " << h1EEupUncl << endl ;
    cout << " h1EEunclDown " << h1EEdownUncl << endl ;

   gStyle->SetOptStat(0);
//   gStyle->SetOptFit(1111);
   gStyle->SetErrorX(0.5);     

   PlotWithRatioEEMuMuGJetsComparison(c1, h1EE, h2EE, h1MM, h2MM, h1GG,h2GG, h1EEup, h1EEdown,  h1EEjerup, h1EEjerdown, h1EEupUncl, h1EEdownUncl, EjeX, histograma, Type);
}


TGraphAsymmErrors JESerror(TH1F *Background, TH1F *Backgrounddown, TH1F *Backgroundup, TH1F * hdata){
    
    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");
	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	
    const Int_t nvar = 28;
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
	const Int_t nvar=28;// (const Int_t)Background->GetNbinsX ();
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
	
    const Int_t nvar = 28;
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
	
    const Int_t nvar = 28;
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



