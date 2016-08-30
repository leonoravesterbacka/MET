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
TGraphAsymmErrors staterror(TH1F *, TH1F * );
TGraphAsymmErrors JESerror(TH1F *, TH1F *, TH1F *, TH1F * );
TGraphAsymmErrors TOTerror(TH1F *, TH1F *, TH1F *, TH1F *, TH1F *, TH1F *  );

void PlotWithRatioEEMuMuGJetsComparison(TCanvas *c1, TGraphErrors *mcdyee,  TGraphErrors *datadyee,  TGraphErrors *mcdymm,  TGraphErrors *datadymm,TGraphErrors *mcdygg,  TGraphErrors *datadygg,  TGraphErrors *mcdyeeup, TGraphErrors *mcdyeedown, TGraphErrors *mcdyeeupUncl, TGraphErrors *mcdyeedownUncl,  TString EjeX, TString histogramaname) {

TH1::SetDefaultSumw2();
  c1->cd(); 
  Float_t xbins[] = { 0, 15, 30, 51, 65, 80, 110, 140,170, 200,250, 330};
    double nbins=11.;
 //     cout << "plotting resolution as a function of qt" << endl;
 //   if (histogramaname.Contains("over")){
//    Float_t xbins[] = { 20,28,36, 44,52,60, 68,76, 84,92,100, 120, 150, 175, 200, 225,250, 275, 305, 335, 365, 395, 430, 500}; 
 //   double nbins=23.;
 //     cout << "plotting scale " << endl;
 //   }
///    else if (histogramaname.Contains("nVert")){
 //  Float_t xbins[] = { 0, 6, 8, 10, 12 ,14, 16, 20,  40};
 //     double nbins=8.;
  //    cout << "plotting resolution as a function of nVert" << endl;
  //  }                                                                                                                                         



   TH1F * hdatadyee= new TH1F("","", nbins,xbins);
   TH1F * hmcdyee= new TH1F("hmcdyee","hmcdyee", nbins,xbins);
   TH1F * hmcdyeeup= new TH1F("hmcdyeeup","hmcdyeeup", nbins, xbins);
   TH1F * hmcdyeedown= new TH1F("hmcdyeedown","hmcdyeedown", nbins,xbins);
   TH1F * hmcdyeeupUncl= new TH1F("hmcdyeeupUncl","hmcdyeeupUncl", nbins, xbins);
   TH1F * hmcdyeedownUncl= new TH1F("hmcdyeedownUncl","hmcdyeedownUncl", nbins, xbins);
   TH1F * hdatadymm= new TH1F("hdatadymm","hdatadymm", nbins,xbins);
   TH1F * hmcdymm= new TH1F("hmcdymm","hdatadymm", nbins,xbins);
   TH1F * hdatadygg= new TH1F("hdatadygg","hdatadygg", nbins,xbins);
   TH1F * hmcdygg= new TH1F("hmcdygg","hdatadygg", nbins,xbins);
   
   for (int i = 0, n = nbins; i < n; ++i) {
       hdatadyee->SetBinContent(i+1,datadyee->GetY()[i]);
       hmcdyee->SetBinContent(i+1,mcdyee->GetY()[i]);
       hmcdyeedown->SetBinContent(i+1,mcdyeedown->GetY()[i]);
       hmcdyeeup->SetBinContent(i+1,mcdyeeup->GetY()[i]);
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

if (histogramaname=="met_uPara_over_zll_pt") ylabel="   -<u_{||}>/<q_{T}> ";
if (histogramaname=="met_metx_over_nVert") ylabel="   <E_{X}^{miss}> ";
if (histogramaname=="met_mety_over_nVert") ylabel="   <E_{Y}^{miss}> ";
if (histogramaname=="met_uParavs_zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerpvs_zll_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_met_sumEt-zll_pt") ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_met_sumEt-zll_pt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_nVert")ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_nVert")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_metx_vs_nVert")ylabel="#sigma ( E^{miss} X) [GeV]";
if (histogramaname=="met_mety_vs_nVert")ylabel="#sigma ( E^{miss} Y) [GeV]";
if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.58,1.2);  
   
    hdatadymm->Draw("e1, X0 same");
    hdatadyee->Draw("AXIS");
    hdatadygg->Draw("e1, X0 same");
    hdatadymm->Draw("e1, X0 same");
    hdatadyee->Draw("e1, X0 same");
    //hdatadyee->Draw("e1, X0 same");
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hdatadyee->GetXaxis ()->GetXmin (), 1, hdatadyee->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw();
    }                                                                                                              
    
    TLegend* leg(0);
    leg = new TLegend(0.15,0.65,0.35,0.85);
    if (histogramaname.Contains("over")){
        leg = new TLegend(0.45,0.15,0.85,0.35);
    }
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillColor(0);
    leg->AddEntry(hdatadymm, "Data Z#rightarrow #mu#mu", "lp");
    //leg->AddEntry(hmcdymm, "MC Z#rightarrow #mu#mu ", "lp");
    leg->AddEntry(hdatadyee, "Data Z#rightarrow ee", "lp");
    leg->AddEntry(hdatadygg, "Data #gamma, q_{T} > 50 GeV", "lp");

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
    ratiodyee->GetYaxis()->SetRangeUser(0.7, 1.3);
    if (histogramaname.Contains("over")) 
        ratiodyee->GetYaxis()->SetRangeUser(0.9,1.1);
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
    TGraphAsymmErrors toterree= TOTerror(hmcdyee,hmcdyeedown,hmcdyeeup, hmcdyeedownUncl,hmcdyeeupUncl,hdatadyee);
    TGraphAsymmErrors erree= JESerror(hmcdyee,hmcdyeedown,hmcdyeeup,hdatadyee);
    TGraphAsymmErrors staterree= staterror(hmcdyee,hdatadyee);

   TAxis *erree_yaxis = erree.GetYaxis ();
   erree_yaxis->SetRangeUser (0, 2);
   erree.SetTitle (0);
   erree.SetFillColor (kRed);
   erree.SetFillStyle (3002);                                   
   
   TAxis *staterree_yaxis = staterree.GetYaxis ();
   staterree_yaxis->SetRangeUser (0, 2);
   staterree.SetTitle (0);
   staterree.SetFillColor (kGreen);
   staterree.SetFillStyle (3002);  
   
   TAxis *toterree_yaxis = toterree.GetYaxis ();
   toterree_yaxis->SetRangeUser (0, 2);
   toterree.SetTitle (0);
   toterree.SetFillColor (kGray+3);
   toterree.SetFillStyle (3002);                                   
   
   ratiodyee->Draw("e1, X0"); //"e2"
   ratiodygg->Draw("e1, X0 same "); //"e2"
   ratiodymm->Draw("e1, X0 same"); //"e2"
   ratiodyee->Draw("e1, X0 same"); //"e2"
   toterree.Draw ("2 same");
  // erree.Draw ("2 same");
  // staterree.Draw("2 same");
   TLine *lineR =  new TLine ( ratiodyee->GetXaxis ()->GetXmin (), 1, ratiodyee->GetXaxis ()->GetXmax (), 1);
   lineR->SetLineColor (kBlue + 1); lineR->SetLineWidth (2); lineR->SetLineStyle (2); lineR->Draw();     
        
    TLegend* legratio(0);
   
    legratio = new TLegend(0.12,0.8,0.8,0.95);
    if (histogramaname.Contains("over")){   
        legratio = new TLegend(0.2,0.8,0.89,0.95);
    }
    legratio->SetFillColor(0);
    legratio->SetBorderSize(0);
    legratio-> SetNColumns(3);
    legratio->AddEntry(&toterree, "Uncl E_{T}^{miss}+ JES + Stat","f");
    legratio->AddEntry(&erree, "JES + Stat","f");
    legratio->AddEntry(&staterree, "Stat","f");
    //legratio->Draw();
//    ratiodygg->Draw("e1, X0 same "); //"e2"
//    ratiodymm->Draw("e1, X0 same"); //"e2"
//    ratiodyee->Draw("e1, X0 same"); //"e2"

    pad1->cd();
    TLatex latex;                  
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.07);
    latex.DrawLatex(0.21, 0.91, "#bf{CMS}");

//    if ((histogramaname.Contains("sumEt"))  or (histogramaname.Contains("nVert")) ) {
        TLatex latexz;
        latexz.SetNDC();
        latexz.SetTextAngle(0);
        latexz.SetTextFont(42);
        latexz.SetTextAlign(31);
        latexz.SetTextSize(0.05);
//   if  (histogramaname.Contains("over")){ 
//        latexz.DrawLatex(0.78, 0.1, "p_{T}^{#gamma} > 50 GeV");
//    }
//    else {
//        latexz.DrawLatex(0.35, 0.58, "p_{T}^{#gamma} > 50 GeV");
//    }
    
    TLatex latexb;
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(31);
    latexb.SetTextSize(0.05);
    latexb.DrawLatex(0.41, 0.91, "#it{Preliminary}");
                                                    
    TLatex latexc;
    latexc.SetNDC();
    latexc.SetTextAngle(0);
    latexc.SetTextFont(42);
    latexc.SetTextAlign(31);
    latexc.SetTextSize(0.05);
    latexc.DrawLatex(0.90, 0.91, "12.9 fb^{-1} (13 TeV, 2016)");              
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
    c1->Print("~/www/met/comparisons/Aug04/"+histogramaname+".png");
    c1->Print("~/www/met/comparisons/Aug04/"+histogramaname+".root");
    c1->Print("~/www/met/comparisons/Aug04/"+histogramaname+".pdf");

}


void testRatioEEMuMuGJetsComparison(TString histograma, TString histograma2, TString EjeX) {                                               


TCanvas *c1 = new TCanvas("c1","example",600,700);
if (histograma.Contains("over")){
    
    
 TFile file1EE ("DYEZllScaleTestwoOF.root");
 TFile file2EE ("DataEZllScaleTestwoOF.root");
 TFile file1MM ("DYMZllScaleTestwoOF.root"); 
 TFile file2MM ("DataMZllScaleTestwoOF.root");
 TFile file1GG ("GJetsgjetsScalewoOF.root"); 
 TFile file2GG ("DatagjetsScalewoOF.root");
 TFile file1EEup     ("DY_up_jes_DYEZllScaleTestwoOF.root"); 
 TFile file1EEdown   ("DY_down_jes_DYEZllScaleTestwoOF.root"); 
 TFile file1EEupUncl ("DY_up_uncl_DYEZllScaleTestwoOF.root");
 TFile file1EEdownUncl ("DY_down_uncl_DYEZllScaleTestwoOF.root");
 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
 h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
 h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 

}
else{
    TFile file1EE ("DYEZllResolutionwoOF.root");
    TFile file2EE ("DataEZllResolutionwoOF.root");
    TFile file1MM ("DYMZllResolutionwoOF.root"); 
    TFile file2MM ("DataMZllResolutionwoOF.root");
    TFile file1GG ("GJetsgjetsResolutionwoOF.root");
    TFile file2GG ("DatagjetsResolutionwoOF.root");
    TFile file1EEup ("DY_up_jes_DYEZllResolutionwoOF.root"); 
    TFile file1EEdown ("DY_down_jes_DYEZllResolutionwoOF.root"); 
    TFile file1EEupUncl ("DY_up_uncl_DYEZllResolutionwoOF.root");
    TFile file1EEdownUncl ("DY_down_uncl_DYEZllResolutionwoOF.root"); 
    h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
    h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
    h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
    h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
    h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
    h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
    h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
    h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
    h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
    h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 

}
if (histograma.Contains("nVert")){
    TFile file1EE ("DYEZllResolutionXY.root");
    TFile file2EE ("DataEZllResolutionXY.root");
    TFile file1MM ("DYMZllResolutionXY.root"); 
    TFile file2MM ("DataMZllResolutionXY.root");
    TFile file1GG ("GJetsgjetsResolution.root");
    TFile file2GG ("DatagjetsResolution.root");
    TFile file1EEup ("DY_up_jes_DYEZllResolutionXY.root"); 
    TFile file1EEdown ("DY_down_jes_DYEZllResolutionXY.root"); 
    TFile file1EEupUncl ("DY_up_uncl_DYEZllResolutionXY.root");
    TFile file1EEdownUncl ("DY_down_uncl_DYEZllResolutionXY.root"); 
    h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
    h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
    h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
    h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
    h1GG = (TGraphErrors*) file1GG.Get(histograma2); 
    h2GG = (TGraphErrors*) file2GG.Get(histograma2); 
    h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
    h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
    h1EEupUncl = (TGraphErrors*)   file1EEupUncl.Get("E_"+ histograma); 
    h1EEdownUncl = (TGraphErrors*)   file1EEdownUncl.Get("E_"+ histograma); 

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
   gStyle->SetErrorX(0.5);

   PlotWithRatioEEMuMuGJetsComparison(c1, h1EE, h2EE, h1MM, h2MM, h1GG,h2GG, h1EEup, h1EEdown, h1EEupUncl, h1EEdownUncl, EjeX, histograma);
}


TGraphAsymmErrors JESerror(TH1F *Background, TH1F *Backgrounddown, TH1F *Backgroundup, TH1F * hdata){
    
    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");
	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	
    const Int_t nvar = 25;
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
	const Int_t nvar=25;// (const Int_t)Background->GetNbinsX ();
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


TGraphAsymmErrors TOTerror(TH1F *Background, TH1F *Backgrounddown, TH1F *Backgroundup,  TH1F *BackgroundUncldown, TH1F *BackgroundUnclup, TH1F * hdata){
    
    TH1F *den1 = (TH1F *) Background->Clone ("bkgden1");
	TH1F *den2 = (TH1F *) Background->Clone ("bkgden2");
	
    const Int_t nvar = 25;
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
                                                                                                                                                                                                                                                                                                                                                                                                 





