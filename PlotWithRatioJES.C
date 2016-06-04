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

//TH1::SetDefaultSumw2();

TGraphErrors *h1= new TGraphErrors;
TGraphErrors *h2= new TGraphErrors;
TGraphErrors *h1up= new TGraphErrors;
TGraphErrors *h1down = new TGraphErrors;
TGraphErrors *h1EE = new TGraphErrors;
TGraphErrors *h2EE = new TGraphErrors;
TGraphErrors *h1EEup = new TGraphErrors;
TGraphErrors *h1EEdown = new TGraphErrors;
TGraphErrors *h1MM = new TGraphErrors;
TGraphErrors *h2MM = new TGraphErrors;
TGraphErrors *h1MMup = new TGraphErrors;
TGraphErrors *h1MMdown = new TGraphErrors;
TGraphErrors *h1GJ = new TGraphErrors;
TGraphErrors *h2GJ = new TGraphErrors;
TGraphErrors *h1GJup = new TGraphErrors;
TGraphErrors *h1GJdown = new TGraphErrors;
void testRatioEEMuMuComparison(TString, TString );
void testRatioGammaJets(TString , TString , TString );
TGraphAsymmErrors staterror(TH1F *, TH1F * );
TGraphAsymmErrors JESerror(TH1F *, TH1F *, TH1F *, TH1F * );

void PlotWithRatioEEMuMuComparison(TCanvas *c1, TGraphErrors *mcdyee,  TGraphErrors *datadyee,  TGraphErrors *mcdymm,  TGraphErrors *datadymm,TGraphErrors *mcdyeeup, TGraphErrors *mcdyeedown,  TGraphErrors *mcdymmup, TGraphErrors *mcdymmdown,  TString EjeX, TString histogramaname) {

TH1::SetDefaultSumw2();
  c1->cd(); 
   //Since all the histograms should have the same binning, take the info from one of them
   double xmin=mcdyee->GetX()[0];
   double xmax=mcdyee->GetX()[mcdyee->GetN()-1];
   double nbins=mcdyee->GetN();
  
   TH1F * hdatadyee= new TH1F("","", nbins,xmin,xmax);
   TH1F * hmcdyee= new TH1F("hmcdyee","hmcdyee", nbins,xmin,xmax);
   TH1F * hmcdyeeup= new TH1F("hmcdyeeup","hmcdyeeup", nbins,xmin,xmax);
   TH1F * hmcdyeedown= new TH1F("hmcdyeedown","hmcdyeedown", nbins,xmin,xmax);
   TH1F * hdatadymm= new TH1F("hdatadymm","hdatadymm", nbins,xmin,xmax);
   TH1F * hmcdymm= new TH1F("hmcdymm","hdatadymm", nbins,xmin,xmax);
   TH1F * hmcdymmup= new TH1F("hmcdymmup","hmcdymmup", nbins,xmin,xmax);
   TH1F * hmcdymmdown= new TH1F("hmcdymmdown","hmcdymmdown", nbins,xmin,xmax) ; 
   for (int i = 0, n = datadyee->GetN(); i < n; ++i) {
       
       hdatadyee->SetBinContent(i+1,datadyee->GetY()[i]);
       hmcdyee->SetBinContent(i+1,mcdyee->GetY()[i]);
       hmcdyeedown->SetBinContent(i+1,mcdyeedown->GetY()[i]);
       hmcdyeeup->SetBinContent(i+1,mcdyeeup->GetY()[i]);
 	   hdatadyee->SetBinError(i+1,datadyee->GetErrorY(i));
       hmcdyee->SetBinError(i+1,mcdyee->GetErrorY(i));          
       hdatadymm->SetBinContent(i+1,datadymm->GetY()[i]);
       hmcdymm->SetBinContent(i+1,mcdymm->GetY()[i]);
       hmcdymmdown->SetBinContent(i+1,mcdymmdown->GetY()[i]);
       hmcdymmup->SetBinContent(i+1,mcdymmup->GetY()[i]);
       hdatadymm->SetBinError(i+1,datadymm->GetErrorY(i));
       hmcdymm->SetBinError(i+1,mcdymm->GetErrorY(i));                  
    }

    hdatadyee->SetMarkerStyle(22);
    hdatadyee->SetMarkerColor(46);
    hdatadyee->SetLineColor(46);
    hmcdyee->SetMarkerColor(46);
    hmcdyee->SetMarkerStyle(22); 
    hdatadymm->SetMarkerStyle(23);
    hdatadymm->SetMarkerColor(4);
    hdatadymm->SetLineColor(4);
    hmcdymm->SetMarkerColor(4);
    hmcdymm->SetMarkerStyle(23);
    TH1F *ratiodyee=(TH1F*)hdatadyee->Clone(); 
    TH1F *ratiodymm=(TH1F*)hdatadymm->Clone(); 
   
    ratiodyee->Divide(hmcdyee);
    ratiodymm->Divide(hmcdymm);
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.014);
    pad1->SetLeftMargin(0.1);
    pad1->Draw();
    pad1->cd();
    hdatadyee->GetYaxis()->SetRangeUser(0.5, 1.5);
    if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.8,1.2);
    hdatadyee->GetXaxis()->SetLabelFont(43); //font in pixels
    hdatadyee->GetXaxis()->SetLabelSize(1); //in pixels
    hdatadyee->GetYaxis()->SetLabelFont(43);
    hdatadyee->GetYaxis()->SetLabelSize(16);
    hdatadyee->GetXaxis()->SetTitleFont(43); //font in pixels
    hdatadyee->GetXaxis()->SetTitleSize(20); //in pixels  
    hdatadyee->GetYaxis()->SetTitleFont(43); //font in pixels
    hdatadyee->GetYaxis()->SetTitleSize(20); //in pixels
    hdatadyee->GetXaxis()->SetTitleOffset(2.7);
    //hdatadyee->GetYaxis()->SetTitleOffset(1.9);
    hdatadyee->SetFillColor(54);
    hdatadyee->SetFillStyle(3345);                                                                                                                                    

    hdatadyee->Draw("e1");
    pad1->Update();
    TString ylabel="";
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hdatadyee->GetXaxis ()->GetXmin (), 1, hdatadyee->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw();
    }                                                                                                                                                                                                                        

if (histogramaname=="met_uPara_over_qt") ylabel="-<u_{||}/q_{T}> ";
if (histogramaname=="met_uPara_vs_zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerp_vs_zll_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_met_sumEt") ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_met_sumEt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_nVert")ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_nVert")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname.Contains("uPara_vs")) hdatadyee->GetYaxis()->SetRangeUser(0.,40);
if (histogramaname.Contains("uPerp")) hdatadyee->GetYaxis()->SetRangeUser(9.7,30);
if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.6,1.2);  
   
    hdatadymm->Draw("e1, X0 same");
    hdatadyee->Draw("AXIS");
    hdatadyee->Draw("e1, X0 same");
    hdatadymm->Draw("e1, X0 same");
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hdatadyee->GetXaxis ()->GetXmin (), 1, hdatadyee->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw();
    }                                                                                                              
    
    TLegend* leg(0);
    leg = new TLegend(0.15,0.65,0.35,0.85);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->AddEntry(hdatadymm, "Z#rightarrow #mu#mu", "lp");
    leg->AddEntry(hdatadyee, "Z#rightarrow ee", "lp");

    leg->Draw();
    c1->Modified();
                                                                   
    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.04);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    ratiodyee->GetYaxis()->SetTitle("Data/MC");
    ratiodyee->GetXaxis()->SetTitle(EjeX);
    ratiodyee->GetYaxis()->SetRangeUser(0.5, 1.5);
    if (histogramaname.Contains("over")) ratiodyee->GetYaxis()->SetRangeUser(0.8,1.2);
    ratiodyee->SetMarkerStyle(22);
    ratiodyee->SetMarkerColor(46);
    ratiodyee->SetLineColor(46);
    ratiodymm->SetMarkerStyle(23);
    ratiodymm->SetMarkerColor(4);
    ratiodymm->SetLineColor(4); 
    ratiodyee->GetXaxis()->SetLabelFont(43); //font in pixels
    ratiodyee->GetXaxis()->SetLabelSize(16); //in pixels
    ratiodyee->GetYaxis()->SetLabelFont(43);
    ratiodyee->GetYaxis()->SetLabelSize(16);
    ratiodyee->GetXaxis()->SetTitleFont(43); //font in pixels
    ratiodyee->GetXaxis()->SetTitleSize(20); //in pixels  
    ratiodyee->GetYaxis()->SetTitleFont(43); //font in pixels
    ratiodyee->GetYaxis()->SetTitleSize(20); //in pixels
    ratiodyee->GetXaxis()->SetTitleOffset(2.7);
    ratiodyee->GetYaxis()->SetTitleOffset(1.5);
    ratiodyee->SetFillColor(54);
    ratiodyee->SetFillStyle(3345);                                                                
  
    //==========  JES ============================================                                                                                                            
    
    //============================================================
    TGraphAsymmErrors erree= JESerror(hmcdyee,hmcdyeedown,hmcdyeeup,hdatadyee);
    TGraphAsymmErrors staterree= staterror(hmcdyee,hdatadyee);

   TAxis *erree_yaxis = erree.GetYaxis ();
   erree_yaxis->SetRangeUser (0, 2);
   erree.SetTitle (0);
   
   erree.SetFillColor (kRed);
   erree.SetFillStyle (3003);
   
   TAxis *staterree_yaxis = staterree.GetYaxis ();
   staterree_yaxis->SetRangeUser (0, 3);
   staterree.SetTitle (0);
   staterree.SetFillColor (kGreen);
   staterree.SetFillStyle (3003);  
   ratiodyee->Draw("e1, X0"); //"e2"
   ratiodymm->Draw("e1, X0 same"); //"e2"
   erree.Draw ("2 same");
   staterree.Draw("2 same");
   
   TLine *lineR =  new TLine ( ratiodyee->GetXaxis ()->GetXmin (), 1, ratiodyee->GetXaxis ()->GetXmax (), 1);
   lineR->SetLineColor (kBlue + 1); lineR->SetLineWidth (2); lineR->SetLineStyle (2); lineR->Draw();     
      
   TLegend* legratio(0);
   legratio = new TLegend(0.12,0.71,0.34,0.93);
   legratio->SetFillColor(0);
   legratio->SetBorderSize(0);
   legratio->AddEntry(&staterree, "Stat","f");
   legratio->AddEntry(&erree, "JES+Stat","f");
   legratio->Draw();

   pad1->cd();
   pad1->RedrawAxis();
   TLatex latex;                  
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextFont(42);
   latex.SetTextAlign(31);
   latex.SetTextSize(0.05);
   latex.DrawLatex(0.18, 0.91, "CMS");

   TLatex latexb;
   latexb.SetNDC();
   latexb.SetTextAngle(0);
   latexb.SetTextFont(42);
   latexb.SetTextAlign(31);
   latexb.SetTextSize(0.04);
   latexb.DrawLatex(0.34, 0.91, "Preliminary");
                                                   
   TLatex latexc;
   latexc.SetNDC();
   latexc.SetTextAngle(0);
   latexc.SetTextFont(42);
   latexc.SetTextAlign(31);
   latexc.SetTextSize(0.05);
   latexc.DrawLatex(0.90, 0.91, " 2.3 fb^{-1} (13 TeV)");              
   
   TLatex latexd;
   latexd.SetNDC();
   latexd.SetTextAngle(0);
   latexd.SetTextFont(42);
   latexd.SetTextAlign(31);
   latexd.SetTextSize(0.05);
   latexd.SetTextAngle(90);      
   latexd.DrawLatex(0.05, 0.87, ylabel);      


   pad1->Modified();
   pad2->Modified();
   pad1->Update();
   pad2->Update();
   pad1->cd();
   c1->Update();
   c1->cd();
   c1->Print("~/www/met/comparisons/May12/EEMuMu/"+histogramaname+".png");
   c1->Print("~/www/met/comparisons/May12/EEMuMu/"+histogramaname+".root");


}

void PlotWithRatio(TCanvas *c1, TGraphErrors *mcdy,  TGraphErrors *datady, TGraphErrors *mcdyup, TGraphErrors *mcdydown,  TString EjeX, TString histogramaname, TString channel) { 

TH1::SetDefaultSumw2();

  c1->cd(); 

  double xmin=mcdy->GetX()[0];
  double xmax=mcdy->GetX()[mcdy->GetN()-1];
  double nbins=mcdy->GetN();
  
  TH1F * hmcdy= new TH1F("","", nbins,xmin,xmax);
  TH1F * hdatady= new TH1F("","", nbins,xmin,xmax);
  TH1F * hmcdyup= new TH1F("hmcdyup","hmcdyup", nbins,xmin,xmax);
  TH1F * hmcdydown= new TH1F("hmcdydown","hmcdydown", nbins,xmin,xmax);

    for (int i = 0, n = datady->GetN(); i < n; ++i) {
            hdatady->SetBinContent(i+1,datady->GetY()[i]);
            hmcdy->SetBinContent(i+1,mcdy->GetY()[i]);
            hmcdydown->SetBinContent(i+1,mcdydown->GetY()[i]);
            hmcdyup->SetBinContent(i+1,mcdyup->GetY()[i]);
      	    hdatady->SetBinError(i+1,datady->GetErrorY(i));
            hmcdy->SetBinError(i+1,mcdy->GetErrorY(i));
    }

    TH1F *ratiody=(TH1F*)hdatady->Clone(); 
    ratiody->Divide(hmcdy);
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.011);
    
    pad1->Draw();
    pad1->cd();
    if (histogramaname.Contains("over")) hdatady->GetYaxis()->SetRangeUser(0.6,1.2);
    hdatady->SetMarkerStyle(22);
    hdatady->SetMarkerColor(46);
    hdatady->SetLineColor(46);
    hmcdy->SetMarkerStyle(23);
    hmcdy->SetMarkerColor(4);
    hmcdy->SetLineColor(4);    
    hmcdy->GetXaxis()->SetLabelFont(43); //font in pixels
    hmcdy->GetXaxis()->SetLabelSize(1); //in pixels
    hmcdy->GetYaxis()->SetLabelFont(43);
    hmcdy->GetYaxis()->SetLabelSize(16);
    hmcdy->GetXaxis()->SetTitleFont(43); //font in pixels
    hmcdy->GetXaxis()->SetTitleSize(20); //in pixels  
    hmcdy->GetYaxis()->SetTitleFont(43); //font in pixels
    hmcdy->GetYaxis()->SetTitleSize(20); //in pixels
    hmcdy->GetXaxis()->SetTitleOffset(3.5);
    hmcdy->GetYaxis()->SetTitleOffset(1.5);
    hmcdy->SetFillColor(54);
    hmcdy->SetFillStyle(3345);                                                                    
    
    datady->GetHistogram()->GetYaxis()->SetTitle(EjeX);
    
    TString ylabel="";
if (histogramaname=="met_uPara_over_qt") ylabel="-<u_{||}/q_{T}> [GeV]";
if (histogramaname=="met_uPara_vs_zll_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerp_vs_zll_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_met_sumEt-zll_pt") ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_met_sumEt-zll_pt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_vs_nVert")ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="met_uPerp_vs_nVert")ylabel="#sigma ( u_{#perp}  ) [GeV]";            
if (histogramaname=="gamma_met_uPara_vs_gamma_pt") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="gamma_met_uPerp_vs_gamma_pt") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="gamma_met_uPara_vs_met_sumEt") ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="gamma_met_uPerp_vs_met_sumEt")ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="gamma_met_uPara_vs_nVert")ylabel="#sigma ( u_{||}  ) [GeV]";
if (histogramaname=="gamma_met_uPerp_vs_nVert")ylabel="#sigma ( u_{#perp}  ) [GeV]";         
if (histogramaname=="met_uPara_over_qt") ylabel="-<u_{||}/q_{T}> [GeV]";
    
    hmcdy->GetYaxis()->SetTitle(ylabel);
    hmcdy->GetYaxis()->SetRangeUser(9.5,40);
    if (histogramaname.Contains("uPerp")) hmcdy->GetYaxis()->SetRangeUser(9.7,25);
    if (histogramaname.Contains("sumEt")) hmcdy->GetYaxis()->SetRangeUser(0,40);
    if (histogramaname.Contains("nVert")) hmcdy->GetYaxis()->SetRangeUser(0,40);
    if (histogramaname.Contains("over")) hmcdy->GetYaxis()->SetRangeUser(0.5,1.3);
    hmcdy->Draw("e1, X0 ");
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hmcdy->GetXaxis ()->GetXmin (), 1, hmcdy->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw(); 
    }                                                                                                   
    hdatady->Draw("e1, X0 same");
    hmcdy->Draw("e1, X0 same");
    TLegend* leg(0);
    leg = new TLegend(0.15,0.65,0.35,0.85);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    if (channel=="M"){
        leg->AddEntry(hmcdy, "Z#rightarrow #mu#mu MC", "lp");
        leg->AddEntry(hdatady, "Z#rightarrow #mu#mu Data", "lp");
    }
    if (channel=="E"){
        leg->AddEntry(hmcdy, "Z#rightarrow ee MC", "lp");
        leg->AddEntry(hdatady, "Z#rightarrow ee Data", "lp");
    }                                                          
                                                                     
    if (channel=="G"){
        leg->AddEntry(hmcdy, "#gamma + jets", "lp");     
        leg->AddEntry(hmcdy, "#gamma + jets", "lp");
    }                                                          

    //leg->AddEntry(hdatady, "Data", "lp");
    leg->Draw();
    c1->Modified();
                                                                   
    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.04);
    pad2->SetBottomMargin(0.2);

    pad2->Draw();
    pad2->cd();
    
    ratiody->GetYaxis()->SetTitle("Data/MC");
    ratiody->GetXaxis()->SetTitle(EjeX);
    ratiody->SetTitle("");
    ratiody->GetYaxis()->SetRangeUser(0.5, 1.5);
    //ratiody->GetXaxis()->SetRangeUser(xmin, xmax);
    if (histogramaname.Contains("over")) ratiody->GetYaxis()->SetRangeUser(0.9,1.1);
    ratiody->SetMarkerStyle(22);
    ratiody->SetMarkerColor(4);
    ratiody->SetLineColor(4);
    ratiody->GetXaxis()->SetLabelFont(43); //font in pixels
    ratiody->GetXaxis()->SetLabelSize(16); //in pixels
    ratiody->GetYaxis()->SetLabelFont(43);
    ratiody->GetYaxis()->SetLabelSize(16);
    ratiody->GetXaxis()->SetTitleFont(43); //font in pixels
    ratiody->GetXaxis()->SetTitleSize(20); //in pixels  
    ratiody->GetYaxis()->SetTitleFont(43); //font in pixels
    ratiody->GetYaxis()->SetTitleSize(20); //in pixels
    ratiody->GetXaxis()->SetTitleOffset(3.2);
    ratiody->GetYaxis()->SetTitleOffset(1.5);
    ratiody->SetFillColor(54);
    ratiody->SetFillStyle(3345);
  
    
    //==========  JES ============================================                                                                                                            
    
    
    //============================================================
    
    TGraphAsymmErrors err= JESerror(hmcdy,hmcdydown,hmcdyup,hdatady);
    TGraphAsymmErrors staterr= staterror(hmcdy,hdatady);

    
    TAxis *err_yaxis = err.GetYaxis ();

    err_yaxis->SetRangeUser (0, 2);
    err.SetTitle (0);
    
    err.SetFillColor (kRed);
    err.SetFillStyle (3003);
    
    TAxis *staterr_yaxis = staterr.GetYaxis ();
    staterr_yaxis->SetRangeUser (0, 3);
    staterr.SetTitle (0);

    staterr.SetFillColor (kGreen);
    staterr.SetFillStyle (3003);
    
    ratiody->Draw("e1 X0"); //"e2"
    err.Draw ("2 same");
    staterr.Draw("2 same");
     
    ratiody->SetMarkerSize (0.8);
    ratiody->SetMarkerStyle (20);
    TLine *lineR =  new TLine ( ratiody->GetXaxis ()->GetXmin (), 1, ratiody->GetXaxis ()->GetXmax (), 1);

    lineR->SetLineColor (kBlue + 1); lineR->SetLineWidth (2); lineR->SetLineStyle (2); lineR->Draw();     
        
    TLegend* legratio(0);
    legratio = new TLegend(0.12,0.73,0.32,0.93);
    legratio->SetFillColor(0);
    legratio->SetBorderSize(0);
    legratio->AddEntry(&staterr, "Stat","f");
    legratio->AddEntry(&err, "JES+Stat","f");
    legratio->Draw();

    pad1->cd();
    TLatex latex;                  
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.18, 0.91, "CMS");

    TLatex latexb;
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(31);
    latexb.SetTextSize(0.04);
    latexb.DrawLatex(0.34, 0.91, "Preliminary");
                                                    
    TLatex latexc;
    latexc.SetNDC();
    latexc.SetTextAngle(0);
    latexc.SetTextFont(42);
    latexc.SetTextAlign(31);
    latexc.SetTextSize(0.05);
    latexc.DrawLatex(0.90, 0.91, " 0.6 fb^{-1} (13 TeV)");              
    if ((histogramaname.Contains("sumEt"))  or (histogramaname.Contains("nVert")) ) {   
        TLatex latexz;                                                                  
        latexz.SetNDC();                                                                
        latexz.SetTextAngle(0);                                                         
        latexz.SetTextFont(42);                                                         
        latexz.SetTextAlign(31);                                                        
        latexz.SetTextSize(0.04);                                                       
        latexz.DrawLatex(0.38, 0.60, "Z p_{T} > 50 GeV");                        
    }                                                                                   
    TLatex latexe;                                                   
    latexe.SetNDC();
    latexe.SetTextAngle(0);
    latexe.SetTextFont(42);
    latexe.SetTextAlign(31);
    latexe.SetTextSize(0.05);
    latexe.DrawLatex(0.72, 0.78, "No JEC Residuals");       
    pad1->Modified();
    pad2->Modified();
    pad1->Update();
    pad2->Update();

    c1->Update();
    c1->cd();
    if (channel == "G") c1->Print("~/www/met/comparisons/May11/GJets/"+histogramaname+".png");
    if (channel == "G") c1->Print("~/www/met/comparisons/May11/GJets/"+histogramaname+".root");
    if (channel == "E") c1->Print("~/www/met/comparisons/June02/EE/"+histogramaname+".png");
    if (channel == "E") c1->Print("~/www/met/comparisons/June02/EE/"+histogramaname+".root");
    if (channel == "M") c1->Print("~/www/met/comparisons/June02/MM/"+histogramaname+".png");
    if (channel == "M") c1->Print("~/www/met/comparisons/June02/MM/"+histogramaname+".root");


}




void testRatio(TString histograma,TString histograma2, TString EjeX, TString channel) {                                                 

cout << " channel" <<  channel << endl;
TCanvas *c1 = new TCanvas("c1","example",600,700);
if (channel == "MM"){
TFile file1 ("Not_BKG_SubtractionDYtgraphsMZll80X.root"); // MC DY
TFile file2 ("Not_BKG_SubtractionDatatgraphsMZll80X.root");// data DY
TFile file1up ("Not_BKG_SubtractionDY_up_tgraphs_jes_DYMZll80X.root"); // MC DY up
TFile file1down ("Not_BKG_SubtractionDY_down_tgraphs_jes_DYMZll80X.root"); // MC DY down
 h1 = (TGraphErrors*) file1.Get("M_"+ histograma2); //mc mumu 
 h2 = (TGraphErrors*) file2.Get("M_"+ histograma); //mc ee 
 h1up = (TGraphErrors*)  file1up.Get("M_"+ histograma2); //mc mumu 
 h1down = (TGraphErrors*) file1down.Get("M_"+ histograma2); //mc mumu

   gStyle->SetOptStat(0);
   gStyle->SetErrorX(0.5);
   PlotWithRatio(c1, h1, h2, h1up, h1down, EjeX, histograma, "M");
}
if (channel == "EE"){
TFile file1 ("Not_BKG_SubtractionDYtgraphsEZll80X.root"); // MC DY
TFile file2 ("Not_BKG_SubtractionDatatgraphsEZll80X.root");// data DY
TFile file1up ("Not_BKG_SubtractionDY_up_tgraphs_jes_DYEZll80X.root"); // MC DY up
TFile file1down ("Not_BKG_SubtractionDY_down_tgraphs_jes_DYEZll80X.root"); // MC DY down
 h1 = (TGraphErrors*) file1.Get("E_"+ histograma2); //mc mumu 
 h2 = (TGraphErrors*) file2.Get("E_"+ histograma); //mc ee 
 h1up = (TGraphErrors*)  file1up.Get("E_"+ histograma2); //mc mumu 
 h1down = (TGraphErrors*) file1down.Get("E_"+ histograma2); //mc mumu

   gStyle->SetOptStat(0);
   gStyle->SetErrorX(0.5);
   PlotWithRatio(c1, h1, h2, h1up, h1down, EjeX, histograma, "E");
}
if (channel == "GG"){
TFile file1 ("Not_BKG_SubtractionGJetstgraphsgammagjets.root"); // MC DY
TFile file2 ("Not_BKG_SubtractionDatatgraphsgammagjets.root");// data DY
TFile file1up ("Not_BKG_SubtractionGJets_up_tgraphs_jes_GJetsgammagjets.root"); // MC DY up
TFile file1down ("Not_BKG_SubtractionGJets_down_tgraphs_jes_GJetsgammagjets.root"); // MC DY down
 h1 = (TGraphErrors*) file1.Get(histograma); //mc mumu 
 h2 = (TGraphErrors*) file2.Get(histograma); //mc ee 
 h1up = (TGraphErrors*)  file1up.Get(histograma); //mc mumu 
 h1down = (TGraphErrors*) file1down.Get(histograma); //mc mumu

   gStyle->SetOptStat(0);
   gStyle->SetErrorX(0.5);
   PlotWithRatio(c1, h1, h2, h1up, h1down, EjeX, histograma, "G");
}

}
                                                                                                                                                              


void testRatioEEMuMuComparison(TString histograma, TString EjeX) {                                               


TCanvas *c1 = new TCanvas("c1","example",600,700);

TFile file1EE ("Not_BKG_SubtractionDYtgraphsEZllNEW.root");
TFile file2EE ("Not_BKG_SubtractionDatatgraphsEZllNEW.root");
TFile file1MM ("Not_BKG_SubtractionDYtgraphsMZllNEW.root"); 
TFile file2MM ("Not_BKG_SubtractionDatatgraphsMZllNEW.root");
TFile file1EEup ("Not_BKG_SubtractionDY_up_tgraphs_jes_DYEZllNEW.root"); 
TFile file1EEdown ("Not_BKG_SubtractionDY_down_tgraphs_jes_DYEZllNEW.root"); 
TFile file1MMup ("Not_BKG_SubtractionDY_up_tgraphs_jes_DYMZllNEW.root"); 
TFile file1MMdown ("Not_BKG_SubtractionDY_down_tgraphs_jes_DYMZllNEW.root"); 

 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 

 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1MMup = (TGraphErrors*)   file1MMup.Get("M_"+ histograma); 
 h1MMdown = (TGraphErrors*) file1MMdown.Get("M_"+ histograma);   
  
    cout << " h1EEup " << h1EEup << endl ;
    cout << " h1EEdown " << h1EEdown << endl ;
    cout << " h1EE " << h1EE << endl ;
    cout << " h2EE " << h2EE << endl ;            
    cout << " h1MMup " << h1MMup << endl ;
    cout << " h1MMdown " << h1MMdown << endl ;       
    cout << " h1MM " << h1MM << endl ;               
    cout << " h2MM " << h2MM << endl ;               
    gStyle->SetOptStat(0);
   gStyle->SetErrorX(0.5);
   PlotWithRatioEEMuMuComparison(c1, h1EE, h2EE, h1MM, h2MM, h1EEup, h1EEdown, h1MMup, h1MMdown, EjeX, histograma);
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


