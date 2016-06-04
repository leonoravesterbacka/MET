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
TGraphErrors *h1GG = new TGraphErrors;
TGraphErrors *h2GG = new TGraphErrors;
TGraphErrors *h1GGup = new TGraphErrors;
TGraphErrors *h1GGdown = new TGraphErrors;
void testRatioEEMuMuGJetsComparison(TString, TString,TString );
TGraphAsymmErrors staterror(TH1F *, TH1F * );
TGraphAsymmErrors JESerror(TH1F *, TH1F *, TH1F *, TH1F * );

void PlotWithRatioEEMuMuGJetsComparison(TCanvas *c1, TGraphErrors *mcdyee,  TGraphErrors *datadyee,  TGraphErrors *mcdymm,  TGraphErrors *datadymm,TGraphErrors *mcdygg,  TGraphErrors *datadygg,  TGraphErrors *mcdyeeup, TGraphErrors *mcdyeedown,  TGraphErrors *mcdymmup, TGraphErrors *mcdymmdown,TGraphErrors *mcdyggup, TGraphErrors *mcdyggdown,   TString EjeX, TString histogramaname) {

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
   TH1F * hdatadygg= new TH1F("hdatadygg","hdatadygg", nbins,xmin,xmax);
   TH1F * hmcdygg= new TH1F("hmcdygg","hdatadygg", nbins,xmin,xmax);
   TH1F * hmcdyggup= new TH1F("hmcdyggup","hmcdyggup", nbins,xmin,xmax);
   TH1F * hmcdyggdown= new TH1F("hmcdyggdown","hmcdyggdown", nbins,xmin,xmax) ; 
   
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
       hdatadygg->SetBinContent(i+1,datadygg->GetY()[i]);
       hmcdygg->SetBinContent(i+1,mcdygg->GetY()[i]);
       hmcdyggdown->SetBinContent(i+1,mcdyggdown->GetY()[i]);
       hmcdyggup->SetBinContent(i+1,mcdyggup->GetY()[i]);
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
    hmcdymm->SetMarkerColor(4);
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
    hdatadyee->GetYaxis()->SetLabelFont(43);
    hdatadyee->GetYaxis()->SetLabelSize(16);
    hdatadyee->GetYaxis()->SetTitleFont(43); //fon
    hdatadyee->GetYaxis()->SetTitleSize(20); //in 
   
    hdatadymm->Draw("e1, X0 same");
    hdatadyee->Draw("AXIS");
    hdatadygg->Draw("e1, X0 same");
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
    leg->AddEntry(hdatadygg, "#gamma + jets", "lp");

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
    if (histogramaname.Contains("over")) ratiodyee->GetYaxis()->SetRangeUser(0.9,1.1);
    ratiodyee->SetMarkerStyle(22);
    ratiodyee->SetMarkerColor(46);
    ratiodyee->SetLineColor(46);
    ratiodymm->SetMarkerStyle(23);
    ratiodymm->SetMarkerColor(4);
    ratiodymm->SetLineColor(4); 
    ratiodygg->SetMarkerStyle(20);
    ratiodygg->SetMarkerColor(8);
    ratiodygg->SetLineColor(8);     
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
   ratiodygg->Draw("e1, X0 same "); //"e2"
   ratiodymm->Draw("e1, X0 same"); //"e2"
   ratiodyee->Draw("e1, X0 same"); //"e2"
   erree.Draw ("2 same");
   staterree.Draw("2 same");
   TLine *lineR =  new TLine ( ratiodyee->GetXaxis ()->GetXmin (), 1, ratiodyee->GetXaxis ()->GetXmax (), 1);
   lineR->SetLineColor (kBlue + 1); lineR->SetLineWidth (2); lineR->SetLineStyle (2); lineR->Draw();     
        
    TLegend* legratio(0);
    legratio = new TLegend(0.12,0.73,0.32,0.93);
    legratio->SetFillColor(0);
    legratio->SetBorderSize(0);
    legratio->AddEntry(&staterree, "Stat","f");
    legratio->AddEntry(&erree, "JES+Stat","f");
    legratio->Draw();

    pad1->cd();
    TLatex latex;                  
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.18, 0.91, "CMS");

    if ((histogramaname.Contains("sumEt"))  or (histogramaname.Contains("nVert")) ) {
        TLatex latexz;
        latexz.SetNDC();
        latexz.SetTextAngle(0);
        latexz.SetTextFont(42);
        latexz.SetTextAlign(31);
        latexz.SetTextSize(0.04);
        latexz.DrawLatex(0.38, 0.60, "Z/#gamma p_{T} > 50 GeV");
    }
    
    
    
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
    TLatex latexd;
    latexd.SetNDC();
    latexd.SetTextAngle(0);
    latexd.SetTextFont(42);
    latexd.SetTextAlign(31);
    latexd.SetTextSize(0.05);
    latexd.SetTextAngle(90); 
    
    
    TLatex latexe;                                          
    latexe.SetNDC();
    latexe.SetTextAngle(0);
    latexe.SetTextFont(42);
    latexe.SetTextAlign(31);
    latexe.SetTextSize(0.05);
    latexe.DrawLatex(0.72, 0.78, "No JEC Residuals");       




    latexd.DrawLatex(0.05, 0.87, ylabel);      
    pad1->Modified();
    pad2->Modified();
    pad1->Update();
    pad2->Update();
    pad1->cd();
    c1->Update();
    c1->cd();
    c1->Print("~/www/met/comparisons/June02/"+histogramaname+".png");
    c1->Print("~/www/met/comparisons/June02/"+histogramaname+".root");

}


void testRatioEEMuMuGJetsComparison(TString histograma, TString histogramagjets, TString EjeX) {                                               


TCanvas *c1 = new TCanvas("c1","example",600,700);

TFile file1EE ("Not_BKG_SubtractionDYtgraphsEZll80X.root");
TFile file2EE ("Not_BKG_SubtractionDatatgraphsEZll80X.root");
TFile file1MM ("Not_BKG_SubtractionDYtgraphsMZll80X.root"); 
TFile file2MM ("Not_BKG_SubtractionDatatgraphsMZll80X.root");
TFile file1GG ("Not_BKG_SubtractionGJetstgraphsgammagjets80X.root"); 
TFile file2GG ("Not_BKG_SubtractionDatatgraphsgammagjets80X.root");
TFile file1EEup ("Not_BKG_SubtractionDY_up_tgraphs_jes_DYEZll80X.root"); 
TFile file1EEdown ("Not_BKG_SubtractionDY_down_tgraphs_jes_DYEZll80X.root"); 
TFile file1MMup ("Not_BKG_SubtractionDY_up_tgraphs_jes_DYMZll80X.root"); 
TFile file1MMdown ("Not_BKG_SubtractionDY_down_tgraphs_jes_DYMZll80X.root"); 
TFile file1GGup ("Not_BKG_SubtractionGJets_up_tgraphs_jes_GJetsgammagjets80X.root"); 
TFile file1GGdown ("Not_BKG_SubtractionGJets_down_tgraphs_jes_GJetsgammagjets80X.root"); 

 h1EE = (TGraphErrors*) file1EE.Get("E_"+ histograma); 
 h2EE = (TGraphErrors*) file2EE.Get("E_"+ histograma); 
 h1MM = (TGraphErrors*) file1MM.Get("M_"+histograma); 
 h2MM = (TGraphErrors*) file2MM.Get("M_"+histograma); 
 h1GG = (TGraphErrors*) file1GG.Get(histogramagjets); 
 h2GG = (TGraphErrors*) file2GG.Get(histogramagjets); 

 h1EEup = (TGraphErrors*)   file1EEup.Get("E_"+ histograma); 
 h1EEdown = (TGraphErrors*) file1EEdown.Get("E_"+ histograma); 
 h1MMup = (TGraphErrors*)   file1MMup.Get("M_"+ histograma); 
 h1MMdown = (TGraphErrors*) file1MMdown.Get("M_"+ histograma);   
 h1GGup = (TGraphErrors*)   file1GGup.Get( histogramagjets); 
 h1GGdown = (TGraphErrors*) file1GGdown.Get(histogramagjets);   
  
    cout << " h1EEup " << h1EEup << endl ;
    cout << " h1EEdown " << h1EEdown << endl ;
    cout << " h1EE " << h1EE << endl ;
    cout << " h2EE " << h2EE << endl ;            
    cout << " h1MMup " << h1MMup << endl ;
    cout << " h1MMdown " << h1MMdown << endl ;       
    cout << " h1MM " << h1MM << endl ;               
    cout << " h2MM " << h2MM << endl ;               
    cout << " h1GGup " << h1GGup << endl ;
    cout << " h1GGdown " << h1GGdown << endl ;       
    cout << " h1GG " << h1GG << endl ;               
    cout << " h2GG " << h2GG << endl ;               

    gStyle->SetOptStat(0);
   gStyle->SetErrorX(0.5);
   PlotWithRatioEEMuMuGJetsComparison(c1, h1EE, h2EE, h1MM, h2MM, h1GG,h2GG, h1EEup, h1EEdown, h1MMup, h1MMdown, h1GGup, h1GGdown, EjeX, histograma);
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


