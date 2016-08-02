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
TGraphErrors *h1upUncl = new TGraphErrors;
TGraphErrors *h1downUncl = new TGraphErrors;
TGraphErrors *h1puppi = new TGraphErrors;
TGraphErrors *h2puppi = new TGraphErrors;
void testRatioEEMuMuGJetsComparison(TString, TString, TString, TString);                                            
TGraphAsymmErrors staterror(TH1F *, TH1F * );
TGraphAsymmErrors JESerror(TH1F *, TH1F *, TH1F *, TH1F * );
TGraphAsymmErrors TOTerror(TH1F *, TH1F *, TH1F *, TH1F *, TH1F *, TH1F *  );

void PlotWithRatioEEMuMuGJetsComparison(TCanvas *c1, TGraphErrors *mcdyee,  TGraphErrors *datadyee, TGraphErrors *mcdypuppi,  TGraphErrors *datadypuppi, TGraphErrors *mcdyeeup, TGraphErrors *mcdyeedown, TGraphErrors *mcdyeeupUncl, TGraphErrors *mcdyeedownUncl,  TString EjeX, TString histogramaname, TString channel) {

TH1::SetDefaultSumw2();
  c1->cd();
 //   Float_t xbins[] = { 0, 15, 30, 51, 65, 80, 110, 140,170, 200,250, 330};
 //    double nbins=11.;
 //     cout << "plotting resolution as a function of qt" << endl;
   // if (histogramaname.Contains("over")){
    Float_t xbins[] = { 20,28,36, 44,52,60, 68,76, 84,92,100, 120, 150, 175, 200, 225,250, 275, 305, 335, 365, 395, 430, 500}; 
     double nbins=23.;
 //     cout << "plotting scale " << endl;
 //   }
//    else if (histogramaname.Contains("nVert")){
  //  Float_t xbins[] = { 0, 6, 8, 10, 12 ,14, 16, 20,  27};
  //    double nbins=8.;
  //    cout << "plotting resolution as a function of nVert" << endl;
  //  }                                                                                                                                         



   TH1F * hdatadyee= new TH1F("","", nbins,xbins);
   TH1F * hdatadypuppi= new TH1F("hdatadypuppi","hdatadypuppi", nbins,xbins);
   TH1F * hmcdyee= new TH1F("hmcdyee","hmcdyee", nbins,xbins);
   TH1F * hmcdypuppi= new TH1F("hmcdypuppi","hmcdypuppi", nbins,xbins);
   TH1F * hmcdyeeup= new TH1F("hmcdyeeup","hmcdyeeup", nbins, xbins);
   TH1F * hmcdyeedown= new TH1F("hmcdyeedown","hmcdyeedown", nbins,xbins);
   TH1F * hmcdyeeupUncl= new TH1F("hmcdyeeupUncl","hmcdyeeupUncl", nbins, xbins);
   TH1F * hmcdyeedownUncl= new TH1F("hmcdyeedownUncl","hmcdyeedownUncl", nbins, xbins);
   
   for (int i = 0, n = nbins; i < n; ++i) {
       hdatadyee->SetBinContent(i+1,datadyee->GetY()[i]);
       hdatadypuppi->SetBinContent(i+1,datadypuppi->GetY()[i]);
       hmcdyee->SetBinContent(i+1,mcdyee->GetY()[i]);
       hmcdypuppi->SetBinContent(i+1,mcdypuppi->GetY()[i]);
       hmcdyeedown->SetBinContent(i+1,mcdyeedown->GetY()[i]);
       hmcdyeeup->SetBinContent(i+1,mcdyeeup->GetY()[i]);
       hmcdyeedownUncl->SetBinContent(i+1,mcdyeedownUncl->GetY()[i]);
       hmcdyeeupUncl->SetBinContent(i+1,mcdyeeupUncl->GetY()[i]);
 	   hdatadyee->SetBinError(i+1,datadyee->GetErrorY(i));
       hmcdyee->SetBinError(i+1,mcdyee->GetErrorY(i));          
       hmcdypuppi->SetBinError(i+1,mcdypuppi->GetErrorY(i));                  
       hdatadypuppi->SetBinError(i+1,datadypuppi->GetErrorY(i));                  
    }

    hdatadyee->SetMarkerStyle(23);
    hdatadyee->SetMarkerColor(4);
    hdatadyee->SetLineColor(4);   
    hdatadypuppi->SetMarkerStyle(22);
    hdatadypuppi->SetMarkerColor(46);
    hdatadypuppi->SetLineColor(46);   
    hmcdyee->SetMarkerColor(4);
    hmcdyee->SetMarkerStyle(22); 
    TH1F *ratiodyee=(TH1F*)hdatadyee->Clone(); 
    TH1F *ratiodypuppi=(TH1F*)hdatadypuppi->Clone(); 
   
    ratiodyee->Divide(hmcdyee);
    ratiodypuppi->Divide(hmcdypuppi);
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.014);
    pad1->SetLeftMargin(0.1);
    pad1->Draw();
    pad1->cd();
    hdatadyee->GetYaxis()->SetRangeUser(0, 40);
    if (histogramaname.Contains("over")) hdatadyee->GetYaxis()->SetRangeUser(0.8,1.2);
    hdatadyee->GetXaxis()->SetLabelFont(43); //font in pixels
    hdatadyee->GetXaxis()->SetLabelSize(1); //in pixels
    hdatadyee->GetYaxis()->SetLabelFont(43);
    hdatadyee->GetYaxis()->SetLabelSize(16);
    hdatadyee->GetYaxis()->SetTitleFont(43); //font in pixels
    hdatadyee->GetYaxis()->SetTitleSize(20); //in pixels
    hdatadyee->GetXaxis()->SetTitleOffset(3.3);
    //hdatadyee->Draw("e1");
    pad1->Update();
    TString ylabel="";
if (histogramaname=="met_uParavs_nVert") ylabel="#sigma ( u_{||} ) [GeV]";
if (histogramaname=="met_uPerpvs_nVert") ylabel="#sigma ( u_{#perp}  ) [GeV]";
if (histogramaname=="met_uPara_over_zll_pt") ylabel="   -<u_{||}/q_{T}> ";
if (histogramaname.Contains("over_zll")) hdatadyee->GetYaxis()->SetRangeUser(0.58,1.2);  
if (histogramaname.Contains("over_nVert")) hdatadyee->GetYaxis()->SetRangeUser(-5,5);  
   
    hdatadyee->Draw("e1, X0");
    //hdatadyee->Draw("AXIS");
    hdatadypuppi->Draw("e1, X0 same");
    if (histogramaname.Contains("over")){   
        TLine *line =  new TLine ( hdatadyee->GetXaxis ()->GetXmin (), 1, hdatadyee->GetXaxis ()->GetXmax (), 1);
        line->SetLineColor (kBlack); line -> SetLineStyle(2); line->Draw();
    }                                                                                                              
    
    TLegend* leg(0);
    leg = new TLegend(0.15,0.65,0.35,0.85);
    if (histogramaname.Contains("over")){
        leg = new TLegend(0.45,0.15,0.8,0.35);
    }
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillColor(0);
    if (channel == "EE"){
        leg->AddEntry(hdatadyee, "PF E_{T}^{miss}     Z#rightarrow ee", "lp");
        leg->AddEntry(hdatadypuppi, "Puppi E_{T}^{miss} Z#rightarrow ee  ", "lp");
    }
    else{
        leg->AddEntry(hdatadyee, "PF E_{T}^{miss}     Z#rightarrow #mu#mu", "lp");
        leg->AddEntry(hdatadypuppi, "Puppi E_{T}^{miss} Z#rightarrow #mu#mu  ", "lp");
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
    ratiodyee->GetYaxis()->SetRangeUser(0.7, 1.3);
    if (histogramaname.Contains("over")) 
        ratiodyee->GetYaxis()->SetRangeUser(0.9,1.1);
    ratiodyee->SetMarkerStyle(23);
    ratiodyee->SetMarkerColor(4);
    ratiodyee->SetLineColor(4);
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
   toterree.SetFillColor (kBlue);
   toterree.SetFillStyle (3002);                                   
   
   ratiodyee->Draw("e1, X0"); //"e2"
   //ratiodygg->Draw("e1, X0 same "); //"e2"
   ratiodypuppi->Draw("e1, X0 same"); //"e2"
   ratiodyee->Draw("e1, X0 same"); //"e2"
   toterree.Draw ("2 same");
    erree.Draw ("2 same");
   staterree.Draw("2 same");
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
    legratio->Draw();
    //ratiodygg->Draw("e1, X0 same "); //"e2"
    ratiodypuppi->Draw("e1, X0 same"); //"e2"
    ratiodyee->Draw("e1, X0 same"); //"e2"

    pad1->cd();
    TLatex latex;                  
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.07);
    latex.DrawLatex(0.21, 0.91, "#bf{CMS}");

    if ((histogramaname.Contains("sumEt"))  or (histogramaname.Contains("nVert")) ) {
        TLatex latexz;
        latexz.SetNDC();
        latexz.SetTextAngle(0);
        latexz.SetTextFont(42);
        latexz.SetTextAlign(31);
        latexz.SetTextSize(0.05);
        latexz.DrawLatex(0.36, 0.58, "q_{T} > 50 GeV");
        latexz.DrawLatex(0.8, 0.13, "Response Corrected");
    }
    
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
    c1->Print("~/www/met/comparisons/Aug02/"+histogramaname+channel+".png");
    c1->Print("~/www/met/comparisons/Aug02/"+histogramaname+channel+".root");
    c1->Print("~/www/met/comparisons/Aug02/"+histogramaname+channel+".pdf");
    c1->Print("~/www/met/comparisons/Aug02/"+histogramaname+channel+".C");

}


void testRatioEEMuMuGJetsComparison(TString histograma, TString histograma2, TString EjeX, TString channel) {                                               

TCanvas *c1 = new TCanvas("c1","example",600,700);
cout << "histograma"<< histograma << endl;
if (histograma.Contains("over")){
if (channel == "EE"){
    TFile file1MC ("DYEZllScaleV3.root");
    TFile file2Data ("DataEZllScaleV3.root");
    TFile file1up ("DY_up_jes_DYEZllScaleV3.root"); 
    TFile file1down ("DY_down_jes_DYEZllScaleV3.root"); 
    TFile file1upUncl ("DY_up_uncl_DYEZllScaleV3.root");
    TFile file1downUncl ("DY_down_uncl_DYEZllScaleV3.root");

    h1 = (TGraphErrors*) file1MC.Get("E_"+ histograma); 
    h2 = (TGraphErrors*) file2Data.Get("E_"+ histograma); 
    h1puppi = (TGraphErrors*) file1MC.Get("E_"+histograma2); 
    h2puppi = (TGraphErrors*) file2Data.Get("E_"+histograma2); 
    h1up = (TGraphErrors*)   file1up.Get("E_"+ histograma); 
    h1down = (TGraphErrors*) file1down.Get("E_"+ histograma); 
    h1upUncl = (TGraphErrors*)   file1upUncl.Get("E_"+ histograma); 
    h1downUncl = (TGraphErrors*)   file1downUncl.Get("E_"+ histograma); 

}
else if (channel== "MM"){                                                                                                                       
     TFile file1MC ("DYMZllScaleTestwoOF.root");
     TFile file2Data ("DataMZllScaleTestwoOF.root");
     TFile file1PUPPIMC ("DYMZllScaleTestwoOF.root");
     TFile file2PUPPIData ("DataMZllScaleTestwoOF.root");
     TFile file1up ("DY_up_jes_DYMZllScaleTestwoOF.root"); 
     TFile file1down ("DY_down_jes_DYMZllScaleTestwoOF.root"); 
     TFile file1upUncl ("DY_up_uncl_DYMZllScaleTestwoOF.root");
     TFile file1downUncl ("DY_down_uncl_DYMZllScaleTestwoOF.root");
    
     h1 = (TGraphErrors*) file1MC.Get("M_"+ histograma); 
     h2 = (TGraphErrors*) file2Data.Get("M_"+ histograma); 
     h1puppi = (TGraphErrors*) file1PUPPIMC.Get("M_"+histograma2); 
     h2puppi = (TGraphErrors*) file2PUPPIData.Get("M_"+histograma2); 
     h1up = (TGraphErrors*)   file1up.Get("M_"+ histograma); 
     h1down = (TGraphErrors*) file1down.Get("M_"+ histograma); 
     h1upUncl = (TGraphErrors*)   file1upUncl.Get("M_"+ histograma); 
     h1downUncl = (TGraphErrors*)   file1downUncl.Get("M_"+ histograma); 
 }                                               
}                                                                                   
else{

if (channel == "EE"){
    TFile file1MC ("DYEZllResolutionV3.root");
    TFile file2Data ("DataEZllResolutionV3.root");
    TFile file1up ("DY_up_uncl_DYEZllResolutionV3.root"); 
    TFile file1down ("DY_down_jes_DYEZllResolutionV3.root"); 
    TFile file1upUncl ("DY_up_uncl_DYEZllResolutionV3.root");
    TFile file1downUncl ("DY_down_uncl_DYEZllResolutionV3.root");

    h1 = (TGraphErrors*) file1MC.Get("E_"+ histograma); 
    h2 = (TGraphErrors*) file2Data.Get("E_"+ histograma); 
    h1puppi = (TGraphErrors*) file1MC.Get("E_"+histograma2); 
    h2puppi = (TGraphErrors*) file2Data.Get("E_"+histograma2); 
    h1up = (TGraphErrors*)   file1up.Get("E_"+ histograma); 
    h1down = (TGraphErrors*) file1down.Get("E_"+ histograma); 
    h1upUncl = (TGraphErrors*)   file1upUncl.Get("E_"+ histograma); 
    h1downUncl = (TGraphErrors*)   file1downUncl.Get("E_"+ histograma); 

}
else if (channel== "MM"){                                                               
     TFile file1MC ("DYMZllResolutionPuppi.root");
     TFile file2Data ("DataMZllResolutionPuppi.root");
     TFile file1PUPPIMC ("DYMZllResolutionPuppi.root");
     TFile file2PUPPIData ("DataMZllResolutionPuppi.root");
     TFile file1up ("DY_up_jes_DYMZllResolutionPuppi.root"); 
     TFile file1down ("DY_down_jes_DYMZllResolutionPuppi.root"); 
     TFile file1upUncl ("DY_up_uncl_DYMZllResolutionPuppi.root");
     TFile file1downUncl ("DY_down_uncl_DYMZllResolutionPuppi.root");
     h1 = (TGraphErrors*) file1MC.Get("M_"+ histograma); 
     h2 = (TGraphErrors*) file2Data.Get("M_"+ histograma); 
     h1puppi = (TGraphErrors*) file1PUPPIMC.Get("M_"+histograma2); 
     h2puppi = (TGraphErrors*) file2PUPPIData.Get("M_"+histograma2); 
     h1up = (TGraphErrors*)   file1up.Get("M_"+ histograma); 
     h1down = (TGraphErrors*) file1down.Get("M_"+ histograma); 
     h1upUncl = (TGraphErrors*)   file1upUncl.Get("M_"+ histograma); 
     h1downUncl = (TGraphErrors*)   file1downUncl.Get("M_"+ histograma); 
 }                                               
}                                                                                   

    cout << " up " << h1up << endl ;
    cout << " down" << h1down << endl ;
    cout << " unclUp " << h1upUncl << endl ;
    cout << " unclDown " << h1downUncl << endl ;
    cout << " h1 " << h1 << endl ;               
    cout << " h2 " << h2 << endl ;               
    cout << " h1puppi " << h1puppi << endl ;               
    cout << " h2puppi " << h2puppi << endl ;               

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.5);
    //PlotWithRatioEEMuMuGJetsComparison(c1, h1, h1, h2, h2, h1up, h1down, h1upUncl, h1downUncl, EjeX, histograma, channel);
    PlotWithRatioEEMuMuGJetsComparison(c1, h1, h2, h1puppi, h2puppi, h1up, h1down, h1upUncl, h1downUncl, EjeX, histograma, channel);
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
                                                                                                                                                                                                                                                                                                                                                                                                 





