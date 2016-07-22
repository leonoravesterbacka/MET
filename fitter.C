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
TGraphErrors *hDataX = new TGraphErrors;
TGraphErrors *hDataY = new TGraphErrors;
TGraphErrors *hMCX = new TGraphErrors;
TGraphErrors *hMCY = new TGraphErrors;
void testRatioEEMuMuGJetsComparison(TString, TString, TString);                                            


void fitter(TString histogramX, TString histogramY, TString EjeX) {                                               
    
  //  TCanvas *cDataX = new TCanvas("cdatax","example",600,700);
  //  TCanvas *cDataY = new TCanvas("cdatay","example",600,700);
  //  TCanvas *cMCX = new TCanvas("cmcx","example",600,700);
   TCanvas *cMCY = new TCanvas("cmcy","example",600,700);
    TFile fileData ("DataMZllScaleXY.root");                      
    TFile fileMC ("DYMZllScaleXY.root");                      
    
  //  hDataX = (TGraphErrors*) fileData.Get("M_"+ histogramX); 
  //  hDataY = (TGraphErrors*) fileData.Get("M_"+ histogramY); 
  //  hMCX = (TGraphErrors*) fileMC.Get("M_"+ histogramX); 
    hMCY = (TGraphErrors*) fileMC.Get("M_"+ histogramY); 

 //   hDataX->Draw("a*");
 //   hDataX->Fit("pol1");
 //   TF1 *fDataX = hDataX->GetFunction("pol1");

 //   hDataY->Draw("a*");
 //   hDataY->Fit("pol1");
 //   TF1 *fDataY = hDataY->GetFunction("pol1");

 //   hMCX->Draw("a*");
 //   hMCX->Fit("pol1");
 //   TF1 *fMCX = hMCX->GetFunction("pol1");
    //fMCX->SetLineWidth(1);
  //                                                     
   hMCY->Draw("a*");
   hMCY->Fit("pol1");
   TF1 *fMCY = hMCY->GetFunction("pol1");
 //   cDataX->Print("~/www/met/fits/DataMetXFits.png");  
    //cDataY->Print("~/www/met/fits/DataMetYFits.png");   
    cMCY->Print("~/www/met/fits/MCMetYFits.png");  
   // cMCX->Print("~/www/met/fits/MCMetXFits.png");  
}
