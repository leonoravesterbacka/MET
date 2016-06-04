#include "TString.h"
#include "TSystem.h"
#include "PlotWithRatioJES.C"
void runPlotWithRatioJES(TString channel){
gSystem->CompileMacro("PlotWithRatioJES.C", channel);

  //  if (channel == "GG"){
  //      testRatio("met_uPara_over_qt", "#gamma p_{T} [GeV]", channel);
  //      testRatio("gamma_met_uPara_vs_gamma_pt", "#gamma p_{T} [GeV]",channel) ;
  //      testRatio("gamma_met_uPerp_vs_gamma_pt","#gamma p_{T} [GeV]", channel);                
  //      testRatio("gamma_met_uPara_vs_nVert", " # Vertices",channel) ;
  //      testRatio("gamma_met_uPerp_vs_nVert"," # Vertices", channel);                
  //      testRatio("gamma_met_uPara_vs_met_sumEt", "#Sigma E_{T} [GeV]",channel) ;
  //      testRatio("gamma_met_uPerp_vs_met_sumEt","#Sigma E_{T} [GeV]", channel);                
  //  }
  //  if ((channel == "EE") || channel == "MM")  {
       // testRatio("met_uParaRaw_over_qt", "met_uParaPuppi_over_qt", "Z p_{T} [GeV]", "MM");
       // testRatio("met_uPara_vs_zll_pt","met_uParaPuppi_vs_zll_pt",  "Z p_{T} [GeV]",channel) ;
       // testRatio("met_uPerp_vs_zll_pt","met_uPerpPuppi_vs_zll_pt", "Z p_{T} [GeV]", channel);       
       // testRatio("met_uPerp_vs_nVert","met_uPerpPuppi_vs_nVert", "# Vertices", channel);
       // testRatio("met_uPara_vs_nVert","met_uParaPuppi_vs_nVert", "# Vertices ", channel);    

    if ((channel == "EE") || channel == "MM")  {
      testRatio("met_uPara_over_qt", "met_uPara_over_qt", "Z p_{T} [GeV]", channel);
      testRatio("met_uPara_vs_zll_pt","met_uPara_vs_zll_pt",  "Z p_{T} [GeV]",channel); 
      testRatio("met_uPerp_vs_zll_pt","met_uPerp_vs_zll_pt", "Z p_{T} [GeV]", channel);
      testRatio("met_uPerp_vs_nVert" ,"met_uPerp_vs_nVert", "# Vertices", channel);
      testRatio("met_uPara_vs_nVert" ,"met_uPara_vs_nVert", "# Vertices ", channel);    
      testRatio("met_uPerp_vs_met_sumEt-zll_pt" ,"met_uPerp_vs_met_sumEt-zll_pt", "#Sigma E_{T} [GeV]", channel);
      testRatio("met_uPara_vs_met_sumEt-zll_pt" ,"met_uPara_vs_met_sumEt-zll_pt", "#Sigma E_{T} [GeV]", channel);
    }

   // }

   // testRatioEEMuMuComparison("met_uPara_over_qt", " p_{T} [GeV]");
//    testRatioEEMuMuComparison("met_uPara_vs_zll_pt", "Z p_{T} [GeV]") ;
//    testRatioEEMuMuComparison("met_uPerp_vs_zll_pt","Z p_{T} [GeV]");
//    testRatioEEMuMuComparison("met_uPerp_vs_nVert","# Vertices");
//    testRatioEEMuMuComparison("met_uPara_vs_nVert","# Vertices");
//    testRatioEEMuMuComparison("met_uPara_vs_met_sumEt", "#Sigma E_{T} [GeV]");
//    testRatioEEMuMuComparison("met_uPerp_vs_met_sumEt","#Sigma E_{T} [GeV]");              



}

