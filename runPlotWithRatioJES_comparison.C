#include "TString.h"
#include "TSystem.h"
#include "PlotWithRatioJES_comparisonEEMUMU.C"
void runPlotWithRatioJES(TString channel){
gSystem->CompileMacro("PlotWithRatioJES_comparisonEEMUMU.C", channel);


    testRatioEEMuMuGJetsComparison("met_uPara_over_qt", "met_uPara_over_qt", "q_{T} [GeV]");
     // testRatioEEMuMuComparison("met_uPara_over_qt", , " p_{T} [GeV]");
   // testRatioEEMuMuGJetsComparison("met_uPara_vs_zll_pt","gamma_met_uPara_vs_gamma_pt",  "Boson p_{T} [GeV]") ;
   // testRatioEEMuMuGJetsComparison("met_uPerp_vs_zll_pt","gamma_met_uPerp_vs_gamma_pt", "Boson p_{T} [GeV]");
   // testRatioEEMuMuGJetsComparison("met_uPerp_vs_nVert", "gamma_met_uPerp_vs_nVert", "# Vertices");
   // testRatioEEMuMuGJetsComparison("met_uPara_vs_nVert", "gamma_met_uPara_vs_nVert", "# Vertices");
   // testRatioEEMuMuGJetsComparison("met_uPara_vs_met_sumEt-zll_pt", "gamma_met_uPara_vs_met_sumEt-gamma_pt",  "#Sigma E_{T} [GeV]");
   // testRatioEEMuMuGJetsComparison("met_uPerp_vs_met_sumEt-zll_pt", "gamma_met_uPerp_vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]");              
//


}

