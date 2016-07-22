#include "TString.h"
#include "TSystem.h"
#include "plotter.C"
void runPlotWithRatioJES(TString channel){
gSystem->CompileMacro("plotter.C", channel);


    //testRatioEEMuMuGJetsComparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]");
    //testRatioEEMuMuGJetsComparison("met_uPara_over_qt", "met_uPara_over_qt", "q_{T} [GeV]");
    //testRatioEEMuMuGJetsComparison("met_metx_over_nVert", "met_metx_over_nVert", "# Vertices");
    testRatioEEMuMuGJetsComparison("met_uPara_vs_zll_pt","gamma_met_uPara_vs_gamma_pt",  "q_{T} [GeV]") ;
    testRatioEEMuMuGJetsComparison("met_uPerp_vs_zll_pt","gamma_met_uPerp_vs_gamma_pt",  "q_{T} [GeV]");   
    //testRatioEEMuMuGJetsComparison("met_metx_vs_nVert","gamma_met_metx_vs_nVert", "# Vertices ");   
    //testRatioEEMuMuGJetsComparison("met_mety_vs_nVert","gamma_met_mety_vs_nVert", "# Vertices ");   
    //testRatioEEMuMuGJetsComparison("met_uParaPuppi_vs_zll_pt","gamma_met_uPara_vs_gamma_pt",  "Boson p_{T} [GeV]") ;
    //testRatioEEMuMuGJetsComparison("met_uPerpPuppi_vs_zll_pt","gamma_met_uPerp_vs_gamma_pt", "Boson p_{T} [GeV]");   
    //testRatioEEMuMuGJetsComparison("met_uPerp_vs_nVert", "met_uPerpPuppi_vs_nVert", "# Vertices");
    //testRatioEEMuMuGJetsComparison("met_uPara_vs_nVert", "met_uParaPuppi_vs_nVert", "# Vertices");  


}

