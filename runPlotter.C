#include "TString.h"
#include "TSystem.h"
#include "plotter.C"
void runPlotWithRatioJES(TString channel){
gSystem->CompileMacro("plotter.C", channel);
    // this one takes the tgraphs "met_uPara_vs_zll_pt" or which ever, and plots it. The first argument is for the zll case, and the second for the photon, since the names of the tgraphs end up being different. this makes the scale or resolution plot comparing the three bosons. But make the resolution and scale separately, because there is an option in the plotter.C that needs to be changed (the binning) between the two different plots. Otherwise it does everything automatically. The file where the tgraphs are is specified in the plotter.C. There is a similar on of these codes for plotting the Puppi plots, since for that, we want just dy mumu pf met vs puppi met.  
   //scale 
//     Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1Mean");
//     Comparison("met_uParaPuppi_over_zll_pt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Puppi_Type1Mean");   
//     Comparison("met_uPara__vs_nVert", "M_met_uParaPuppi__vs_nVert", "N_{vtx}", "Type1_PFvsPuppiMM");
//     Comparison("met_uPerp__vs_nVert", "M_met_uPerpPuppi__vs_nVert", "N_{vtx}", "Type1_PFvsPuppiMM");   
//     Comparison("met_uPara__vs_nVert", "E_met_uParaPuppi__vs_nVert", "N_{vtx}", "Type1_PFvsPuppiEE");   
//     Comparison("met_uPerp__vs_nVert", "E_met_uPerpPuppi__vs_nVert", "N_{vtx}", "Type1_PFvsPuppiEE");   
//     Comparison("met_uParaPuppi__vs_nVert", "met_uParaPuppi__vs_nVert", "N_{vtx}", "Puppi_Type1");
//     Comparison("met_uPerpPuppi__vs_nVert", "met_uPerpPuppi__vs_nVert", "N_{vtx}", "Puppi_Type1");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "N_{vtx}", "PF_Type1nVert");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "N_{vtx}", "PF_Type1nVert");
//     Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qT");   
//     Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qT");   
     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEt");       
     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEt");       
//Supplementary Material
//     Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1FitVsMean");
//     Comparison("met_uPerp__vs_nVert", "met_uPerp__vs_nVert", "N_{vtx}", "PF_Type1nVertFitVsRMS");
//     Comparison("met_uPara__vs_nVert", "met_uPara__vs_nVert", "N_{vtx}", "PF_Type1nVertFitVsRMS");
   
   
   
   
   
   //   Comparison("met_uPara_over_zll_pt", "met_uPara_over_zll_pt", "q_{T} [GeV]", "PF_Type1MatchEle");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_zll_pt", "q_{T} [GeV]", "PF_Type1JetEF");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_zll_pt", "q_{T} [GeV]", "PF_Type1_0jet");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_zll_pt", "q_{T} [GeV]", "PF_Type1_1jet");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1NewBinning");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_zll_pt", "q_{T} [GeV]", "PF_Type1IdTest");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1Mean");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1FitVsMean");
//   Comparison("met_uPara_over_nVert", "met_uPara_over_nVert", "Number of vertices ", "PF_Type1_Scale_nVert");
//   Comparison("met_uPara_over_met_sumEt-zll_pt", "met_uPara_over_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV] ", "PF_Type1_Scale_sumEt");
//   Comparison("met_uPara_over_zll_pt", "met_uPara_over_qt", "q_{T} [GeV]", "PF_Type1MC");
//`   Comparison("met_uParaRaw_over_zll_pt", "met_uParaRaw_over_qt", "q_{T} [GeV]", "PF_Raw");
//   Comparison("met_uParaPuppi_over_zll_pt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Puppi_Type1");
//   Comparison("met_uParaPuppi_over_zll_pt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Puppi_Type1NewBinning");
//   Comparison("met_uParaPuppi_over_zll_pt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Puppi_Type1MC");
//   Comparison("met_uParaPuppi_over_zll_pt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Puppi_Type1Mean");
//   Comparison("met_uParaPuppi_over_nVert", "met_uParaPuppi_over_nVert", "Number of vertices ", "Puppi_Type1_Scale_nVert");
//   Comparison("met_uParaTk_over_zll_pt", "met_uParaTkCHS_over_zll_pt", "q_{T} [GeV]", "tk");
//   Comparison("met_uParaPuppi_over_qt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Type1Puppi_OldvsNew");
//   Comparison("met_uPara_over_zll_pt", "E_met_uParaPuppi_over_zll_pt", "q_{T} [GeV]", "Type1_PFvsPuppiEE");
//   Comparison("met_uPara_over_zll_pt", "M_met_uParaPuppi_over_zll_pt", "q_{T} [GeV]", "Type1_PFvsPuppiMM");
//   Comparison("met_uParaRaw_over_zll_pt", "E_met_uParaPuppiRaw_over_zll_pt", "q_{T} [GeV]", "Raw_PFvsPuppiEE");
//   Comparison("met_uParaRaw_over_zll_pt", "M_met_uParaPuppiRaw_over_zll_pt", "q_{T} [GeV]", "Raw_PFvsPuppiMM");
    //Comparison("E_met_uParaPuppi_over_zll_pt", "E_met_uParaPuppi_over_zll_pt", "q_{T} [GeV]", "Type1_PuppiJetCorr");
   // Comparison("E_met_uParaPuppiRaw_over_zll_pt", "E_met_uParaPuppiRaw_over_zll_pt", "q_{T} [GeV]", "Raw_PuppiJetCorr");
   // Comparison("met_uPara_over_qt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Type1_PuppiGJetsFinal");
   // Comparison("met_uParaRaw_over_qt", "met_uParaPuppiRaw_over_qt", "q_{T} [GeV]", "Raw_PuppiGJetsFinal");
    //Comparison("met_uPara_over_qt", "met_uParaPuppi_over_qt", "q_{T} [GeV]", "Type1_PFvsPuppiGJets");
    //Comparison("met_uParaRaw_over_qt", "met_uParaPuppiRaw_over_qt", "q_{T} [GeV]", "Raw_PFvsPuppiGJets");
    //Comparison("met_uParaRaw_over_qt", "met_uParaPuppiRaw_over_qt", "q_{T} [GeV]", "Raw_PFvsPuppiGJets");
  // //Comparison("met_uParaPuppi_over_zll_pt", "met_uParaPuppiRaw_over_zll_pt", "q_{T} [GeV]", "Puppi_RawvsType1");
   //Comparison("met_uPara_over_zll_pt", "met_uParaRaw_over_zll_pt", "q_{T} [GeV]", "PF_RawvsType1");
   //Comparison("met_uPara_over_qt", "met_uPara_over_qt", "q_{T} [GeV]", "gammaType1");                              
   //resolution
//nVert plots    
//     Comparison("met_uPerp__vs_nVert", "met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertFitVsRMS");
//     Comparison("met_uPara__vs_nVert", "met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertFitVsRMS");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVert");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVert");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertMC");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertMC");
//     Comparison("met_uParaPuppi__vs_nVert", "met_uParaPuppi__vs_nVert", "Number of vertices", "Puppi_Type1");
//     Comparison("met_uPerpPuppi__vs_nVert", "met_uPerpPuppi__vs_nVert", "Number of vertices", "Puppi_Type1");
//     Comparison("met_uParaRaw_vs_nVert", "gamma_met_uParaRaw__vs_nVert", "Number of vertices", "PF_Raw_Res");
//     Comparison("met_uPerpRaw_vs_nVert", "gamma_met_uPerpRaw__vs_nVert", "Number of vertices", "PF_Raw_Res");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertCentral");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertCentral");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertNoSC");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertNoSC");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertZeta");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertZeta");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertnJet50");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertnJet50");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertnJet60");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertnJet60");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertnJet70");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertnJet70");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertZetanJet50");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertZetanJet50");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertNoSC");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertNoSC");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertbkgSub");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVert");
//     Comparison("met_uPerpRaw_vs_nVert", "gamma_met_uPerpRaw__vs_nVert", "Number of vertices", "PF_Type1nVertRaw");
//     Comparison("met_uParaRaw_vs_nVert", "gamma_met_uParaRaw__vs_nVert", "Number of vertices", "PF_Type1nVertRaw");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertAlpha");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertAlpha");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertJetEta");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertJetEta");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertNoJetEta");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertNoJetEta");
//     Comparison("met_uPara__vs_nVert", "gamma_met_uPara__vs_nVert", "Number of vertices", "PF_Type1nVertRunBCD");
//     Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp__vs_nVert", "Number of vertices", "PF_Type1nVertRunBCD");
//qT plots    
//     Comparison("met_uPara__vs_zll_pt", "met_uPara__vs_zll_pt", "q_{T} [GeV]", "PF_Type1qTFitVsRMS");       
//     Comparison("met_uPerp__vs_zll_pt", "met_uPerp__vs_zll_pt", "q_{T} [GeV]", "PF_Type1qTFitVsRMS");       
//     Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTRMS");       
//     Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTRMS");       
//     Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTAlpha");       
//     Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTAlpha");       
//     Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qT");       
//     Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qT");       
//     Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTNoSC");       
//     Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTNoSC");       
//     Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTMC");       
//     Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTMC");       
//     Comparison("met_uPerpPuppi__vs_zll_pt", "gamma_met_uPerpPuppi__vs_gamma_pt", "q_{T} [GeV]", "Puppi_Type1");       
//     Comparison("met_uParaPuppi__vs_zll_pt", "gamma_met_uParaPuppi__vs_gamma_pt", "q_{T} [GeV]", "Puppi_Type1");       
//     Comparison("met_uParaRaw_vs_zll_pt", "gamma_met_uParaRaw__vs_gamma_pt", "q_{T} [GeV]", "PF_Raw_Res");       
//     Comparison("met_uPerpRaw_vs_zll_pt", "gamma_met_uPerpRaw__vs_gamma_pt", "q_{T} [GeV]", "PF_Raw_Res");       
//       Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTCentral");       
//       Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp__vs_gamma_pt", "q_{T} [GeV]", "PF_Type1qTCentral");       
//sumEt plots
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtRMS");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtRMS");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtAlpha");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtAlpha");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEt");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEt");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtMC");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtMC");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1MC");       
//     Comparison("met_uPara__vs_sumEPara", "met_uPara__vs_sumEPara", "#Sigma E_{||} [GeV]", "PF_Type1SumEPara");       
//     Comparison("met_uPerp__vs_sumEPara", "met_uPerp__vs_sumEPara", "#Sigma E_{||} [GeV]", "PF_Type1SumEPara");       
//     Comparison("met_uPara__vs_sumEPerp", "met_uPara__vs_sumEPerp", "#Sigma E_{#perp} [GeV]", "PF_Type1SumEPerp");       
//     Comparison("met_uPerp__vs_sumEPerp", "met_uPerp__vs_sumEPerp", "#Sigma E_{#perp} [GeV]", "PF_Type1SumEPerp");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1MC");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1");       
//     Comparison("met_uParaPuppi__vs_met_sumEt-zll_pt", "gamma_met_uParaPuppi__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "Puppi_Type1");       
//     Comparison("met_uPerpPuppi__vs_met_sumEt-zll_pt", "gamma_met_uPerpPuppi__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "Puppi_Type1");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtCentral");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtCentral");       
//     Comparison("met_uPara__vs_met_sumEt-zll_pt", "gamma_met_uPara__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtNoSC");       
//     Comparison("met_uPerp__vs_met_sumEt-zll_pt", "gamma_met_uPerp__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Type1SumEtNoSC");       
   //  Comparison("met_uParaRaw_vs_met_sumEt-zll_pt", "gamma_met_uParaRaw__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Raw_Res");       
   //  Comparison("met_uPerpRaw_vs_met_sumEt-zll_pt", "gamma_met_uPerpRaw__vs_met_sumEt-gamma_pt", "#Sigma E_{T} [GeV]", "PF_Raw_Res");       
     
     // //Comparison("met_uParaTk__vs_nVert", "met_uParaTkCHS__vs_nVert", "Number of vertices", "tkRes");
    // //Comparison("met_uPerpTk__vs_nVert", "met_uPerpTkCHS__vs_nVert", "Number of vertices", "tkRes");
    // Comparison("met_uPerp__vs_nVert", "gamma_met_uPerp_vs_nVert_scaleCorr", "Number of vertices", "PF_Type1");
    // Comparison("gamma_met_uPerp_vs_nVert_scaleCorr", "gamma_met_uPerp_vs_nVert_scaleCorr", "Number of vertices", "gammaVert");
    // Comparison("gamma_met_uPara_vs_nVert_scaleCorr", "gamma_met_uPara_vs_nVert_scaleCorr", "Number of vertices", "gammaVert");
    // //Comparison("met_uPerp__vs_zll_pt", "gamma_met_uPerp_vs_gamma_pt_scaleCorr", "q_{T} [GeV]", "gammaqt");
    // //Comparison("met_uPara__vs_zll_pt", "gamma_met_uPara_vs_gamma_pt_scaleCorr", "q_{T} [GeV]", "gammaqt");
//     Comparison("met_uPara__vs_nVert", "M_met_uParaPuppi__vs_nVert", "Number of vertices", "Type1_PFvsPuppiMM");
//     Comparison("met_uPerp__vs_nVert", "M_met_uPerpPuppi__vs_nVert", "Number of vertices", "Type1_PFvsPuppiMM");
//     Comparison("met_uPara__vs_nVert", "E_met_uParaPuppi__vs_nVert", "Number of vertices", "Type1_PFvsPuppiEE");
//     Comparison("met_uPerp__vs_nVert", "E_met_uPerpPuppi__vs_nVert", "Number of vertices", "Type1_PFvsPuppiEE");
    // Comparison("gamma_met_uPara_vs_nVert_scaleCorr", "gamma_met_uPara_vs_nVert", "Number of vertices", "gammaType1");
    // Comparison("gamma_met_uPerp_vs_nVert_scaleCorr", "gamma_met_uPerp_vs_nVert", "Number of vertices", "gammaType1");                







}

