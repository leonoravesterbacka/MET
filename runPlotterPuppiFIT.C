#include "TString.h"
#include "TSystem.h"
#include "plotterPuppiFIT.C"
void runPlotWithRatioJES(TString channel){
gSystem->CompileMacro("plotterPuppiFIT.C", channel);
     Comparison("met_uParaPuppi__vs_nVert", "met_uParaPuppi__vs_nVert", "Number of vertices", "Puppi_Type1FIT");       
     Comparison("met_uPerpPuppi__vs_nVert", "met_uPerpPuppi__vs_nVert", "Number of vertices", "Puppi_Type1FIT");       

}

