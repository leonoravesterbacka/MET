##############################################################
##############################################################
######                                                 #######
######  |##\     /##|      |########| |##############| #######
######  |###\   /###|      |##|             |##|       #######
######  |####\_/####|      |##|             |##|       #######
######  |##|\###/|##|      |######|         |##|       #######
######  |##|     |##|      |##|             |##|       #######
######  |##|     |##|      |##|             |##|       #######
######  |##|     |##|      |##|             |##|       #######
######  |##|     |##|      |########|       |##|       #######
######                                                 #######             
##############################################################
##############################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, TMath,  SetOwnership
import math, sys, optparse, array, time, copy

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as Rounder        

def makeUPara(met, phi,boson, boson_phi):
    uPara = ( "((( -"+ met + "*cos("+phi +" ) "+ " - " +boson + "*cos( " + boson_phi +"))* " + boson + "*cos( " +boson_phi +" )+(- "+ met + "*sin( "+ phi+ " )- " +boson+"*sin(" + boson_phi +" ))* " +boson +"*sin( " +boson_phi +" ))/" +boson+ " + " +boson+")" )
    print uPara
    return uPara
def makeUPerp(met, phi,boson, boson_phi):
    uPerp = ( "((( -"+ met + "*cos("+phi +" ) "+ " - " +boson + "*cos( " + boson_phi +"))* " + boson + "*sin( " +boson_phi +" )-(- "+ met + "*sin( "+ phi+ " )- " +boson+"*sin(" + boson_phi +" ))* " +boson +"*cos( " +boson_phi +" ))/" +boson +")" )                              
    print uPerp
    return uPerp                                                                                                                                                                                                                                                                      
def makeMET(met):
# yes this is the smartest function 
    justMET = met 
    return justMET
if __name__ == "__main__":

    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()
    doDY = False
    doee = False
    doReReco = True
    print 'Going to load DATA and MC trees...'
    if doDY:
        if doee:
            ttDatasets = ['TTJets_DiLepton_total_ee', 'TBarToLep_tch_ee', 'TBar_tWch_ee', 'T_tWch_ee', 'TToLep_sch_ee','TToLep_tch_ee', 'TTJets_SingleLeptonFromT_ee', 'TTJets_SingleLeptonFromTbar_ee']
            dyDatasets = ['DYJetsToLL_M50_ee']
            ewkDatasets = ['WWTo2L2Nu_ee', 'WZTo2L2Q_ee', 'WZTo3LNu_ee', 'ZZTo4L_ee', 'ZZTo2L2Nu_ee', 'ZZTo2L2Q_ee', 'WWW_ee', 'WWZ_ee', 'WZZ_ee', 'ZZZ_ee' ]
            #daDatasets = ['DoubleEG_Run2016B_PromptReco_v2', 'DoubleEG_Run2016C_PromptReco_v2', 'DoubleEG_Run2016D_PromptReco_v2', 'DoubleEG_Run2016E_PromptReco_v2', 'DoubleEG_Run2016F_PromptReco_v2', 'DoubleEG_Run2016G_PromptReco_v2', 'DoubleEG_Run2016H_PromptReco_v2', 'DoubleEG_Run2016H_PromptReco_v3']  
            daDatasets = ['DoubleEG_Run2016B_ReReco', 'DoubleEG_Run2016C_ReReco', 'DoubleEG_Run2016D_ReReco', 'DoubleEG_Run2016E_ReReco', 'DoubleEG_Run2016F_ReReco', 'DoubleEG_Run2016G_ReReco', 'DoubleEG_Run2016H_PromptReco_v2', 'DoubleEG_Run2016H_PromptReco_v3']  
            channel = 'EE'
            lumi = 36.4
            run = 'RunB-H' 
        else:
            ttDatasets = ['TTJets_DiLepton_total_mm',  'TBarToLep_tch_mm', 'TBar_tWch_mm', 'T_tWch_mm', 'TToLep_sch_mm','TToLep_tch_mm', 'TTJets_SingleLeptonFromT_mm', 'TTJets_SingleLeptonFromTbar_mm']
            dyDatasets = ['DYJetsToLL_M50_mm']
            ewkDatasets = ['WWTo2L2Nu_mm', 'WZTo2L2Q_mm', 'WZTo3LNu_mm', 'ZZTo4L_mm', 'ZZTo2L2Nu_mm', 'ZZTo2L2Q_mm',  'WWW_mm', 'WWZ_mm', 'WZZ_mm', 'ZZZ_mm' ]
            daDatasets = ['DoubleMuon_Run2016H_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v3']  
            #daDatasets = ['DoubleMuon_Run2016B_ReReco', 'DoubleMuon_Run2016C_ReReco', 'DoubleMuon_Run2016D_ReReco', 'DoubleMuon_Run2016E_ReReco', 'DoubleMuon_Run2016F_ReReco', 'DoubleMuon_Run2016G_ReReco', 'DoubleMuon_Run2016H_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v3']  
            channel = 'MM'
            lumi = 8.75
            run = 'RunH' 
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0)
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
        treeEWK = Sample.Tree(helper.selectSamples(opts.sampleFile, ewkDatasets, 'EWK'), 'EWK'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
        mcTrees = [  treeTT, treeEWK,  treeDY]   
        boson = 'zll_pt'
        boson_phi = 'zll_phi'
    else:
        qcdDatasets = ['QCD_HT200to300_ext', 'QCD_HT300to500_ext', 'QCD_HT500to700_ext',  'QCD_HT700to1000_ext','QCD_HT1000to1500_ext', 'QCD_HT1500to2000_ext', 'QCD_HT2000toInf_ext']
        gjetsDatasets = ['GJets_HT40to100_ext','GJets_HT100to200_ext', 'GJets_HT200to400_ext', 'GJets_HT400to600_ext','GJets_HT600toInf_ext' ]
        ewkDatasets =  ['TGJets_ext', 'TTGJets', 'WGToLNuG_amcatnlo_ext', 'ZGTo2LG_ext', 'ZGTo2NuG', 'WJetsToLNu_HT100to200_ext', 'WJetsToLNu_HT200to400_ext','WJetsToLNu_HT400to600_ext','WJetsToLNu_HT600to800_ext', 'WJetsToLNu_HT800to1200_ext','WJetsToLNu_HT1200to2500_ext', 'WJetsToLNu_HT2500toInf_ext']
      # for debugging purposes, here one can split up the various samples, just uncomment, and pick some nice colors for the different samples in the samples.dat
        if doReReco:
            #daDatasets = ['SinglePhoton_Run2016B_03Feb2017_v2', 'SinglePhoton_Run2016C_03Feb2017', 'SinglePhoton_Run2016D_03Feb2017', 'SinglePhoton_Run2016E_03Feb2017', 'SinglePhoton_Run2016F_03Feb2017', 'SinglePhoton_Run2016G_03Feb2017']
            daDatasets = ['SinglePhoton_Run2016B_03Feb2017_v2', 'SinglePhoton_Run2016C_03Feb2017', 'SinglePhoton_Run2016D_03Feb2017', 'SinglePhoton_Run2016E_03Feb2017', 'SinglePhoton_Run2016F_03Feb2017', 'SinglePhoton_Run2016G_03Feb2017', 'SinglePhoton_Run2016H_03Feb2017_v2']
            run_str ='ReMINIAOD'
        if doReReco == False:
            daDatasets = ['SinglePhoton_Run2016B_PromptReco_v2','SinglePhoton_Run2016C_PromptReco_v2', 'SinglePhoton_Run2016D_PromptReco_v2', 'SinglePhoton_Run2016E_PromptReco_v2', 'SinglePhoton_Run2016F_PromptReco_v1', 'SinglePhoton_Run2016G_PromptReco_v1', 'SinglePhoton_Run2016H_PromptReco_v2', 'SinglePhoton_Run2016H_PromptReco_v3']
            run_str ='Prompt'
        channel = '' 
        treeQCD =   Sample.Tree(helper.selectSamples(opts.sampleFile, qcdDatasets, 'QCD'), 'QCD', 0)
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'GJETS'),'GJETS', 0)
        treeEWK =   Sample.Tree(helper.selectSamples(opts.sampleFile, ewkDatasets, 'EWK'), 'EWK', 0)
       # treeWG =   Sample.Tree(helper.selectSamples(opts.sampleFile, wgDatasets, 'WG'), 'WG', 0)
       # treeWJ =   Sample.Tree(helper.selectSamples(opts.sampleFile, wjDatasets, 'WJ'), 'WJ', 0)
       # treeZG =   Sample.Tree(helper.selectSamples(opts.sampleFile, zgDatasets, 'ZG'), 'ZG', 1)
       # treeTTG =   Sample.Tree(helper.selectSamples(opts.sampleFile, ttgDatasets, 'TTG'), 'TTG', 0)
        treeDA =    Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
        #mcTrees = [  treeQCD, treeGJETS]  
        #mcTrees = [  treeTTG, treeZG, treeWJ, treeWG, treeQCD, treeGJETS]  
        #mcTrees = [treeQCD, treeGJETS]  
        mcTrees = [  treeEWK, treeQCD, treeGJETS]  
        boson = 'gamma_pt'
        boson_phi = 'gamma_phi'
        #lumi = 0.2
        lumi = 35.9
        
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    cuts = CutManager.CutManager()
    color = helper.color

    lumi_str =channel+ 'lumi'+str(lumi)

    etabins = [-10,-9.5, -9, -8.5, -8,  -7.5, -7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5,1,  1.5, 2, 2.5,3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    hoverebins = [0.01, 0.015, 0.02, 0.025, 0.03,0.035, 0.040, 0.045, 0.05,0.055,   0.06, 0.065, 0.07, 0.075,  0.08,0.085, 0.09, 0.095, 0.1]
    chibins = [ 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0]
    sigmabins = [0., 0.001, 0.002, 0.003, 0.004,0.005,0.006, 0.007, 0.008, 0.009,  0.01 ]
    r9bins = [0.8 , 0.81, 0.82,0.83,  0.84,0.85, 0.86, 0.87,0.88, 0.89, 0.9,0.91, 0.92,0.93, 0.94,0.95, 0.96,0.97, 0.98,0.99, 1.0]
    phibins = [ -3., -2.9,-2.8, -2.7, -2.6,-2.5,  -2.4,-2.3, -2.2,-2.1, -2,-1.9, -1.8,-1.7, -1.6,-1.5, -1.4,-1.3, -1.2,-1.1, -1.,-0.9, -0.8,-0.7, -0.6, -0.5,-0.4,-0.3, -0.2,-0.1, 0.,0.1,  0.2,0.3,  0.4,0.5,  0.6,0.7,  0.8,0.9 ,  1.,1.1,  1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9,  2., 2.1, 2.2,2.3,  2.4,2.5,  2.6, 2.7, 2.8,2.9,  3.]
    regions = []
    met = []
    gammabin = range(50, 300, 2) 
    binsig = range(0, 100, 5) 
    binx = range(0, 400, 5) 
    bin2 = range(-200, 200, 5)
    bin1 = range(0, 200, 5)
    sumbin = range(0, 30000, 1000)
    nbin = range(0, 50, 2)
    fullbin = range(0, 10000, 1)
    doHighMetValues = False
    # doVariables below is a list of all the variables which should be plotted, select some or all and go for a coffee.
    if doDY:
        doVariables = ['met_pt_0jet']
        #doVariables = ['nVert', 'zll_pt', 'met_pt', 'metPuppi_pt', 'met_pt_0jet','met_pt_1jet', 'met_rawPt', 'met_rawPt_0jet',  'met_rawPt_1jet','metPuppi_rawPt',  'metPuppi_rawPt_0jet',  'metPuppi_rawPt_1jet', 'metPuppi_pt',  'metPuppi_pt_0jet',  'metPuppi_pt_1jet']
        #doVariables = [ 'met_x', 'met_x_0jet','met_x_1jet',  'met_y','met_y_0jet', 'met_y_1jet','met_sumEt', 'met_sumEt_0jet','met_sumEt_1jet',  'met_rawSumEt',  'met_rawSumEt_0jet',  'met_rawSumEt_1jet']
        #doVariables = [ 'met_uPerp',  'met_uPara','metRaw_uPerp',  'metRaw_uPara',  'metPuppi_uPerp', 'metPuppi_uPara']
        #doVariables = [ 'met_uPerp', 'met_uPerp_0jet', 'met_uPerp_1jet', 'met_uPara',  'met_uPara_0jet', 'met_uPara_1jet','metRaw_uPerp', 'metRaw_uPerp_0jet','metRaw_uPerp_1jet', 'metRaw_uPara', 'metRaw_uPara_0jet','metRaw_uPara_1jet','metPuppi_uPerp', 'metPuppi_uPerp_0jet','metPuppi_uPerp_1jet','metPuppi_uPara','metPuppi_uPara_0jet','metPuppi_uPara_1jet',  'metPuppiRaw_uPerp', 'metPuppiRaw_uPerp_0jet','metPuppiRaw_uPerp_1jet','metPuppiRaw_uPara','metPuppiRaw_uPara_0jet','metPuppiRaw_uPara_1jet']
        #doVariables = [ 'met_sig', 'met_sig_1jet', 'metPuppi_sig', 'metPuppi_sig_1jet']
        #doVariables = ['nVert', 'zll_pt', 'met_x', 'met_x_0jet','met_x_1jet',  'met_y','met_y_0jet', 'met_y_1jet','met_sumEt', 'met_sumEt_0jet','met_sumEt_1jet',  'met_rawSumEt',  'met_rawSumEt_0jet',  'met_rawSumEt_1jet', 'met_uPerp', 'met_uPerp_0jet', 'met_uPerp_1jet', 'met_uPara',  'met_uPara_0jet', 'met_uPara_1jet', 'met_pt', 'met_pt_0jet','met_pt_1jet', 'met_rawPt', 'met_rawPt_0jet',  'met_rawPt_1jet', 'metPuppi_pt',  'metPuppi_pt_0jet',  'metPuppi_pt_1jet','metRaw_uPerp', 'metRaw_uPara',  'metPuppi_uPerp', 'metPuppi_uPara', 'metPuppi_rawPt', 'metPuppiRaw_uPerp', 'metPuppiRaw_uPara', 'met_sig','met_sig_1jet']
        #binnings    = [nbin, bin1, bin1, bin1, bin1, bin1, bin1, bin1, bin1, bin1, bin1, bin1]
        #binnings    = [bin2, bin2, bin2, bin2, bin2, bin2, sumbin,sumbin,sumbin,sumbin,sumbin,sumbin]
        #binnings    = [bin2, bin2, bin2, bin2, bin2, bin2]
        #binnings    = [binsig, binsig, binsig, binsig]
        #binnings    = [bin2, bin2, bin2, bin2, bin2, bin2, bin2, bin2, bin2, bin2, bin2, bin2,bin2, bin2, bin2, bin2, bin2, bin2,bin2, bin2, bin2, bin2, bin2, bin2]
        binnings    = [bin1]
        dependence = 'zll_pt'
        ZtoLL = Region.region('Zto',
                           cuts.leps(), doVariables, dependence, 'met', binnings, True)
        regions.append(ZtoLL)                                                                                   
    else:
        #doVariables = [ 'met_uPerp']
        #doVariables = ['nVert', 'gamma_pt', 'met_pt', 'metPuppi_pt',  'met_rawPt', 'metPuppi_rawPt',  'metPuppi_pt']
        #doVariables = [ 'gamma_pt', 'nVert',  'met_pt', 'metPuppi_pt',  'met_rawPt', 'metPuppi_rawPt',  'metPuppi_pt']
        #doVariables = [ 'met_pt', 'met_rawPt','metPuppi_rawPt']
        #doVariables = [ 'met_x',  'met_y','met_sumEt',  'met_rawSumEt']
        #doVariables = [  'met_uPara']
        doVariables = [ 'met_uPerp', 'met_uPara','metRaw_uPerp',  'metRaw_uPara','metPuppi_uPerp', 'metPuppi_uPara']
        #doVariables = [ 'metPuppi_uPerp','metPuppi_uPara']
        #binnings    = [bin2, bin2]
        #binnings    = [bin2]
        #binnings    = [nbin,gammabin,  binx, binx, binx, binx, binx]
        #binnings    = [bin1, bin1, bin1]
        #binnings    = [bin2, bin2, sumbin,sumbin]
        binnings    = [bin2, bin2, bin2, bin2, bin2, bin2]
        dependence = 'gamma_pt'
        Gamma = Region.region('Gamma',
                           cuts.gammas(), doVariables, dependence, 'met', binnings, True)
        regions.append(Gamma)                                                                                   

    for reg in regions:
        print color.bold+color.red+'=========================================='
        print 'I am doing', reg.name, channel
        print '=========================================='+color.end
                    
        for var in reg.rvars:  
            Variable = []
            noRatio = False 
            justMET = False 
            doUperp = False 
            doUpara = False
            isLog = 1; isSig = 0; fixAxis = 0;option ='met'
            mc_histo = 0; histo_err = 0; mc_up = 0; mc_down = 0; mc_jerup = 0; mc_jerdown = 0; mc_unclUp = 0; mc_unclDown = 0; mc_stack = r.THStack();
            cut = reg.cuts
            if doDY:
                leg = [0.7, 0.6, 0.88, 0.9]
            else:
                leg = [0.65, 0.6, 0.81, 0.9]
            jetCut = "(nJet40 >= 0)"
            print color.blue+'************************************************************************************************'+color.end
            print 'loading variable %s '%(var)
            # for all the variables with some met in them, we need to run over the jec/unclustered energy to get the errors, that is why there is some nasty repeating of lines below, and both the emt_pt and met_phi is needed for the uPara and uPerp, and the justMET, doUperp and doUpara are set, along with some options used in the plotting, definitions of the option can be found in the include/Canvas.py where all the plotting style is set
            if var == 'met_pt':
                varTitle    = 'E_{T}^{miss} [GeV]';justMET = True;
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt', 'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt','met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
            elif var == 'met_pt_0jet':
                varTitle    = 'E_{T}^{miss} [GeV] (0 jet)';justMET = True;jetCut = (" (nJet40 == 0)")
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt', 'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt','met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
            elif var == 'met_pt_1jet':
                varTitle    = 'E_{T}^{miss} [GeV] (1 jet)';justMET = True;  jetCut = (" (nJet40 >= 1)") 
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt', 'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt','met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
            elif var == 'metPuppi_pt_0jet':
                varTitle    = 'Puppi E_{T}^{miss} [GeV] (0 jet)';justMET = True;jetCut = (" (nJet40 == 0)")
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt','metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
            elif var == 'metPuppi_pt_1jet':
                varTitle    = 'Puppi E_{T}^{miss} [GeV] (1 jet)';justMET = True;  jetCut = (" (nJet40 >= 1)") 
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt','metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
            elif var == 'met_x':
                varTitle    = 'E_{x}^{miss} [GeV]';justMET = True;  
                met = ['met_pt*sin(met_phi)', 'met_jecUp_pt*sin(met_jecUp_phi)', 'met_jecDown_pt*sin(met_jecDown_phi)', 'met_shifted_JetResUp_pt*sin(met_shifted_JetResUp_phi)', 'met_shifted_JetResDown_pt*sin(met_shifted_JetResDown_phi)','met_shifted_UnclusteredEnUp_pt*sin(met_shifted_UnclusteredEnUp_phi)', 'met_shifted_UnclusteredEnDown_pt*sin(met_shifted_UnclusteredEnDown_phi)']
            elif var == 'met_x_0jet':
                varTitle    = 'E_{x}^{miss} [GeV] (0 jet)';justMET = True;jetCut = (" (nJet40 == 0)")  
                met = ['met_pt*sin(met_phi)', 'met_jecUp_pt*sin(met_jecUp_phi)', 'met_jecDown_pt*sin(met_jecDown_phi)', 'met_shifted_JetResUp_pt*sin(met_shifted_JetResUp_phi)', 'met_shifted_JetResDown_pt*sin(met_shifted_JetResDown_phi)','met_shifted_UnclusteredEnUp_pt*sin(met_shifted_UnclusteredEnUp_phi)', 'met_shifted_UnclusteredEnDown_pt*sin(met_shifted_UnclusteredEnDown_phi)']
            elif var == 'met_x_1jet':
                varTitle    = 'E_{x}^{miss} [GeV] (1 jet)';justMET = True;jetCut = (" (nJet40 >= 1)") 
                met = ['met_pt*sin(met_phi)', 'met_jecUp_pt*sin(met_jecUp_phi)', 'met_jecDown_pt*sin(met_jecDown_phi)', 'met_shifted_JetResUp_pt*sin(met_shifted_JetResUp_phi)', 'met_shifted_JetResDown_pt*sin(met_shifted_JetResDown_phi)','met_shifted_UnclusteredEnUp_pt*sin(met_shifted_UnclusteredEnUp_phi)', 'met_shifted_UnclusteredEnDown_pt*sin(met_shifted_UnclusteredEnDown_phi)']
            elif var == 'met_y':
                varTitle    = 'E_{y}^{miss} [GeV]';justMET = True;
                met = ['met_pt*cos(met_phi)', 'met_jecUp_pt*cos(met_jecUp_phi)', 'met_jecDown_pt*cos(met_jecDown_phi)', 'met_shifted_JetResUp_pt*cos(met_shifted_JetResUp_phi)', 'met_shifted_JetResDown_pt*cos(met_shifted_JetResDown_phi)','met_shifted_UnclusteredEnUp_pt*cos(met_shifted_UnclusteredEnUp_phi)', 'met_shifted_UnclusteredEnDown_pt*cos(met_shifted_UnclusteredEnDown_phi)']
            elif var == 'met_y_0jet':
                varTitle    = 'E_{y}^{miss} [GeV] (0 jet)';jetCut = (" (nJet40 == 0)");justMET = True;
                met = ['met_pt*cos(met_phi)', 'met_jecUp_pt*cos(met_jecUp_phi)', 'met_jecDown_pt*cos(met_jecDown_phi)', 'met_shifted_JetResUp_pt*cos(met_shifted_JetResUp_phi)', 'met_shifted_JetResDown_pt*cos(met_shifted_JetResDown_phi)','met_shifted_UnclusteredEnUp_pt*cos(met_shifted_UnclusteredEnUp_phi)', 'met_shifted_UnclusteredEnDown_pt*cos(met_shifted_UnclusteredEnDown_phi)']
            elif var == 'met_y_1jet':
                varTitle    = 'E_{y}^{miss} [GeV] (1 jet)';jetCut = (" (nJet40 >= 1)"); justMET = True;
                met = ['met_pt*cos(met_phi)', 'met_jecUp_pt*cos(met_jecUp_phi)', 'met_jecDown_pt*cos(met_jecDown_phi)', 'met_shifted_JetResUp_pt*cos(met_shifted_JetResUp_phi)', 'met_shifted_JetResDown_pt*cos(met_shifted_JetResDown_phi)','met_shifted_UnclusteredEnUp_pt*cos(met_shifted_UnclusteredEnUp_phi)', 'met_shifted_UnclusteredEnDown_pt*cos(met_shifted_UnclusteredEnDown_phi)']
            elif var == 'met_sumEt':
                varTitle    = '#Sigma E_{T} [GeV]';justMET = True;
                met = ['met_sumEt', 'met_jecUp_sumEt', 'met_jecDown_sumEt', 'met_shifted_JetResUp_sumEt', 'met_shifted_JetResDown_sumEt','met_shifted_UnclusteredEnUp_sumEt', 'met_shifted_UnclusteredEnDown_sumEt']
            elif var == 'met_rawPt': 
                varTitle    = 'Raw E_{T}^{miss} [GeV]';justMET = True;
                met = ['met_rawPt', 'met_jecUp_rawPt', 'met_jecDown_rawPt', 'met_shifted_JetResUp_rawPt', 'met_shifted_JetResDown_rawPt','met_shifted_UnclusteredEnUp_rawPt', 'met_shifted_UnclusteredEnDown_rawPt']
            elif var == 'met_rawPt_0jet': 
                varTitle    = 'Raw E_{T}^{miss} [GeV] (0 jet)';justMET = True;jetCut = (" (nJet40 == 0)")
                met = ['met_rawPt', 'met_jecUp_rawPt', 'met_jecDown_rawPt', 'met_shifted_JetResUp_rawPt', 'met_shifted_JetResDown_rawPt','met_shifted_UnclusteredEnUp_rawPt', 'met_shifted_UnclusteredEnDown_rawPt']
            elif var == 'met_rawPt_1jet':
                varTitle    = 'Raw E_{T}^{miss} [GeV] (1 jet)';justMET = True; jetCut = (" (nJet40 >= 1)")
                met = ['met_rawPt', 'met_jecUp_rawPt', 'met_jecDown_rawPt', 'met_shifted_JetResUp_rawPt', 'met_shifted_JetResDown_rawPt','met_shifted_UnclusteredEnUp_rawPt', 'met_shifted_UnclusteredEnDown_rawPt']
            elif var == 'met_rawSumEt':
                varTitle    = 'Raw #Sigma E_{T} [GeV]' ;justMET = True;
                met = ['met_rawSumEt', 'met_jecUp_rawSumEt', 'met_jecDown_rawSumEt', 'met_shifted_JetResUp_rawSumEt', 'met_shifted_JetResDown_rawSumEt','met_shifted_UnclusteredEnUp_rawSumEt', 'met_shifted_UnclusteredEnDown_rawSumEt']
            elif var == 'metPuppi_rawPt':
                varTitle    = 'Puppi Raw E_{T}^{miss} [GeV]' ; justMET = True;option = "puppi" 
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt','metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
            elif var == 'metPuppi_pt':
                varTitle    = 'Puppi E_{T}^{miss} [GeV]' ; justMET = True;option = "puppi" 
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt','metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
            elif var == 'metPuppi_rawPt_0jet':
                varTitle    = 'Puppi Raw E_{T}^{miss} [GeV] (0 jet)' ;justMET = True;option = "puppi";jetCut = (" (nJet40 == 0)")
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt','metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
            elif var == 'metPuppi_rawPt_1jet':
                varTitle    = 'Puppi Raw E_{T}^{miss} [GeV] (1 jet)';justMET = True;option = "puppi";jetCut = (" (nJet40 >= 1)")
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt','metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
            elif var == 'met_uPerp':
                varTitle    = 'u_{#perp}'+"  "+ ' [GeV]' ;justMET = False;doUperp = True;
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt', 'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt', 'met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
                phi = ['met_phi', 'met_jecUp_phi', 'met_jecDown_phi','met_shifted_JetResUp_phi','met_shifted_JetResDown_phi',  'met_shifted_UnclusteredEnUp_phi', 'met_shifted_UnclusteredEnDown_phi']
            elif var == 'met_uPerp_0jet':
                varTitle    = 'u_{#perp}'+"  "+ ' [GeV] (0 jet)' ;justMET = False;doUperp = True;jetCut = (" (nJet40 == 0)")
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt', 'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt', 'met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
                phi = ['met_phi', 'met_jecUp_phi', 'met_jecDown_phi','met_shifted_JetResUp_phi','met_shifted_JetResDown_phi',  'met_shifted_UnclusteredEnUp_phi', 'met_shifted_UnclusteredEnDown_phi']
            elif var == 'met_uPerp_1jet':
                varTitle    = 'u_{#perp}'+"  "+ ' [GeV] (1 jet)';justMET = False;doUperp = True;jetCut = (" (nJet40 >= 1)")
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt', 'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt', 'met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
                phi = ['met_phi', 'met_jecUp_phi', 'met_jecDown_phi','met_shifted_JetResUp_phi','met_shifted_JetResDown_phi',  'met_shifted_UnclusteredEnUp_phi', 'met_shifted_UnclusteredEnDown_phi']
            elif var == 'met_uPara':
                varTitle =  'u_{||} + q_{T} [GeV]';justMET = False;doUpara = True; 
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt',  'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt', 'met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
                phi = ['met_phi', 'met_jecUp_phi', 'met_jecDown_phi', 'met_shifted_JetResUp_phi','met_shifted_JetResDown_phi', 'met_shifted_UnclusteredEnUp_phi', 'met_shifted_UnclusteredEnDown_phi']
            elif var == 'met_uPara_0jet':
                varTitle =  'u_{||} + q_{T} [GeV] (0 jet)';justMET = False;doUpara = True;jetCut = (" (nJet40 >= 0)")
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt',  'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt', 'met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
                phi = ['met_phi', 'met_jecUp_phi', 'met_jecDown_phi', 'met_shifted_JetResUp_phi','met_shifted_JetResDown_phi', 'met_shifted_UnclusteredEnUp_phi', 'met_shifted_UnclusteredEnDown_phi']
            elif var == 'met_uPara_1jet':
                varTitle =  'u_{||} + q_{T} [GeV] (1 jet)';justMET = False;doUpara = True;jetCut = (" (nJet40 >= 1)")
                met = ['met_pt', 'met_jecUp_pt', 'met_jecDown_pt',  'met_shifted_JetResUp_pt', 'met_shifted_JetResDown_pt', 'met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnDown_pt']
                phi = ['met_phi', 'met_jecUp_phi', 'met_jecDown_phi', 'met_shifted_JetResUp_phi','met_shifted_JetResDown_phi', 'met_shifted_UnclusteredEnUp_phi', 'met_shifted_UnclusteredEnDown_phi']
            elif var == 'metRaw_uPerp':
                varTitle    = 'Raw u_{#perp}'+"  "+ ' [GeV]';justMET = False;doUperp = True;
                met=['met_rawPt','met_jecUp_rawPt','met_jecDown_rawPt','met_shifted_JetResUp_rawPt','met_shifted_JetResDown_rawPt','met_shifted_UnclusteredEnUp_rawPt','met_shifted_UnclusteredEnDown_rawPt']
                phi=['met_rawPhi','met_jecUp_rawPhi','met_jecDown_rawPhi','met_shifted_JetResUp_rawPhi','met_shifted_JetResDown_rawPhi','met_shifted_UnclusteredEnUp_rawPhi','met_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metRaw_uPerp_0jet':
                varTitle    = 'Raw u_{#perp}'+"  "+ ' [GeV] (0 jet)';justMET = False;doUperp = True;jetCut = (" (nJet40 == 0)")
                met=['met_rawPt','met_jecUp_rawPt','met_jecDown_rawPt','met_shifted_JetResUp_rawPt','met_shifted_JetResDown_rawPt','met_shifted_UnclusteredEnUp_rawPt','met_shifted_UnclusteredEnDown_rawPt']
                phi=['met_rawPhi','met_jecUp_rawPhi','met_jecDown_rawPhi','met_shifted_JetResUp_rawPhi','met_shifted_JetResDown_rawPhi','met_shifted_UnclusteredEnUp_rawPhi','met_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metRaw_uPerp_1jet':
                varTitle    = 'Raw u_{#perp}'+"  "+ ' [GeV] (1 jet)';justMET = False;doUperp = True;jetCut = (" (nJet40 >= 1)")
                met=['met_rawPt','met_jecUp_rawPt','met_jecDown_rawPt','met_shifted_JetResUp_rawPt','met_shifted_JetResDown_rawPt','met_shifted_UnclusteredEnUp_rawPt','met_shifted_UnclusteredEnDown_rawPt']
                phi=['met_rawPhi','met_jecUp_rawPhi','met_jecDown_rawPhi','met_shifted_JetResUp_rawPhi','met_shifted_JetResDown_rawPhi','met_shifted_UnclusteredEnUp_rawPhi','met_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metRaw_uPara':                                                                                                                                                                                     
                varTitle =  'Raw u_{||} + q_{T} [GeV]';justMET = False;doUpara = True;
                met = ['met_rawPt', 'met_jecUp_rawPt', 'met_jecDown_rawPt','met_shifted_JetResUp_rawPt','met_shifted_JetResDown_rawPt', 'met_shifted_UnclusteredEnUp_rawPt', 'met_shifted_UnclusteredEnDown_rawPt']
                phi = ['met_rawPhi', 'met_jecUp_rawPhi', 'met_jecDown_rawPhi','met_shifted_JetResUp_rawPhi','met_shifted_JetResDown_rawPhi',  'met_shifted_UnclusteredEnUp_rawPhi', 'met_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metRaw_uPara_0jet':
                varTitle =  'Raw u_{||} + q_{T} [GeV] (0 jet)';justMET = False;doUpara = True;jetCut = (" (nJet40 == 0)")
                met = ['met_rawPt', 'met_jecUp_rawPt', 'met_jecDown_rawPt','met_shifted_JetResUp_rawPt','met_shifted_JetResDown_rawPt', 'met_shifted_UnclusteredEnUp_rawPt', 'met_shifted_UnclusteredEnDown_rawPt']
                phi = ['met_rawPhi', 'met_jecUp_rawPhi', 'met_jecDown_rawPhi','met_shifted_JetResUp_rawPhi','met_shifted_JetResDown_rawPhi',  'met_shifted_UnclusteredEnUp_rawPhi', 'met_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metRaw_uPara_1jet':
                varTitle =  'Raw u_{||} + q_{T} [GeV] (1 jet)';justMET = False;doUpara = True;jetCut = (" (nJet40 >= 1)")
                met = ['met_rawPt', 'met_jecUp_rawPt', 'met_jecDown_rawPt','met_shifted_JetResUp_rawPt','met_shifted_JetResDown_rawPt', 'met_shifted_UnclusteredEnUp_rawPt', 'met_shifted_UnclusteredEnDown_rawPt']
                phi = ['met_rawPhi', 'met_jecUp_rawPhi', 'met_jecDown_rawPhi','met_shifted_JetResUp_rawPhi','met_shifted_JetResDown_rawPhi',  'met_shifted_UnclusteredEnUp_rawPhi', 'met_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppiRaw_uPerp':
                varTitle    = 'Raw Puppi u_{#perp}'+"  "+ ' [GeV]';justMET = False; doUperp = True;option = "puppi"
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt' , 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
                phi = ['metPuppi_rawPhi', 'metPuppi_jecUp_rawPhi', 'metPuppi_jecDown_rawPhi', 'metPuppi_shifted_JetResUp_rawPhi', 'metPuppi_shifted_JetResDown_rawPhi', 'metPuppi_shifted_UnclusteredEnUp_rawPhi', 'metPuppi_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppiRaw_uPerp_0jet':
                varTitle    = 'Raw Puppi u_{#perp}'+"  "+ ' [GeV] (0 jet)';justMET = False; doUperp = True;option = "puppi";jetCut = (" (nJet40 == 0)")
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt' , 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
                phi = ['metPuppi_rawPhi', 'metPuppi_jecUp_rawPhi', 'metPuppi_jecDown_rawPhi', 'metPuppi_shifted_JetResUp_rawPhi', 'metPuppi_shifted_JetResDown_rawPhi', 'metPuppi_shifted_UnclusteredEnUp_rawPhi', 'metPuppi_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppiRaw_uPerp_1jet':                                                                                                                    
                varTitle    = 'Raw Puppi u_{#perp}'+"  "+ ' [GeV] (1 jet)';justMET = False; doUperp = True;option = "puppi";jetCut = (" (nJet40 >= 1)")
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt' , 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
                phi = ['metPuppi_rawPhi', 'metPuppi_jecUp_rawPhi', 'metPuppi_jecDown_rawPhi', 'metPuppi_shifted_JetResUp_rawPhi', 'metPuppi_shifted_JetResDown_rawPhi', 'metPuppi_shifted_UnclusteredEnUp_rawPhi', 'metPuppi_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppiRaw_uPara':
                varTitle =  'Raw Puppi u_{||} + q_{T} [GeV]';justMET = False; doUpara = True;option = "puppi"
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
                phi = ['metPuppi_rawPhi', 'metPuppi_jecUp_rawPhi', 'metPuppi_jecDown_rawPhi', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPhi', 'metPuppi_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppi_uPara':
                varTitle =  'Puppi u_{||} + q_{T} [GeV]';justMET = False; doUpara = True;option = "puppi"
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
                phi = ['metPuppi_phi', 'metPuppi_jecUp_phi', 'metPuppi_jecDown_phi', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'metPuppi_shifted_UnclusteredEnDown_phi']
            elif var == 'metPuppi_uPerp':                                                                                                                                                                                           
                varTitle =  'Puppi u_{#perp} + q_{T} [GeV]';justMET = False; doUperp = True;option = "puppi"
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']

                phi = ['metPuppi_phi', 'metPuppi_jecUp_phi', 'metPuppi_jecDown_phi', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'metPuppi_shifted_UnclusteredEnDown_phi']
            elif var == 'metPuppiRaw_uPara_0jet':
                varTitle =  'Raw Puppi u_{||} + q_{T} [GeV] (0 jet)';justMET = False; doUpara = True;option = "puppi";jetCut = (" (nJet40 == 0)")
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
                phi = ['metPuppi_rawPhi', 'metPuppi_jecUp_rawPhi', 'metPuppi_jecDown_rawPhi', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPhi', 'metPuppi_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppiRaw_uPara_1jet':
                varTitle =  'Raw Puppi u_{||} + q_{T} [GeV] (1 jet)';justMET = False; doUpara = True;option = "puppi";jetCut = (" (nJet40 >= 1)")
                met = ['metPuppi_rawPt', 'metPuppi_jecUp_rawPt', 'metPuppi_jecDown_rawPt', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPt', 'metPuppi_shifted_UnclusteredEnDown_rawPt']
                phi = ['metPuppi_rawPhi', 'metPuppi_jecUp_rawPhi', 'metPuppi_jecDown_rawPhi', 'metPuppi_shifted_JetResUp_rawPt', 'metPuppi_shifted_JetResDown_rawPt', 'metPuppi_shifted_UnclusteredEnUp_rawPhi', 'metPuppi_shifted_UnclusteredEnDown_rawPhi']
            elif var == 'metPuppi_uPara_0jet':
                varTitle =  'Puppi u_{||} + q_{T} [GeV] (0 jet)';justMET = False; doUpara = True;option = "puppi";jetCut = (" (nJet40 == 0)")
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
                phi = ['metPuppi_phi', 'metPuppi_jecUp_phi', 'metPuppi_jecDown_phi', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'metPuppi_shifted_UnclusteredEnDown_phi']
            elif var == 'metPuppi_uPara_1jet':
                varTitle =  'Puppi u_{||} + q_{T} [GeV] (1 jet)';justMET = False; doUpara = True;option = "puppi";jetCut = (" (nJet40 >= 1)")
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
                phi = ['metPuppi_phi', 'metPuppi_jecUp_phi', 'metPuppi_jecDown_phi', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'metPuppi_shifted_UnclusteredEnDown_phi']  
            elif var == 'metPuppi_uPerp_0jet':
                varTitle =  'Puppi u_{#perp} + q_{T} [GeV] (0 jet)';justMET = False; doUperp = True;option = "puppi";jetCut = (" (nJet40 == 0)")
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
                phi = ['metPuppi_phi', 'metPuppi_jecUp_phi', 'metPuppi_jecDown_phi', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'metPuppi_shifted_UnclusteredEnDown_phi']
            elif var == 'metPuppi_uPerp_1jet':
                varTitle =  'Puppi u_{#perp} + q_{T} [GeV] (1 jet)';justMET = False; doUperp = True;option = "puppi";jetCut = (" (nJet40 >= 1)")
                met = ['metPuppi_pt', 'metPuppi_jecUp_pt', 'metPuppi_jecDown_pt', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnDown_pt']
                phi = ['metPuppi_phi', 'metPuppi_jecUp_phi', 'metPuppi_jecDown_phi', 'metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'metPuppi_shifted_UnclusteredEnDown_phi']  
            #miscellaneous plots without jec or unclustered errors
            elif var == 'nVert':
                varTitle = 'number of vertices'
                met = ['nVert']   
                option = 'nvert'
                isLog = 0
                noRatio = False
                leg = [0.68, 0.74, 0.82, 0.9]
            elif var == 'met_phi':
                varTitle = 'E_{T}^{miss} #phi'
                met = ['met_phi']       
            elif var == 'zll_pt':
                varTitle    = 'q_{T} [GeV]'
                met = ['zll_pt']               
                option = 'qt'                       
                noRatio = False
            elif var == 'gamma_r9':
                varTitle    = 'R9'
                met = ['gamma_r9']               
                option = 'chi'                       
            elif var == 'gamma_hOverE':
                varTitle    = 'H/E'
                met = ['gamma_hOverE']      
                option = 'chi'          
            elif var == 'gamma_sigmaIetaIeta':
                varTitle    = '#sigma i#eta i#eta '
                met = ['gamma_sigmaIetaIeta']      
                option = 'chi'          
            elif var == 'zll_mass':
                varTitle    = 'm_{ll} [GeV]'
                met = ['zll_mass']    
                option = 'mass'
            elif var == 'lep_pt':
                varTitle    = 'lep p_{T} [GeV]'
                met = ['lep_pt']                    
            elif var == 'lep_eta1':
                varTitle    = 'Leading lepton #eta'
                met = ['lep_eta[0]']                                   
            elif var == 'lep_eta2':
                varTitle    = 'Subleading lepton #eta'
                met = ['lep_eta[1]']                    
            elif var == 'met_sig':
                varTitle    = 'E_{T}^{miss} Significance [GeV]'
                met = ['met_sig']    
                isSig = 1
                option = 'sig0'
                leg = [0.68, 0.6, 0.82, 0.88]
            elif var == 'met_sig_1jet':
                varTitle    = 'E_{T}^{miss} Significance [GeV]'
                met = ['met_sig']   
                isSig = 1
                option = 'sig1'
                jetCut = ("(nJet40 >= 1)")
                leg = [0.68, 0.6, 0.82, 0.88]
            elif var == 'metPuppi_sig':
                varTitle    = 'E_{T}^{miss} Significance [GeV]'
                met = ['metPuppi_sig']            
            elif var == 'metPuppi_sig_1jet':
                varTitle    = 'Puppi E_{T}^{miss} Significance [GeV]'
                met = ['metPuppi_sig']            
                option = 'sig1'
                jetCut = ("(nJet40 >= 1)")
            elif var == 'gamma_pt':
                varTitle    = 'q_{T} [GeV]'
                met = ['gamma_pt']      
                option ='qt'
                noRatio = False
                leg = [0.66, 0.7, 0.8, 0.88]
            elif var == 'chi2':
                varTitle = "#chi^{2} Probability"
                met =  ["TMath::Prob(met_sig,2)"]
                isLog = 0
                option = 'chi'
            elif var == 'chi2_1jet':
                varTitle = "#chi^{2} Probability , p_{T}^{jet1} > 50 GeV"
                met =  [" TMath::Prob(met_sig,2)"]
                cut = ("((jet1_pt >= 0))&&" + reg.cuts)
                isLog = 0
                option = 'chi'
            tmp_histo = 0
            histo_err = 0
            tmp_full = 0
            for m in met:
                Variable = ""
                # here, either make the met, uPara or uPerp, and loop over all five types of met...
                if doUpara:
                    Variable = makeUPara(met[met.index(m)], phi[met.index(m)], boson, boson_phi)
                elif doUperp:
                    Variable = makeUPerp(met[met.index(m)], phi[met.index(m)], boson, boson_phi)
                elif justMET:                                                           
                    Variable = makeMET(met[met.index(m)])
                else:
                    # ... or just take the first element in the list (for example for nvert, pt, eta plots etc.)
                    Variable = met[0]
                if met.index(m) == 0:
                    data_hist = treeDA.getTH1F(lumi, var+"_"+reg.name+'data'+str(met.index(m)), Variable, reg.bins[reg.rvars.index(var)], 1, 1, cuts.Add(cut, jetCut) , "", varTitle, doReReco)
                    print 'data ', data_hist.Integral()
                    if doHighMetValues:
                        highMet = data_hist.FindLastBinAbove(0, 1)
                        print "high metbin data: ", highMet                 

                for tree in ( mcTrees):
                    ind = 0
                    cuts = CutManager.CutManager()
                    treename = tree.name.lower()
                    print 'with  %s' %treename
                    block = tree.blocks[0]
                    attr =  var+''
                    if met.index(m) == 0:
                        # this option here treats the first variable in the list met, which is either just the regular met, or the only element in a list, for nvert, pt eta, etc. This should get the TH1 for data and mc. For mc, both a full histo and a stack is made. 
                        tmp_full= tree.getTH1F(lumi, var+"_"+reg.name+treename+str(met.index(m)),  Variable, reg.bins[reg.rvars.index(var)], 1, 1, cuts.Add(cut, jetCut) , "", varTitle, doReReco)
                        tmp_full.SetFillColorAlpha(block.color, 0.5)
                        tmp_full.SetTitle(block.name)
                        if treename == 'gjets':
                            tmp_full.SetTitle("#gamma + jets")
                        if treename == 'dy':
                            if doee:
                                tmp_full.SetTitle("Z #rightarrow ee")
                            else:                                    
                                tmp_full.SetTitle("Z #rightarrow #mu#mu")
                        if treename == 'tt':
                            tmp_full.SetTitle("Top")
                        getattr(reg, attr).setHisto(tmp_full, 'MC', 'central')
                        tmp_histo = copy.deepcopy(tmp_full.Clone(var+reg.name))
                        tmp_histo.GetXaxis().SetTitle(varTitle)
                        mc_stack.Add(tmp_histo)
                        if not ind: mc_stack.Draw()
                        ind+=1
                        mc_stack.GetXaxis().SetTitle(tmp_histo.GetXaxis().GetTitle())
                        print treename,  'histo has integral %.2f'%tmp_histo.Integral()
                        if doHighMetValues:
                            highMet = tmp_histo.FindLastBinAbove(0, 1)
                            print "high metbin: ", highMet                 
                            
                        if not mc_histo:
                            mc_histo = copy.deepcopy(tmp_histo)
                        else:
                            mc_histo.Add(tmp_histo, 1.)                
                    if len(met)>1: # when the array is longer than one, do jec and met unclustered errors, and just get the histo, no stack obvi, and this is of course only done for mc.
                        if met.index(m) > 0:
                            histo_err= tree.getTH1F(lumi, var+"_"+reg.name+treename+str(met.index(m)),  Variable, reg.bins[reg.rvars.index(var)], 1, 1, cut , "", varTitle, doReReco)
                            SetOwnership(histo_err, 0 )
                        if met.index(m) == 1:
                            if not mc_up:                                    
                                mc_up = copy.deepcopy(histo_err);SetOwnership(mc_up, 0 )
                            else:
                                mc_up.Add(histo_err)
                            print 'doing mc_up'
                        if met.index(m) == 2:
                            if not mc_down:
                                mc_down = copy.deepcopy(histo_err);SetOwnership(mc_down, 0 )
                            else:
                                mc_down.Add(histo_err)
                            print 'doing mc_down'                                                         
                        if met.index(m) == 3:
                            if not mc_jerup:                                    
                                mc_jerup = copy.deepcopy(histo_err);SetOwnership(mc_jerup, 0 )
                            else:
                                mc_jerup.Add(histo_err)
                            print 'doing mc_JERup'
                        if met.index(m) == 4:
                            if not mc_jerdown:
                                mc_jerdown = copy.deepcopy(histo_err);SetOwnership(mc_jerdown, 0 )
                            else:
                                mc_jerdown.Add(histo_err)
                            print 'doing mc_JERdown'                                                         
                        if met.index(m) == 5:
                            if not mc_unclUp:                                    
                                mc_unclUp = copy.deepcopy(histo_err);SetOwnership(mc_unclUp, 0 )
                            else:
                                mc_unclUp.Add(histo_err)
                            print 'doing mc_unclUp'
                        if met.index(m) == 6:
                            if not mc_unclDown:
                                mc_unclDown = copy.deepcopy(histo_err);SetOwnership(mc_unclDown, 0 )
                            else:
                                mc_unclDown.Add(histo_err)
                            print 'doing mc_unclDown'
                    else: # else, just do stat errors
                        mc_up = mc_histo; mc_down = mc_histo; mc_jerup = mc_histo; mc_jerdown = mc_histo;mc_unclUp = mc_histo; mc_unclDown = mc_histo
            SetOwnership(mc_stack, 0 )                                                                                                                       
            SetOwnership(tmp_histo, 0 )
            SetOwnership(mc_histo, 0 )
            print 'plotting ', var
            if isSig:
                # this is an option for when the met significance plot is made, which make the little chi2 line.
                print "bin cont ", mc_histo.Integral() 
                f = r.TF1("f", str(mc_histo.Integral())+"*2*TMath::Prob(x, 2)", 1, 40) 
                plot_var = Canvas.Canvas("test/%s/%s_%s%s"%(run_str, var,reg.name, channel), "png,root,pdf,C", leg[0], leg[1], leg[2], leg[3])
            else:
                if doReReco:
                    rereco = 'ReMINIAOD'
                else:
                    rereco = 'Prompt'
                plot_var = Canvas.Canvas("test/%s/%s_%s%s_%s"%(run_str, var,reg.name, channel, rereco), "png,root,pdf,C", leg[0], leg[1], leg[2], leg[3])
            plot_var.addStack(mc_stack  , "hist" , 1, 1)
            if isSig:
                plot_var.addHisto(f, "E,SAME"   , "#chi^{2} d.o.f. 2"  , "L", r.kBlack , 1, 0)
            plot_var.addHisto(data_hist, "E,SAME"   , "Data"  , "PL", r.kBlack , 1, 0)
            if noRatio:
                plot_var.save(1,  isLog, lumi,  "", "", varTitle, option)
            else:
                plot_var.saveRatio(1,1, isLog, lumi, data_hist, mc_histo, mc_up, mc_down, mc_jerup, mc_jerdown, mc_unclUp, mc_unclDown, varTitle, option, run_str)
            del plot_var
            del Variable
            time.sleep(0.1)
            print color.blue+'********************************************DONE***************************************************'+color.end
