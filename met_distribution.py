#####################################################################
######                                                              #
###### 88888888888         88                        888888888888   #  
###### 88                  88                                 ,88   #
###### 88                  88                               ,88"    #  
###### 88aaaaa     ,adPPYb,88  ,adPPYb,d8  ,adPPYba,      ,88"      #
###### 88"""""    a8"    `Y88 a8"    `Y88 a8P_____88    ,88"        #
###### 88         8b       88 8b       88 8PP"""""""  ,88"          #
###### 88         "8a,   ,d88 "8a,   ,d88 "8b,   ,aa 88"            #
###### 88888888888 `"8bbdP"Y8  `"YbbdP"Y8  `"Ybbd8"' 888888888888   #
######                       aa,    ,88                             #
######                         "Y8bbdP"                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors
import math, sys, optparse, array, time, copy

import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample
import include.Rounder    as Rounder


if __name__ == "__main__":


    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples2.dat', help='the samples file. default \'samples2.dat\'')
    (opts, args) = parser.parse_args()
    doDY = False
    doee = False
    print 'Going to load DATA and MC trees...'
    if doDY:
        if doee:
            ttDatasets = ['TTJets_DiLepton_ee']
            dyDatasets = ['DYJetsToLL_M50_ee']
            ewkDatasets = ['WWTo2L2Nu_ee', 'WZTo2L2Q_ee', 'WZTo3LNu_ee', 'ZZTo4L_ee', 'WJetsToLNu_ee', 'ZZTo2L2Q_ee']
            daDatasets = ['DoubleEG_Run2016B_PromptReco_v2NEW']  
            channel = 'EE'
        else:
            ttDatasets = ['TTJets_DiLepton_mm']
            dyDatasets = ['DYJetsToLL_M50_mm']
            ewkDatasets = ['WWTo2L2Nu_mm', 'WZTo2L2Q_mm', 'WZTo3LNu_mm', 'ZZTo4L_mm', 'WJetsToLNu_mm', 'ZZTo2L2Q_mm']
            daDatasets = ['DoubleMuon_Run2016B_PromptReco_v2NEW']  
            channel = 'MM'
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, ttDatasets, 'TT'), 'TT'  , 0)
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dyDatasets, 'DY'), 'DY'  , 0)
        treeEWK = Sample.Tree(helper.selectSamples(opts.sampleFile, ewkDatasets, 'EWK'), 'EWK'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
        mcTrees = [ treeTT,treeEWK,  treeDY]                                                                   
    
    else:
        
        qcdDatasets = [ 'QCD_HT200to300_ext','QCD_HT300to500', 'QCD_HT500to700', 'QCD_HT700to1000', 'QCD_HT1000to1500', 'QCD_HT1500to2000', 'QCD_HT2000toInf']
        gjetsDatasets = ['GJets_HT40to100', 'GJets_HT100to200', 'GJets_HT200to400', 'GJets_HT400to600' ,  'GJets_HT600toInf']
        wgDatasets = [ 'WGToLNuG']
        wjDatasets = ['WJetsToLNu_HT100to200_ext', 'WJetsToLNu_HT200to400', 'WJetsToLNu_HT400to600', 'WJetsToLNu_HT600to800', 'WJetsToLNu_HT800to1200', 'WJetsToLNu_HT1200to2500', 'WJetsToLNu_HT2500toInf']
        zgDatasets = ['ZGJets','ZGJets40-130', 'ZGTo2LG']
        ttgDatasets = ['TTGJets']
        daDatasets = ['SinglePhoton_Run2016B_PromptReco_v2Golden']
        channel = '' 
        treeQCD =   Sample.Tree(helper.selectSamples(opts.sampleFile, qcdDatasets, 'QCD'), 'QCD', 0)
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'GJETS'),'GJETS', 0)
        treeWG =   Sample.Tree(helper.selectSamples(opts.sampleFile, wgDatasets, 'WG'), 'WG', 0)
        treeWJ =   Sample.Tree(helper.selectSamples(opts.sampleFile, wjDatasets, 'WJ'), 'WJ', 0)
        treeZG =   Sample.Tree(helper.selectSamples(opts.sampleFile, zgDatasets, 'ZG'), 'ZG', 1)
        treeTTG =   Sample.Tree(helper.selectSamples(opts.sampleFile, ttgDatasets, 'TTG'), 'TTG', 0)
        treeDA =    Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'DA'), 'DATA', 1)
        mcTrees = [  treeTTG, treeZG, treeWJ,treeWG , treeQCD, treeGJETS]                                                                   
        
    print 'Trees successfully loaded...'

    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ####Cuts needed by rmue
    cuts = CutManager.CutManager()

    color = helper.color

    lumi = 0.803015796
    lumi_str = 'lumi'+str(lumi)


    regions = []
    setLog = []
    etabins = [-1.5,-1.4,-1.3, -1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
    phibins = [-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3, -1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5]
    if doDY:
        #doVariables = ['met_zll']
        #doVariables = ['nVert']
        doVariables = [ 'met_uPerp_zll', 'met_uPara_zll']
        #doVariables = ['met_zll','met_uPerp_zll', 'met_uPara_zll', 'nVert']
        #binnings    = [ range(0,200, 5)]
        #binnings    = [ range(0,50, 1)]
        binnings    = [ range(-200,200, 5), range(-200, 200, 5)]
        dependence = 'zll_pt'
        variablename = 'met'
        ZtoLL = Region.region('Zto',
                           cuts.leps(), doVariables, dependence, variablename, binnings, True)
        regions.append(ZtoLL)                                                                                   
    else:
        #doVariables = [ 'nVert']
        #doVariables = [ 'met_uPerp_gamma']
        doVariables = ['met', 'nVert']
        binnings    = [range(0,200, 5), range(0,50, 1)]
        #binnings    = [ range(-200,200, 5)]
        #binnings    = [ range(0,50, 1)]
        dependence = 'gamma_pt'
        variablename = 'met'
        GammaJets = Region.region('GammaJets',
                           cuts.gammas(), doVariables, dependence, variablename, binnings, True)
        regions.append(GammaJets)                                                                                   


    etas = ['central' ]
        

    for reg in regions:
        print color.bold+color.red+'=========================================='
        print 'i am at region', reg.name
        print '=========================================='+color.end
        #for eta in ['central', 'forward']:
        for eta in etas:
                    

            for var in reg.rvars:

                if   var == 'met_uPerp_zll':
                    varTitle    = 'u_{#perp}'+"  "+ ' [GeV]'
                    varVariable = '((-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi) - (-met_pt*sin(met_phi) - zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi))/zll_pt'                                                                                               
                elif var == 'met_uPara_zll':
                    varTitle    = 'u_{||} + Z p_{T} [GeV]'
                    varVariable = '((-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi) + (-met_pt*sin(met_phi) - zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi))/zll_pt+zll_pt'
                elif var == 'met_zll':
                    varTitle    = 'E_{T}^{miss} [GeV]'
                    varVariable = 'met_pt'     
                elif var == 'met_uPerp_gamma':
                    varTitle    = 'u_{ #perp } '+"  "+ ' [GeV]'
                    varVariable = '((-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi) - (-met_pt*sin(met_phi) - gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt'
                elif var == 'met_uPara_gamma':
                    varTitle    = 'u_{||} + #gamma p_{T} [GeV]'
                    varVariable = '((-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi) + (-met_pt*sin(met_phi) - gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt+gamma_pt' 
                elif var == 'nVert':
                    varTitle    = 'nVert'
                    varVariable = 'nVert'         
                elif var == 'met':
                    varTitle    = 'E_{T}^{miss} [GeV]'
                    varVariable = 'met_pt'         
                elif var == 'zll_pt':
                    varTitle    = 'Z p_{T} [GeV]'
                    varVariable = 'zll_pt'         
                elif var == 'zll_mass':
                    varTitle    = 'm_{ll} [GeV]'
                    varVariable = 'zll_mass'         
                
                elif var == 'gamma_pt':
                    varTitle    = '#gamma p_{T} [GeV]'
                    varVariable = 'gamma_pt'         
                elif var == 'met_phi':
                    varTitle    = 'E_{T} #phi'
                    varVariable = 'met_phi'         
                elif var == 'met_eta':
                    varTitle    = 'ME_{T} #eta'
                    varVariable = 'met_eta'         

                print 'loading variable %s in %s'%(var, eta)
                mc_histo = 0
                mc_stack = r.THStack()
                print 'do data', reg.doData
                for tree in ( ([treeDA] +mcTrees) if reg.doData else mcTrees):
                    ind = 0
                    cuts = CutManager.CutManager()
                    treename = tree.name.lower()
                    print '... doing tree %s' %treename
                    dataMC = ('DATA' if tree == treeDA else 'MC'); isData = dataMC == 'DATA';
                    block = tree.blocks[0]
                    attr = var+('' if tree.name in ['DATA', 'MC'] else '_'+treename)
                    tmp_full= tree.getTH1F(lumi, var+"_"+eta+reg.name+treename, varVariable, reg.bins[reg.rvars.index(var)], 1, 1, reg.cuts , "", varTitle)
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
                        tmp_full.SetTitle("top")
                    getattr(reg, attr).setHisto(tmp_full, dataMC, eta)

                    ## don't do any adding for data
                    if isData: 
                        print "data", tmp_full.Integral()
                        continue

                    tmp_histo = copy.deepcopy(tmp_full.Clone(var+eta+reg.name))
                    tmp_histo.GetXaxis().SetTitle(varTitle)
                    mc_stack.Add(tmp_histo)
                    if not ind: mc_stack.Draw()
                    ind+=1
                    mc_stack.GetXaxis().SetTitle(tmp_histo.GetXaxis().GetTitle())
                    print treename,  'histo has integral %.2f'%tmp_histo.Integral()
                    if not mc_histo:
                        mc_histo = copy.deepcopy(tmp_histo)
                    else:
                        mc_histo.Add(tmp_histo, 1.)
                mc_histo.SetMinimum(0);
                setattr(reg, '%s_mc_histo_%s'%(var, eta), mc_histo)
                setattr(reg, '%s_mc_stack_%s'%(var, eta), mc_stack)

                print 'plotting %s in region %s in %s' %(var, reg.name, eta)
                plot_var = Canvas.Canvas("test/%s/%s_%s%s"%(lumi_str, var,reg.name, channel), "png,root", 0.7, 0.6, 0.9, 0.9)
                plot_var.addStack(mc_stack  , "hist" , 1, 1)
                plot_var.addHisto(getattr(reg, var).getHisto('DATA', eta), "E,SAME"   , "Data"  , "PL", r.kBlack , 1, 0)
                plot_var.saveRatio(1, 1, 1, lumi, getattr(reg, var).getHisto('DATA', eta), mc_histo )
                time.sleep(0.1)

