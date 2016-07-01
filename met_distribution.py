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
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()
    doDY = True
    doee = True
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
    if doDY:
        doVariables = ['met_uPerp_zll']
        #doVariables = [ 'met_pt','met_uPerp_zll', 'met_uPara_zll']
        binnings    = [range(-200, 200, 5)]
        #binnings    = [range(0,200, 5), range(-200,200, 5), range(-200,200, 5)]
        dependence = 'zll_pt'
        ZtoLL = Region.region('Zto',
                           cuts.leps(), doVariables, dependence, 'met', binnings, True)
        regions.append(ZtoLL)                                                                                   
    else:
        doVariables = ['met_uPerp_gamma', 'met_uPara_gamma']
        binnings    = [ range(-200,200, 5), range(-200,200, 5)]
        dependence = 'gamma_pt'
        Gamma = Region.region('Gamma',
                           cuts.gammas(), doVariables, dependence, 'met', binnings, True)
        regions.append(Gamma)                                                                                   

    etas = ['central' ]

    for reg in regions:
        print color.bold+color.red+'=========================================='
        print 'i am at region', reg.name
        print '=========================================='+color.end
        for eta in etas:
                    
            for var in reg.rvars:   

                if var == 'met_uPerp_zll':
                    jecs = ['met_uPerp_zll','met_uPerp_zll_Up', 'met_uPerp_zll_Down'] 
                elif var == 'met_uPara_zll':
                    jecs = ['met_uPara_zll','met_uPara_zll_Up', 'met_uPara_zll_Down'] 
                elif var == 'met_pt':
                    jecs = ['met_pt','met_Up', 'met_Down']    
                elif var == 'met_uPerp_gamma':
                    jecs = ['met_uPerp_gamma', 'met_uPerp_gamma_Up', 'met_uPerp_gamma_Down']
                elif var == 'met_uPara_gamma':
                    jecs = ['met_uPara_gamma', 'met_uPara_gamma_Up', 'met_uPara_gamma_Down']
                #the following variables will have stat errors only, that is why I loop over the same variables 3 times, to make only stat errors in Canvas.py
                elif var == 'nVert':
                    jecs  = ['nVert', 'nVertUp', 'nVertDown']
                elif var == 'zll_pt':
                    jecs = ['zll_pt', 'zll_pt_Up', 'zll_pt_Down']         
                elif var == 'zll_mass':
                    jecs  = ['zll_mass', 'zll_mass', 'zll_mass']         
                elif var == 'gamma_pt':
                    jecs = ['gamma_pt', 'gamma_pt', 'gamma_pt']         
                
                mc_histo = 0
                mc_up = 0
                mc_down = 0
                mc_stack = r.THStack()
                for jec in jecs:
                    print 'loading variable %s in %s'%(jec, eta)
                    if jec == 'met_pt':
                        varTitle    = 'E_{T}^{miss} [GeV]'
                        vari = 'met_pt'
                    elif jec == 'met_Up':
                        vari = 'met_jecUp_pt'
                    elif jec == 'met_Down':
                        vari = 'met_jecDown_pt'       
                    elif jec == 'met_uPerp_zll':
                        varTitle    = 'u_{#perp}'+"  "+ ' [GeV]'
                        vari = '((-met_pt*cos(met_phi)-zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi)-(-met_pt*sin(met_phi)-zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi))/zll_pt' 
                    elif jec == 'met_uPerp_zll_Up':
                        vari = '((-met_jecUp_pt*cos(met_jecUp_phi)-zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi)-(-met_jecUp_pt*sin(met_jecUp_phi)-zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi))/zll_pt' 
                    elif jec == 'met_uPerp_zll_Down':
                        vari = '((-met_jecDown_pt*cos(met_jecDown_phi)-zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi)-(-met_jecDown_pt*sin(met_jecDown_phi)-zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi))/zll_pt' 
                    elif jec == 'met_uPara_zll':
                        varTitle =  'u_{||} + Z p_{T} [GeV]'
                        vari = '((-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi)+(-met_pt*sin(met_phi)-zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi))/zll_pt+zll_pt'  
                    elif jec == 'met_uPara_zll_Up':
                        vari = '((-met_jecUp_pt*cos(met_jecUp_phi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi)+(-met_jecUp_pt*sin(met_jecUp_phi)-zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi))/zll_pt+zll_pt'  
                    elif jec == 'met_uPara_zll_Down':
                        vari = '((-met_jecDown_pt*cos(met_jecDown_phi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi)+(-met_jecDown_pt*sin(met_jecDown_phi)-zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi))/zll_pt+zll_pt'  
                    elif jec == 'met_uPerp_gamma':
                        varTitle    = 'u_{#perp}'+"  "+ ' [GeV]'
                        vari = '((-met_pt*cos(met_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi)-(-met_pt*sin(met_phi)-gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt' 
                    elif jec == 'met_uPerp_gamma_Up':
                        vari = '((-met_jecUp_pt*cos(met_jecUp_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi)-(-met_jecUp_pt*sin(met_jecUp_phi)-gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt' 
                    elif jec == 'met_uPerp_gamma_Down':
                        vari = '((-met_jecDown_pt*cos(met_jecDown_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi)-(-met_jecDown_pt*sin(met_jecDown_phi)-gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt' 
                    elif jec == 'met_uPara_gamma':
                        varTitle =  'u_{||} + Z p_{T} [GeV]'
                        vari = '((-met_pt*cos(met_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi)+(-met_pt*sin(met_phi)-gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt+gamma_pt'  
                    elif jec == 'met_uPara_gamma_Up':
                        vari = '((-met_jecUp_pt*cos(met_jecUp_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi)+(-met_jecUp_pt*sin(met_jecUp_phi)-gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt+gamma_pt'  
                    elif jec == 'met_uPara_gamma_Down':
                        vari = '((-met_jecDown_pt*cos(met_jecDown_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi)+(-met_jecDown_pt*sin(met_jecDown_phi)-gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt+gamma_pt'  
                    elif jec == 'nVert':
                        varTitle = 'nVert'
                        vari = 'nVert'  
                    elif jec == 'nVertUp':
                        vari = 'nVert'  
                    elif jec == 'nVertDown':
                        vari = 'nVert'  
                    elif jec == 'zll_pt':
                        varTitle    = 'Z p_{T} [GeV]'
                        vari = 'zll_pt'               
                    elif jec == 'zll_pt_Up':
                        vari = 'zll_pt'              
                    elif jec == 'zll_pt_Down':
                        vari = 'zll_pt'              
                    elif jec == 'zll_mass':
                        varTitle    = 'm_{ll} [GeV]'
                        vari = 'zll_mass'       
                    elif jec == 'gamma_pt':
                        varTitle    = '#gamma p_{T} [GeV]'
                        vari = 'gamma_pt'       
                    tmp_histo = 0
                    tmp_full = 0
                    for tree in ( ([treeDA] +mcTrees) if reg.doData else mcTrees):
                        ind = 0
                        cuts = CutManager.CutManager()
                        treename = tree.name.lower()
                        print '... doing tree %s' %treename
                        dataMC = ('DATA' if tree == treeDA else 'MC'); isData = dataMC == 'DATA';
                        block = tree.blocks[0]
                        attr = (var+'' if tree.name in ['DATA'] else jec+'')
                        tmp_full= tree.getTH1F(lumi, jec+"_"+eta+reg.name+treename, vari, reg.bins[reg.rvars.index(var)], 1, 1, reg.cuts , "", varTitle)
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
                        if isData: 
                            print "data", tmp_full.Integral()
                            continue
                        tmp_histo = copy.deepcopy(tmp_full.Clone(var+eta+reg.name))
                        tmp_histo.GetXaxis().SetTitle(varTitle)
                        if jecs.index(jec) == 0:
                            mc_stack.Add(tmp_histo)
                            if not ind: mc_stack.Draw()
                            ind+=1
                            mc_stack.GetXaxis().SetTitle(tmp_histo.GetXaxis().GetTitle())
                            print treename,  'histo has integral %.2f'%tmp_histo.Integral()
                            if not mc_histo:
                                mc_histo = copy.deepcopy(tmp_histo)
                            else:
                                mc_histo.Add(tmp_histo, 1.)   
                        if jecs.index(jec) == 1:
                            if not mc_up:                                    
                                mc_up = copy.deepcopy(tmp_histo)
                            else:
                                mc_up.Add(tmp_histo, 1.)                    
                        if jecs.index(jec) == 2:
                            if not mc_down:
                                mc_down = copy.deepcopy(tmp_histo)
                            else:
                                mc_down.Add(tmp_histo, 1.)                    
                    #mc_histo.SetMinimum(0);
                    setattr(reg, '%s_mc_histo_%s'%(jec, eta), mc_histo)
                    setattr(reg, '%s_mc_stack_%s'%(jec, eta), mc_stack)
                    setattr(reg, '%s_mc_up_%s'%(jec, eta), mc_up)
                    setattr(reg, '%s_mc_down_%s'%(jec, eta), mc_down)
                print 'plotting %s in region %s in %s' %(var, reg.name, eta)
                plot_var = Canvas.Canvas("test/%s/%s_%s%s"%(lumi_str, var,reg.name, channel), "png,root", 0.7, 0.6, 0.9, 0.9)
                plot_var.addStack(mc_stack  , "hist" , 1, 1)
                plot_var.addHisto(getattr(reg, var).getHisto('DATA', eta), "E,SAME"   , "Data"  , "PL", r.kBlack , 1, 0)
                plot_var.saveRatio(1, 1, 1, lumi, getattr(reg, var).getHisto('DATA', eta), mc_histo, mc_up, mc_down )
                time.sleep(0.1)

