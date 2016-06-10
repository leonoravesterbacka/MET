#####################################################################
#################### MET RESOLUTION##################################
#####################################################################
# this one plots scale of the different met types, in ee/mumu 


import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TLatex, TPaveStats, TGraph, TGraphErrors, TPad, TMultiGraph, TLine, TLegend, RooHistPdf, RooFormulaVar, RooPlot, RooAddPdf, RooFitResult,  RooDataHist, RooDataSet,  RooFit, RooArgList, RooRealVar, RooArgSet, RooLinkedList, RooVoigtian
import math, sys, optparse, copy, re, array
import numpy as np
from array import array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample


def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction,  channel, updo,  cut, var, plot_name):
    plot = Canvas.Canvas("met/gammajets/scale/fit/"+var+"_vs_gamma_pt_"+updo+cut , "png",0.6, 0.7, 0.8, 0.9)
    mean = 0.
    error = 0.
    v_m.setVal(Hist.mean(x) );
    print 'set mean to ', Hist.mean(x) ;
    v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));
    voigt = RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w)
    xFrame = x.frame()
    Hist.plotOn(xFrame)
    if BKGSubtraction:
        bkg_pdf =  RooHistPdf("bkg_pdf","bkg_pdf",RooArgSet(x),bkg_hist)
        lAbkgFrac = RooRealVar("AbkgFrac","AbkgFrac",0.5,0.,1.)
        sigbkgFrac = RooFormulaVar("bkgfrac","@0",RooArgList(lAbkgFrac))
        model = RooAddPdf("modelSB","modelSB",voigt,bkg_pdf,sigbkgFrac)
        result = model.fitTo (Hist, RooFit.Minimizer("Minuit","Migrad"),RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel (-1)) 
        model.plotOn(xFrame)
        model.plotOn(xFrame, RooFit.Components("bkg_pdf"), RooFit.LineColor(r.kRed)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kRed)  ,RooFit.DrawOption("F"))
        model.plotOn(xFrame, RooFit.Components("voigt")  , RooFit.LineColor(r.kGreen)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kGreen)  ,RooFit.DrawOption("L")) 
        xFrame.Draw()
        print 'chi2', xFrame.chiSquare()
        plot.save(0, 1, 0, xFrame.chiSquare())
    else:    
        plotmc = Canvas.Canvas("met/gammajets/scale/mcfit/"+var+"_vs_gamma_pt_"+ updo+cut , "png",0.6, 0.7, 0.8, 0.9)
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
        if (updo ==""):
            plotmc.save(0, 1, 0, xFrame.chiSquare())                                                                                                                                                      
    mean = v_m.getVal()  
    print 'after fit, mean is ', mean
    error = v_m.getError()  
    print 'error is ', error
    return mean, error

if __name__ == "__main__":

    doBKGSubtraction = True
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples2.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    ## make the options globa.. also the lumi
    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'

    bkgDatasets = [ 'QCD_HT200to300_ext', 'QCD_HT300to500', 'QCD_HT500to700','QCD_HT700to1000', 'QCD_HT1000to1500', 'QCD_HT1500to2000', 'QCD_HT2000toInf', 'TTGJets', 'ZGJets','ZGJets40-130', 'WGToLNuG', 'WJetsToLNu_HT200to400', 'WJetsToLNu_HT400to600', 'WJetsToLNu_HT600to800', 'WJetsToLNu_HT800to1200', 'WJetsToLNu_HT1200to2500', 'WJetsToLNu_HT2500toInf']
    gjetsDatasets = [  'ZGTo2LG', 'GJets_HT40to100', 'GJets_HT100to200', 'GJets_HT200to400', 'GJets_HT400to600', 'GJets_HT600toInf']
    daDatasets = ['SinglePhoton_Run2016B_PromptReco_v2NEW']
    channel = 'gamma'
    
    if doBKGSubtraction:
        treeBKG = Sample.Tree(helper.selectSamples(opts.sampleFile, bkgDatasets, 'bkg'), 'bkg'  , 0)
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'da'), 'da', 1)
        plot_name = 'Data'
    else:
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets', 0)
        trees = [treeGJETS] 
        plot_name = 'GJets'
   
    if doBKGSubtraction:
        trees = [treeBKG, treeGJETS, treeDA] 
        updown = [""]   
        direction = [plot_name+'tgraphs']  
        uPara = ['(((-met_pt*sin(met_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt'] 
        uParaPuppi = ['(((-metPuppi_pt*sin(metPuppi_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_pt*cos(metPuppi_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt'] 
    else:
        trees = [treeGJETS]
        direction = [plot_name+'tgraphs',plot_name+'_up_tgraphs_jes_GJets', plot_name+'_down_tgraphs_jes_GJets']  
        updown = ["","_up", "_down"]   
        uPara = ['(((-met_pt*sin(met_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt', 
            '(((-met_jecUp_pt*sin(met_jecUp_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_jecUp_pt*cos(met_jecUp_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt',
            '(((-met_jecDown_pt*sin(met_jecDown_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_jecDown_pt*cos(met_jecDown_phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt']
        uParaPuppi = ['(((-metPuppi_pt*sin(metPuppi_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_pt*cos(metPuppi_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt', 
                '(((-metPuppi_JetEnUp_Pt*sin(metPuppi_JetEnUp_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_JetEnUp_Pt*cos(metPuppi_JetEnUp_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt',
                '(((-metPuppi_JetEnDown_Pt*sin(metPuppi_JetEnDown_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_JetEnDown_Pt*cos(metPuppi_JetEnDown_Phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt)/gamma_pt']

    variable = [uPara]
    #variable = [uPara, uParaPuppi]
    variablename = [ '_uPara_']
    #variablename = [ '_uPara_', '_uParaPuppi_']
    dependence = 'gamma_pt'

    print 'Trees successfully loaded...'
    lumi = 0.803015796

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    
    region = Region.region('met' +variablename[0], 
            [cuts.gammas()],
            'met_uPerp_gamma',
            dependence,
            variablename[0], 
            [range(-100,300,50)],
            False)
    regions.append(region)

    qtbins = [ [20,28],[28, 36],[36, 44],[44, 52],[52,60],[60, 68], [68, 76], [76, 84],[84,92],[92,100], [100, 120], [120, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305]    , [305, 335], [335, 365], [365, 395], [395, 430],  [430, 500]]

    cutsList = []
    binposition = []
    binerror = []
    x = RooRealVar("x", "x", -2., 0.)
    g_w = RooRealVar("g_w", "width Gaus", 10.,0. , 100., "GeV") # sigma
    gamma_Z0 = RooRealVar("gamma_Z0_U", "Z0 width", 2.3, 0., 100., "GeV") # gamma
    v_m = RooRealVar("v_m", "v_m",0,-10.,10.)
   
    for reg in regions:
        cutsList = []
        binposition = []
        binerror = []
        for i in qtbins:
            mini = float(min(qtbins[qtbins.index(i)]))
            maxi = float(max(qtbins[qtbins.index(i)]))
            mid = float((mini+maxi)/2)
            binposition.append(mid)
            binerror.append(0.0)
            cutsList.append(cuts.AddList(reg.cuts +[  reg.dependence +  "<"+ str(maxi) + "&&" + reg.dependence + " >" + str(mini)] ))
        for dire in direction:
            f2 = TFile(dire + "gjetsScale.root", "UPDATE");   
            upd = updown[direction.index(dire)]
            for vari in variable:
                var = vari[direction.index(dire)]
                ratio = []
                data_means = []
                data_errors = []
                mc_means = []
                mc_errors = []
                print 'doing variable: ', dire , variablename[variable.index(vari)]
                for i in cutsList:
                     print i
                     for tree in trees:
                         if doBKGSubtraction:
                             if tree.name == 'bkg' :
                                 bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,100, -2., 0.,  i, '', variablename[variable.index(vari)]+i)    
                                 bkg_Hist = RooDataHist("bkg","bkg"+i,RooArgList(x),bkg_hist)                 
                                 m_bkg = bkg_hist.GetMean()
                                 um_bkg = bkg_hist.GetMean()-bkg_hist.GetRMS()
                                 uM_bkg = bkg_hist.GetMean()+bkg_hist.GetRMS()                                                                                                      
                             if tree.name == 'da':
                                 data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,100, -2., 0.,  i, '', variablename[variable.index(vari)]+i)
                                 data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                 m_da = data_hist.GetMean()
                                 um_da = data_hist.GetMean()-data_hist.GetRMS()
                                 uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                         if tree.name == 'gjets':
                             gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var,100, -2., 0.,  i, '', variablename[variable.index(vari)]+i)
                             gjets_Hist = RooDataHist("gjets","gjets"+i, RooArgList(x), gjets_hist)
                             gjets_BkgHist = RooDataHist()
                             m_gjets = gjets_hist.GetMean() 
                             um_gjets = gjets_hist.GetMean()-gjets_hist.GetRMS()
                             uM_gjets = gjets_hist.GetMean()+gjets_hist.GetRMS()
                     
                     if doBKGSubtraction:
                         data_mean, data_error = constructModel(data_Hist, bkg_Hist, m_da, um_da, uM_da, True, channel,  upd,  str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                         data_means.append(-data_mean)
                         data_errors.append(data_error)
                     else:
                        mc_mean, mc_error = constructModel(gjets_Hist, gjets_BkgHist, m_gjets, um_gjets, uM_gjets, False, channel, upd, str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                        mc_means.append(-mc_mean)
                        mc_errors.append(mc_error)                                                                                                                                                       

                if doBKGSubtraction:
                    graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_means), array("f", binerror), array("f", data_errors))
                    graph.Write ("met" +variablename[variable.index(vari)]+ "over_qt");
                else:
                    graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_means), array("f", binerror), array("f", mc_errors))
                    graph_mc.Write ("met"+variablename[variable.index(vari)] + "over_qt");
        f2.Close()







