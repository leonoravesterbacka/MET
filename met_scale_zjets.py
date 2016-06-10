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


def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm,  updo, cut, var, plot_name):
    plot = Canvas.Canvas("met/zjets/scale/fit/"+eemm+"uPara_vs_zll_pt_"+updo+"_"+cut , "png",0.6, 0.7, 0.8, 0.9)
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
        plotmc = Canvas.Canvas("met/zjets/scale/mcfit/"+eemm+"uPara_vs_zll_pt_"+ updo+cut , "png",0.6, 0.7, 0.8, 0.9)
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

    doBKGSubtraction = False
    doee = True
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples2.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    ## make the options globa.. also the lumi
    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'

    if doBKGSubtraction:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo3LNu_ee', 'WZTo2L2Q_ee']
            tt = ['TTJets_DiLepton_ee', 'WWTo2L2Nu_ee', 'WJetsToLNu_ee']
            da = ['DoubleEG_Run2016B_PromptReco_v2NEW']
            channel = 'E'
        else:
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo3LNu_mm', 'WZTo2L2Q_mm']
            tt = ['TTJets_DiLepton_mm', 'WWTo2L2Nu_mm', 'WJetsToLNu_mm']
            da = ['DoubleMuon_Run2016B_PromptReco_v2NEW']
            channel = 'M'            
        plot_name = 'Data'
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, tt, 'tt'), 'tt', 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, da, 'da'), 'da', 1)
    else:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo3LNu_ee', 'WZTo2L2Q_ee']
            channel = 'E'                                                                 
        else:                                                                             
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo3LNu_mm', 'WZTo2L2Q_mm']
            channel = 'M'                                                                 
        plot_name =  'DY'                                                                 
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)          
        
    if doBKGSubtraction:
        trees = [treeDY, treeTT, treeDA] 
        updown = [""]  
        direction = [plot_name+'tgraphs']  
        uPara = ['(((-met_pt*sin(met_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt'] 
        uParaPuppi = ['(((-metPuppi_pt*sin(metPuppi_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_pt*cos(metPuppi_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt'] 
    else:    
        trees = [treeDY]
        direction = [plot_name+'tgraphs',plot_name+'_up_tgraphs_jes_DY', plot_name+'_down_tgraphs_jes_DY']  
        updown = ["","_up", "_down"]  
        uPara = ['(((-met_pt*sin(met_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt', 
                '(((-met_jecUp_pt*sin(met_jecUp_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_jecUp_pt*cos(met_jecUp_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt',
               '(((-met_jecDown_pt*sin(met_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_jecDown_pt*cos(met_jecDown_phi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt']
        uParaPuppi = ['(((-metPuppi_pt*sin(metPuppi_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_pt*cos(metPuppi_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt', 
                '(((-metPuppi_jecUp_pt*sin(metPuppi_jecUp_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_jecUp_pt*cos(metPuppi_jecUp_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt',
                '(((-metPuppi_jecDown_pt*sin(metPuppi_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_jecDown_pt*cos(metPuppi_jecDown_phi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt)/zll_pt']
    
    variable = [uPara]
    #variable = [uPara, uParaPuppi]
    variablename = [ '_uPara_']
    #variablename = [ '_uPara_', '_uParaPuppi_']
    dependence = 'zll_pt'

    print 'Trees successfully loaded...'
    lumi = 0.803015796

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    
    region = Region.region('met' +variablename[0], 
            [cuts.leps()],
            'met_uPerp_zll',
            dependence,
            variablename[0], 
            [range(-100,300,50)],
            False)
    regions.append(region)
    
    qtbins = [ [20,28],[28, 36],[36, 44],[44, 52],[52,60],[60, 68], [68, 76], [76, 84],[84,92],[92,100], [100, 120], [120, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305], [305, 335], [335, 365], [365, 395], [395, 430],  [430, 500]]

    cutsList = []
    binposition = []
    binerror = []

   
    x = RooRealVar("x", "x", -1.8, -0.2)
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
            f2 = TFile(dire + channel +"ZllScale.root", "UPDATE");   
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
                             if tree.name == 'tt' :
                                 tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var,100, -1.8, -0.2,  i, '', variablename[variable.index(vari)]+i)    
                                 tt_Hist = RooDataHist("tt","tt"+i,RooArgList(x),tt_hist)                 
                                 m_tt = tt_hist.GetMean()
                                 um_tt = tt_hist.GetMean()-tt_hist.GetRMS()
                                 uM_tt = tt_hist.GetMean()+tt_hist.GetRMS()                                                                                                      
                             if tree.name == 'da':
                                 data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,100, -1.8, -0.2,  i, '', variablename[variable.index(vari)]+i)
                                 data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                 m_da = data_hist.GetMean()
                                 um_da = data_hist.GetMean()-data_hist.GetRMS()
                                 uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                         if tree.name == 'dy':
                             dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,100, -1.8, -0.2,  i, '', variablename[variable.index(vari)]+i)
                             dy_Hist = RooDataHist("dy","dy"+i, RooArgList(x), dy_hist)
                             dy_BkgHist = RooDataHist()
                             m_dy = dy_hist.GetMean()
                             um_dy = dy_hist.GetMean()-dy_hist.GetRMS()
                             uM_dy = dy_hist.GetMean()+dy_hist.GetRMS()
                     
                     if doBKGSubtraction:
                         data_mean, data_error = constructModel(data_Hist, tt_Hist, m_da, um_da, uM_da, True, channel, upd,  str(cutsList.index(i)), dire, plot_name)
                         data_means.append(-data_mean)
                         data_errors.append(data_error)
                     else:
                        mc_mean, mc_error = constructModel(dy_Hist, dy_BkgHist, m_dy, um_dy, uM_dy, False, channel, upd, str(cutsList.index(i)), dire, plot_name)
                        mc_means.append(-mc_mean)
                        mc_errors.append(mc_error)                                                                                                                                                       


                if doBKGSubtraction:
                    graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_means), array("f", binerror), array("f", data_errors))
                    graph.Write (channel + "_met" +variablename[variable.index(vari)]+ "over_qt");
                else:
                    graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_means), array("f", binerror), array("f", mc_errors))
                    graph_mc.Write (channel + "_met"+variablename[variable.index(vari)] + "over_qt");
        f2.Close()







