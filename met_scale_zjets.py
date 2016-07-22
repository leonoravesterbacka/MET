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

def makeUpara(met, phi,boson, boson_phi):
    uPara = ( "((( -"+ met + "*cos("+phi +" ) "+ " - " +boson + "*cos( " + boson_phi +"))* " + boson + "*cos( " +boson_phi +" )+(- "+ met + "*sin( "+ phi+ " )- " +boson+"*sin(" + boson_phi +" ))* " +boson +"*sin( " +boson_phi +" ))/" +boson + " + " + boson +")" )
    return uPara
def makeUperp(met, phi,boson, boson_phi):
    uPerp = ( "((( -"+ met + "*cos("+phi +" ) "+ " - " +boson + "*cos( " + boson_phi +"))* " + boson + "*sin( " +boson_phi +" )-(- "+ met + "*sin( "+ phi+ " )- " +boson+"*sin(" + boson_phi +" ))* " +boson +"*cos( " +boson_phi +" ))/" +boson +")" )                             
    return uPerp

def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm,  updo, cut, var, plot_name):
    plot = Canvas.Canvas("met/zjets/scale/fit/"+eemm+"_"+updo+"_"+cut , "png",0.6, 0.7, 0.8, 0.9)
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
    else:    
        plotmc = Canvas.Canvas("met/zjets/scale/mcfit/"+eemm+"_vs_"+ updo+cut , "png",0.6, 0.7, 0.8, 0.9)
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
    mean = v_m.getVal()  
    print 'after fit, mean is ', mean
    error = v_m.getError()  
    if BKGSubtraction:
        plotda.save(0, 1, 0, xFrame.chiSquare(), "mean = %.1f"%(mean))
    else:
        if (upd ==""):
            plotmc.save(0, 1, 0, xFrame.chiSquare(),  "mean = %.1f"%(mean)) 
    return mean, error

if __name__ == "__main__":

    doBKGSubtraction = False
    doee = False
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()
    print 'Going to load DATA and MC trees...'

    if doBKGSubtraction:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo3LNu_ee', 'WZTo2L2Q_ee', 'ZZTo2L2Nu_ee', 'WWW_ee', 'WWZ_ee', 'WZZ_ee', 'ZZZ_ee']
            tt = ['TTJets_DiLepton_ee', 'WWTo2L2Nu_ee']
            da = ['DoubleEG_Run2016B_PromptReco_v2']
            channel = 'E'
        else:
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo3LNu_mm', 'WZTo2L2Q_mm', 'ZZTo2L2Nu_mm', 'WWW_mm', 'WWZ_mm', 'WZZ_mm', 'ZZZ_mm']
            tt = ['TTJets_DiLepton_mm', 'WWTo2L2Nu_mm']
            da = ['DoubleMuon_Run2016B_PromptReco_v2']
            channel = 'M'            
        plot_name = 'Data'
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, tt, 'tt'), 'tt', 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, da, 'da'), 'da', 1)
    else:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo3LNu_ee', 'WZTo2L2Q_ee', 'ZZTo2L2Nu_ee', 'WWW_ee', 'WWZ_ee', 'WZZ_ee', 'ZZZ_ee']
            channel = 'E'                                                                 
            lumi =  5.880219647738 
        else:                                                                             
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo3LNu_mm', 'WZTo2L2Q_mm', 'ZZTo2L2Nu_mm', 'WWW_mm', 'WWZ_mm', 'WZZ_mm', 'ZZZ_mm']
            channel = 'M'                                                                 
            lumi = 5.86459696193
        plot_name =  'DY'                                                                 
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)          
        
    if doBKGSubtraction:
        trees = [treeDY, treeTT, treeDA] 
        updown = [""]  
        direction = [plot_name] 
        uPerp = [makeUperp('met_pt', 'met_phi', 'zll_pt', 'zll_phi')] 
        uPara = [makeUpara('met_pt', 'met_phi', 'zll_pt', 'zll_phi')] 
        metx = ['met_pt*sin(met_phi)']
        mety = ['met_pt*cos(met_phi)']
    else:    
        trees = [treeDY]
        direction = [plot_name,plot_name+'_up_jes_DY', plot_name+'_down_jes_DY', plot_name+'_up_uncl_DY', plot_name+'_down_uncl_DY']  
        updown = ["","_up", "_down", "_unclUp", "_unclDown"]  
        uPerp = [makeUperp('met_pt', 'met_phi', 'zll_pt', 'zll_phi'), makeUperp('met_jecUp_pt', 'met_jecUp_phi', 'zll_pt', 'zll_phi'), makeUperp('met_jecDown_pt', 'met_jecDown_phi', 'zll_pt', 'zll_phi'), makeUperp('met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUperp('met_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        uPara = [makeUpara('met_pt', 'met_phi', 'zll_pt', 'zll_phi'), makeUpara('met_jecUp_pt', 'met_jecUp_phi', 'zll_pt', 'zll_phi'), makeUpara('met_jecDown_pt', 'met_jecDown_phi', 'zll_pt', 'zll_phi'), makeUpara('met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUpara('met_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        metx = ['met_pt*sin(met_phi)',  'met_jecUp_pt*sin(met_jecUp_phi)','met_jecDown_pt*sin(met_jecDown_phi)', 'met_shifted_UnclusteredEnUp_pt*sin(met_shifted_UnclusteredEnUp_phi)','met_shifted_UnclusteredEnDown_pt*sin(met_shifted_UnclusteredEnDown_phi)']
        mety = ['met_pt*cos(met_phi)',  'met_jecUp_pt*cos(met_jecUp_phi)','met_jecDown_pt*cos(met_jecDown_phi)', 'met_shifted_UnclusteredEnUp_pt*cos(met_shifted_UnclusteredEnUp_phi)','met_shifted_UnclusteredEnDown_pt*cos(met_shifted_UnclusteredEnDown_phi)']    
        uPerpPuppi = [makeUperp('metPuppi_pt', 'metPuppi_phi', 'zll_pt', 'zll_phi'), makeUperp('metPuppi_jecUp_pt', 'metPuppi_jecUp_phi', 'zll_pt', 'zll_phi'), makeUperp('metPuppi_jecDown_pt', 'metPuppi_jecDown_phi', 'zll_pt', 'zll_phi'), makeUperp('metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUperp('metPuppi_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        uParaPuppi = [makeUpara('metPuppi_pt', 'metPuppi_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_jecUp_pt', 'metPuppi_jecUp_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_jecDown_pt', 'metPuppi_jecDown_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUpara('metPuppi_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
    #variable = [metx, mety]
    #variable = [uPara]
    variable = [ uParaPuppi]
    #variablename = [ '_metx_', '_mety_']
    #variablename = ['_uPara_']
    variablename = [ '_uParaPuppi_']
    dependence = 'zll_pt'
    dependences = [ 'zll_pt']
    print 'Trees successfully loaded...'
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    region = Region.region('met' +variablename[0], 
            [cuts.leps()],
            'met_uPerp_zll',
            'zll_pt',
            variablename[0], 
            [range(-100,300,50)],
            False)
    regions.append(region)
    qtbins = [ [18,26],[26, 34],[34, 42],[42, 50], [50, 58], [58, 66], [66, 74], [74, 82],[82,90],[90,100], [100, 120], [120, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305], [305, 335], [335, 365], [365, 395], [395, 430],  [430, 500]]
    vtxbins = [[2, 4],[4 ,6],[6, 8], [8,  10],[10, 12],[12, 14],[14 ,16],[16,  18], [18, 20], [20,  22],[22,  24], [24, 26],[26,  28], [28,  30], [30, 32],[32,  34], [34 ,36],[36, 38], [38,  40] ]
    #vtxbins = [[3, 4], [4, 5],[5 ,6],[6, 7],[7, 8], [8,  9], [9, 10],[10, 11], [11, 12],[12, 13],[13, 14], [14, 15],[15 ,16],[16, 17],[17, 18], [18,  19], [19, 20], [20, 21], [21, 22],[22, 23],[23, 24], [24, 25],[25 ,26],[26, 27],[27, 28], [28,  29], [29, 30], [30, 31], [31, 32],[32, 33],[33, 34], [34, 35],[35 ,36],[36, 37],[37, 38], [38,  39], [39, 40], [40, 41], [41, 42],[42, 43], [43, 44], [44, 45] ]
    #vtxbins = [[0, 6],[6, 8],[8, 10], [10, 12],[12 ,14],[14, 16],[16, 20], [20,  40]]

    cutsList = []
    binposition = []
    binerror = []
    g_w = RooRealVar("g_w", "width Gaus", 10.,0. , 100., "GeV") # sigma
    gamma_Z0 = RooRealVar("gamma_Z0_U", "Z0 width", 2.3, 0., 100., "GeV") # gamma
    v_m = RooRealVar("v_m", "v_m",0,-10.,10.)
    for dependence in dependences :                                                                                                                                                                  
        if dependence == 'nVert':
            bins = vtxbins        
        if dependence == 'zll_pt':
            bins = qtbins       
        for dire in direction:
            f2 = TFile(dire + channel +"ZllScaleV1.root", "UPDATE");   
            upd = updown[direction.index(dire)]
            for vari in variable:
                var = vari[direction.index(dire)]
                for reg in regions:
                    cutsList = []
                    binposition = []
                    binerror = []
                    ratio = []
                    data_means = []
                    data_errors = []
                    mc_means = []
                    mc_errors = []
                    print 'doing variable: ', dire , variablename[variable.index(vari)]
                    for i in bins:
                        mini = float(min(bins[bins.index(i)]))
                        maxi = float(max(bins[bins.index(i)]))
                        mid = float((mini+maxi)/2)
                        binposition.append(mid)
                        binerror.append(0.0)
                        if ((dependence == 'nVert') or (dependence == 'met_sumEt-zll_pt')):
                            cutsList.append(cuts.AddList(reg.cuts + [  dependence +  "<"+ str(maxi) + "&&" + dependence + " >" + str(mini) + " &&  zll_pt > 50 "]))
                        else:
                            cutsList.append(cuts.AddList(reg.cuts + [  dependence +  "<"+ str(maxi) + "&&" + dependence + " >" + str(mini)]))
                    for i in cutsList:
                        print i
                        if (dependence == 'nVert'):
                            x = RooRealVar("x", "x", -100,100)
                        else:
                            x = RooRealVar("x", "x", -1.7,-0.3)
                        for tree in trees:
                            if doBKGSubtraction:
                                if tree.name == 'tt' :
                                    if (dependence == 'nVert'):
                                        tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var,50, -100, 100,  i, '', variablename[variable.index(vari)]+i)    
                                    else:
                                        if ((cutsList.index(i) < 5) ):
                                            x = RooRealVar("x", "x", -3.,1.)
                                            if (cutsList.index(i) == 0): 
                                                tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var,350, -3., 1.,  i, '', variablename[variable.index(vari)]+i)    
                                            else:
                                                tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var,190, -3., 1.,  i, '', variablename[variable.index(vari)]+i)    
                                        else:
                                            tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var,110, -1.7, -0.3,  i, '', variablename[variable.index(vari)]+i)    

                                    tt_Hist = RooDataHist("tt","tt"+i,RooArgList(x),tt_hist)                 
                                    m_tt = tt_hist.GetMean()
                                    um_tt = tt_hist.GetMean()-tt_hist.GetRMS()
                                    uM_tt = tt_hist.GetMean()+tt_hist.GetRMS()                                                                                                      
                                if tree.name == 'da':
                                    if (dependence == 'nVert'):
                                        data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,50, -100, 100,  i+"&&( HLT_DoubleMu ==1  || HLT_DoubleEG == 1 || HLT_SingleMu ==1 || HLT_SingleEle ==1 )" , '', variablename[variable.index(vari)]+i)
                                    else:
                                        if ((cutsList.index(i) < 5) ): 
                                            x = RooRealVar("x", "x", -3.,1.)
                                            if (cutsList.index(i) == 0) :
                                                data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,350, -3., 1.,  i+"&&( HLT_DoubleMu ==1 || HLT_DoubleEG ==1 || HLT_SingleMu ==1 || HLT_SingleEle ==1 )" , '', variablename[variable.index(vari)]+i)
                                            else:
                                                data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,190, -3., 1.,  i+"&&( HLT_DoubleMu ==1 || HLT_DoubleEG ==1 || HLT_SingleMu ==1 || HLT_SingleEle ==1 )" , '', variablename[variable.index(vari)]+i)
                                        else:
                                            data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,110, -1.7, -0.3,  i+"&&( HLT_DoubleMu ==1 || HLT_DoubleEG ==1 || HLT_SingleMu ==1 || HLT_SingleEle ==1 )" , '', variablename[variable.index(vari)]+i)


                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            if tree.name == 'dy':
                                if (dependence == 'zll_pt'): 
                                    if doBKGSubtraction:
                                        dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,110, -1.7, -0.3,  i, '', variablename[variable.index(vari)]+i)
                                    else: 
                                        if ((cutsList.index(i) < 5) ):
                                            x = RooRealVar("x", "x", -3.,1.)
                                            if (cutsList.index(i) == 0):
                                                dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,350, -3., 1.,  i, '', variablename[variable.index(vari)]+i)
                                            else:
                                                dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,190, -3., 1.,  i, '', variablename[variable.index(vari)]+i)
                                        else:
                                            dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,110, -1.7, -0.3,  i, '', variablename[variable.index(vari)]+i)
                                else:
                                    dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,50, -100, 100,  i, '', variablename[variable.index(vari)]+i)
                                dy_Hist = RooDataHist("dy","dy"+i, RooArgList(x), dy_hist)
                                dy_BkgHist = RooDataHist()
                                m_dy = dy_hist.GetMean()
                                um_dy = dy_hist.GetMean()-dy_hist.GetRMS()
                                uM_dy = dy_hist.GetMean()+dy_hist.GetRMS()
                        if doBKGSubtraction:
                            data_mean, data_error = constructModel(data_Hist, tt_Hist, m_da, um_da, uM_da, True, channel+variablename[variable.index(vari)]+dependence, upd,  str(cutsList.index(i)), dire, plot_name)
                            if (dependence == 'nVert'):
                                data_means.append(data_mean)
                            else: 
                                data_means.append(-data_mean)
                            data_errors.append(data_error)
                        else:
                            mc_mean, mc_error = constructModel(dy_Hist, dy_BkgHist, m_dy, um_dy, uM_dy, False, channel+variablename[variable.index(vari)] + dependence, upd, str(cutsList.index(i)), dire, plot_name)
                            if (dependence == 'nVert'):
                                mc_means.append(mc_mean)
                            else:
                                mc_means.append(-mc_mean)
                            mc_errors.append(mc_error)                                                                                                                                                       
                    if doBKGSubtraction:
                        graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_means), array("f", binerror), array("f", data_errors))
                        graph.Write (channel + "_met" +variablename[variable.index(vari)]+ "over_"+dependence);
                    else:
                        graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_means), array("f", binerror), array("f", mc_errors))
                        graph_mc.Write (channel + "_met"+variablename[variable.index(vari)] + "over_"+dependence);
        f2.Close()







