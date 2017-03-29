#####################################################################
#################### MET SCALE ######################################
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
    uPara = ( "((( -"+ met + "*cos("+phi +" ) "+ " - " +boson + "*cos( " + boson_phi +"))* " + boson + "*cos( " +boson_phi +" )+(- "+ met + "*sin( "+ phi+ " )- " +boson+"*sin(" + boson_phi +" ))* " +boson +"*sin( " +boson_phi +" ))/" +boson+")/zll_pt" )
    return uPara

def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm,  updo, cut, var, plot_name, lumi):
    plotda = Canvas.Canvas("met/zjets/scale/fit/"+eemm+"_"+updo+"_"+cut , "png",0.6, 0.7, 0.8, 0.9)
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
        plotda.save(0, lumi, 0, xFrame.chiSquare(), "mean = %.2f"%(mean), "u_{||}/q_{T}", "")
    else:
        if (upd ==""):
            plotmc.save(0, lumi, 0, xFrame.chiSquare(),  "mean = %.2f"%(mean), "u_{||}/q_{T}", "") 
    return mean, error

if __name__ == "__main__":

    doBKGSubtraction = False
    doee = True
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()
    print 'Going to load DATA and MC trees...'

    if doBKGSubtraction:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo2L2Q_ee']
            tt = ['TTJets_DiLepton_total_ee' , 'WWTo2L2Nu_ee', 'TBarToLep_tch_ee', 'TBar_tWch_ee', 'T_tWch_ee', 'TToLep_sch_ee','TToLep_tch_ee',  'ZZTo2L2Nu_ee', 'WZTo3LNu_ee', 'WWW_ee', 'WWZ_ee', 'WZZ_ee', 'ZZZ_ee']
            da = ['DoubleEG_Run2016B_ReReco', 'DoubleEG_Run2016C_ReReco', 'DoubleEG_Run2016D_ReReco', 'DoubleEG_Run2016E_ReReco', 'DoubleEG_Run2016F_ReReco', 'DoubleEG_Run2016G_ReReco', 'DoubleEG_Run2016H_PromptReco_v2', 'DoubleEG_Run2016H_PromptReco_v3']  
            channel = 'E'
            lumi =  36.4 
        else:
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo2L2Q_mm' ]
            tt = ['TTJets_DiLepton_total_mm', 'WWTo2L2Nu_mm', 'TBarToLep_tch_mm', 'TBar_tWch_mm', 'T_tWch_mm', 'TToLep_sch_mm','TToLep_tch_mm',  'WWW_mm', 'WWZ_mm', 'WZZ_mm', 'ZZZ_mm', 'WZTo3LNu_mm', 'ZZTo2L2Nu_mm']
            #da = ['DoubleMuon_Run2016B_PromptReco_v2', 'DoubleMuon_Run2016C_PromptReco_v2', 'DoubleMuon_Run2016D_PromptReco_v2', 'DoubleMuon_Run2016E_PromptReco_v2', 'DoubleMuon_Run2016F_PromptReco_v2', 'DoubleMuon_Run2016G_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v3']  
            da = ['DoubleMuon_Run2016B_ReReco', 'DoubleMuon_Run2016C_ReReco', 'DoubleMuon_Run2016D_ReReco', 'DoubleMuon_Run2016E_ReReco', 'DoubleMuon_Run2016F_ReReco', 'DoubleMuon_Run2016G_ReReco',  'DoubleMuon_Run2016H_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v3']  
            channel = 'M'            
            lumi =  36.4
        plot_name = 'Data'
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, tt, 'tt'), 'tt', 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, da, 'da'), 'da', 1)
    else:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo2L2Q_ee']
            channel = 'E'                                                                 
            lumi =  36.4 
        else:                                                                             
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo2L2Q_mm' ]
            channel = 'M'                                                                 
            lumi = 36.4
        plot_name =  'DY'                                                                 
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)          
        
    if doBKGSubtraction:
        trees = [treeDY, treeTT, treeDA] 
        updown = [""]  
        direction = [plot_name] 
        uPara = [makeUpara('met_pt', 'met_phi', 'zll_pt', 'zll_phi')] 
        uParaPuppi = [makeUpara('metPuppi_pt', 'metPuppi_phi', 'zll_pt', 'zll_phi')] 
        metx = ['met_pt*sin(met_phi)']
        mety = ['met_pt*cos(met_phi)']
    else:    
        trees = [treeDY]
        times = 3
        
        
        direction = [plot_name,plot_name+'_up_jes_DY', plot_name+'_down_jes_DY', plot_name+'_up_jer_DY', plot_name+'_down_jer_DY',plot_name+'_up_uncl_DY', plot_name+'_down_uncl_DY']  
        
        updown = ["","_up", "_down",  "_jerUp", "_jerDown", "_unclUp", "_unclDown"]  
        uPara = [makeUpara('met_pt', 'met_phi', 'zll_pt', 'zll_phi'), makeUpara('met_jecUp_pt', 'met_jecUp_phi', 'zll_pt', 'zll_phi'), makeUpara('met_jecDown_pt', 'met_jecDown_phi', 'zll_pt', 'zll_phi'), makeUpara('met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUpara('met_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        uParaPuppi = [makeUpara('metPuppi_pt', 'metPuppi_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_jecUp_pt', 'metPuppi_jecUp_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_jecDown_pt', 'metPuppi_jecDown_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_shifted_JetResUp_pt', 'metPuppi_shifted_JetResUp_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_shifted_JetResDown_pt', 'metPuppi_shifted_JetResDown_phi', 'zll_pt', 'zll_phi'), makeUpara('metPuppi_shifted_UnclusteredEnUp_pt', 'metPuppi_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUpara('metPuppi_shifted_UnclusteredEnDown_pt', 'metPuppi_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        metx = ['met_pt*sin(met_phi)',  'met_jecUp_pt*sin(met_jecUp_phi)','met_jecDown_pt*sin(met_jecDown_phi)', 'met_shifted_UnclusteredEnUp_pt*sin(met_shifted_UnclusteredEnUp_phi)','met_shifted_UnclusteredEnDown_pt*sin(met_shifted_UnclusteredEnDown_phi)']
        mety = ['met_pt*cos(met_phi)',  'met_jecUp_pt*cos(met_jecUp_phi)','met_jecDown_pt*cos(met_jecDown_phi)', 'met_shifted_UnclusteredEnUp_pt*cos(met_shifted_UnclusteredEnUp_phi)','met_shifted_UnclusteredEnDown_pt*cos(met_shifted_UnclusteredEnDown_phi)']    
    variable = [uParaPuppi]
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
    qtbins = [ [18,24],[24, 30],[30, 38],[38, 46], [46, 52], [52, 60], [60, 68], [68, 76], [76, 84],[84,92],[92,100], [100, 115], [115, 130],  [130, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305], [305, 335], [335, 365], [365, 400], [400, 440],  [440, 500]]
    vtxbins = [[2, 4],[4 ,6],[6, 8], [8,  10],[10, 12],[12, 14],[14 ,16],[16,  18], [18, 20], [20,  22],[22,  24], [24, 26],[26,  28], [28,  30], [30, 32],[32,  34], [34 ,36],[36, 38], [38,  40] ]
    #vtxbins = [[3, 4], [4, 5],[5 ,6],[6, 7],[7, 8], [8,  9], [9, 10],[10, 11], [11, 12],[12, 13],[13, 14], [14, 15],[15 ,16],[16, 17],[17, 18], [18,  19], [19, 20], [20, 21], [21, 22],[22, 23],[23, 24], [24, 25],[25 ,26],[26, 27],[27, 28], [28,  29], [29, 30], [30, 31], [31, 32],[32, 33],[33, 34], [34, 35],[35 ,36],[36, 37],[37, 38], [38,  39], [39, 40], [40, 41], [41, 42],[42, 43], [43, 44], [44, 45] ]
    #vtxbins = [[0, 6],[6, 8],[8, 10], [10, 12],[12 ,14],[14, 16],[16, 20], [20,  40]]

    cutsList = []; binposition = []; binerror = []
    g_w = RooRealVar("g_w", "width Gaus", 10.,0. , 100., "GeV") # sigma
    gamma_Z0 = RooRealVar("gamma_Z0_U", "Z0 width", 2.3, 0., 100., "GeV") # gamma
    v_m = RooRealVar("v_m", "v_m",0,-2.,2.)
    for dependence in dependences :                                                                                                                                                                  
        if dependence == 'nVert':
            bins = vtxbins        
        if dependence == 'zll_pt':
            bins = qtbins       
        for dire in direction:
            f2 = TFile(dire + channel +"ZllScaleJan10ReRecoPuppi.root", "UPDATE");   # this is the file where all the tgraphs end up
            upd = updown[direction.index(dire)]
            for vari in variable:
                var = vari[direction.index(dire)]
                for reg in regions:
                    cutsList = [];binposition = [];binerror = [];ratio = [];data_means = [];data_errors = [];mc_means = [];mc_errors = []
                    print 'doing variable: ', dire , variablename[variable.index(vari)]
                    for i in bins:
                        mini = float(min(bins[bins.index(i)]))
                        maxi = float(max(bins[bins.index(i)]))
                        mid = float((mini+maxi)/2)
                        binposition.append(mid)
                        binerror.append(0.0)
                        
                        cutsList.append(cuts.AddList(reg.cuts + [  dependence +  "<"+ str(maxi) + "&&" + dependence + " >" + str(mini) ] ))
                    for i in cutsList:
                        print 'doing variable: ', dire , variablename[variable.index(vari)]
                        print i
                        
                        if ((variablename[variable.index(vari)] == '_uParaPuppi_') or (variablename[variable.index(vari)] == '_uParaRaw_')):
                            binLow = -1.7; binHigh = -0.3; nbins = 180;
                            if ((cutsList.index(i) < 12) ):
                                binLow = -2.3; binHigh = 0.3; nbins = 180;
                            if ((cutsList.index(i) < 7) ):
                                binLow = -3; binHigh = 1; nbins = 190;
                            if ((cutsList.index(i) < 5) ):
                                binLow = -2.5; binHigh = 1; nbins = 220;
                            if ((cutsList.index(i) ==0) ):
                                binLow = -2.5; binHigh = 1.; nbins = 320;
                        x = RooRealVar("x", "x", binLow, binHigh)
                        for tree in trees:
                            # so here loop over the trees, and get the TH1s for bkg, data or signal
                            if doBKGSubtraction:
                                if tree.name == 'tt':
                                    tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var, nbins, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                    tt_Hist = RooDataHist("tt","tt"+i,RooArgList(x),tt_hist)                 
                                    m_tt = tt_hist.GetMean()
                                    um_tt = tt_hist.GetMean()-tt_hist.GetRMS()
                                    uM_tt = tt_hist.GetMean()+tt_hist.GetRMS()                                                                                                      
                                if tree.name == 'da':
                                    data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,nbins, binLow, binHigh,  i , 'noOF', variablename[variable.index(vari)]+i)
                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            else:
                                dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,nbins*times, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)
                                dy_Hist = RooDataHist("dy","dy"+i, RooArgList(x), dy_hist)
                                dy_BkgHist = RooDataHist()
                                m_dy = dy_hist.GetMean()
                                um_dy = dy_hist.GetMean()-dy_hist.GetRMS()
                                uM_dy = dy_hist.GetMean()+dy_hist.GetRMS()
                        if doBKGSubtraction:
                            # and then get the fits here 
                            data_mean, data_error = constructModel(data_Hist, tt_Hist, m_da, um_da, uM_da, True, channel+variablename[variable.index(vari)]+dependence, upd,  str(cutsList.index(i)), dire, plot_name, lumi)
                            # fill some lists with all the fitted results
                            data_means.append(-data_mean)
                            data_errors.append(data_error)
                        else:
                            mc_mean, mc_error = constructModel(dy_Hist, dy_BkgHist, m_dy, um_dy, uM_dy, False, channel+variablename[variable.index(vari)] + dependence, upd, str(cutsList.index(i)), dire, plot_name, lumi)
                            mc_means.append(-mc_mean)
                            print "error ", mc_error
                            if mc_errors > 0.001:
                                mc_errors.append(0.0000001)      
                            else:
                                mc_errors.append(mc_error)                                                                                                                                                       
                    if doBKGSubtraction:
                        # make a little tgraph with the results
                        graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_means), array("f", binerror), array("f", data_errors))
                        graph.Write (channel + "_met" +variablename[variable.index(vari)]+ "over_"+dependence);
                        del tt_hist; del tt_Hist; del data_hist; del data_Hist ; del m_tt; del um_tt; del uM_tt; del m_da; del um_da ; del uM_da
                    else:
                        graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_means), array("f", binerror), array("f", mc_errors))
                        graph_mc.Write (channel + "_met"+variablename[variable.index(vari)] + "over_"+dependence);
                        del m_dy; del um_dy ; del uM_dy
                del nbins; del binLow; del binHigh; 
        f2.Close()                                                                                                                                                                                                                              
                                                                                                                                                                                                                                        
