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



def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm, upd, cut, var, plot_name):
    f = 0
    efwhm = 0
    v_m.setVal(Hist.mean(x) );
    v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));
    
    voigt = RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w)
    xFrame = x.frame()
    Hist.plotOn(xFrame)
    if BKGSubtraction:
        plotda = Canvas.Canvas("met/zjets/resolution/fit/"+dire+var+dependence+eemm +cut , "png",0.6, 0.7, 0.8, 0.9)
        bkg_pdf =  RooHistPdf("bkg_pdf","bkg_pdf",RooArgSet(x),bkg_hist)
        lAbkgFrac = RooRealVar("AbkgFrac","AbkgFrac",0.5,0.,1.)
        sigbkgFrac = RooFormulaVar("bkgfrac","@0",RooArgList(lAbkgFrac))
        model = RooAddPdf("modelSB","modelSB",voigt,bkg_pdf,sigbkgFrac)
        result = model.fitTo (Hist, RooFit.Minimizer("Minuit","Migrad"),RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel (-1)) 
        model.plotOn(xFrame)
        model.plotOn(xFrame, RooFit.Components("bkg_pdf"), RooFit.LineColor(r.kRed)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kRed)  ,RooFit.DrawOption("F"))
        model.plotOn(xFrame, RooFit.Components("voigt")  , RooFit.LineColor(r.kGreen)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kGreen)  ,RooFit.DrawOption("L")) 
        xFrame.Draw()
        plotda.save(0, 1, 0, xFrame.chiSquare())
    else:   
        plotmc = Canvas.Canvas("met/zjets/resolution/mcfit/"+dire+var+dependence+eemm +cut , "png",0.6, 0.7, 0.8, 0.9)
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
        if (upd == ""):
            plotmc.save(0, 1, 0, xFrame.chiSquare())                                                                                                                                                      
    sigma = g_w.getVal()
    gamma = gamma_Z0.getVal()
    esigma = g_w.getError()
    egamma = gamma_Z0.getError()
    Vsg = result.correlation(g_w, gamma_Z0)
    Vgs = result.correlation(gamma_Z0, g_w)
    Vss = result.correlation(g_w, g_w)
    Vgg = result.correlation(gamma_Z0, gamma_Z0)
    f = FWHM(sigma, gamma)
    efwhm = FWHMError_fixed (sigma, gamma, esigma, egamma, Vss, Vsg, Vgs, Vgg)
    print 'fwhm', f
    return f, efwhm

def FWHMError_fixed(sigma, gamma, esigma, egamma, Vss, Vsg, Vgs, Vgg):
    a = 0.5346
    b = 0.2166
    c = 2 * math.sqrt( 2*math.log(2) )
    f_g = c * sigma
    f_l = 2 * gamma
    sq = math.sqrt( b * pow(f_l, 2) + pow(f_g, 2) )
    dfds = c * ( f_g / sq )
    dfdg = 2 * ( a + b * f_l / sq )
    p1 = dfds * dfds * esigma * esigma * pow( Vss, 2 )
    p2 = dfds * dfdg * esigma * egamma * pow( Vsg, 2 )
    p3 = dfdg * dfds * egamma * esigma * pow( Vgs, 2 )
    p4 = dfdg * dfdg * egamma * egamma * pow( Vgg, 2 );
    return math.sqrt ( p1 + p2 + p3 + p4 )

def FWHM(sigma, gamma):
    f_g = 2 * sigma * math.sqrt(2 * math.log(2))
    f_l = 2 * gamma
    return ((0.5346 * 2 * gamma + math.sqrt(0.2166 * f_l * f_l + f_g * f_g))/2.3548) 

if __name__ == "__main__":

    doBKGSubtraction = True
    doee = True
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
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
        destination = 'Not_BKG_Subtraction'
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
        direction = [plot_name]  
        updown = [""]  
        uPerp = ['((-met_pt*sin(met_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt'] 
        uPara = ['((-met_pt*sin(met_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt']
        uPerpPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_pt*cos(metPuppi_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt']
        uParaPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_pt*cos(metPuppi_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt'] 
    else:
        trees = [treeDY]
        #direction = [plot_name+'tgraphs']  
        direction = [plot_name,plot_name+'_up_jes_DY', plot_name+'_down_jes_DY']  
        #updown = [""]  
        updown = ["","_up", "_down"]  
        uPerp = ['((-met_pt*sin(met_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt', 
                '((-met_jecUp_pt*sin(met_jecUp_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_jecUp_pt*cos(met_jecUp_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt',
                '((-met_jecDown_pt*sin(met_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_jecDown_pt*cos(met_jecDown_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt']
        uPara = ['((-met_pt*sin(met_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_pt*cos(met_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt', 
                '((-met_jecUp_pt*sin(met_jecUp_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_jecUp_pt*cos(met_jecUp_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt',
                '((-met_jecDown_pt*sin(met_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_jecDown_pt*cos(met_jecDown_phi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt+ zll_pt']
        uPerpPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_pt*cos(metPuppi_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt', 
                      '((-metPuppi_jecUp_pt*sin(metPuppi_jecUp_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_jecUp_pt*cos(metPuppi_jecUp_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt',
                      '((-metPuppi_jecDown_pt*sin(metPuppi_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_jecDown_pt*cos(metPuppi_jecDown_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt']
        uParaPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_pt*cos(metPuppi_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt', 
                      '((-metPuppi_jecUp_pt*sin(metPuppi_jecUp_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_jecUp_pt*cos(metPuppi_jecUp_phi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt',
                      '((-metPuppi_jecDown_pt*sin(metPuppi_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_jecDown_pt*cos(metPuppi_jecDown_phi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt+ zll_pt']
        
        #updown = [""]  
        updown = ["","_up", "_down"]  
      #  uPerpRaw = ['((-met_rawPt*sin(met_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_rawPt*cos(met_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt', 
      #          '((-met_jecUp_rawPt*sin(met_jecUp_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_jecUp_rawPt*cos(met_jecUp_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt',
      #          '((-met_jecDown_rawPt*sin(met_jecDown_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-met_jecDown_rawPt*cos(met_jecDown_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt']
      #  uParaRaw = ['((-met_rawPt*sin(met_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_rawPt*cos(met_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt', 
      #          '((-met_jecUp_rawPt*sin(met_jecUp_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_jecUp_rawPt*cos(met_jecUp_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt',
      #          '((-met_jecDown_rawPt*sin(met_jecDown_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-met_jecDown_rawPt*cos(met_jecDown_rawPhi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt+ zll_pt']
      #  uPerpPuppiRaw = ['((-metPuppi_rawPt*sin(metPuppi_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_rawPt*cos(metPuppi_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt', 
      #                '((-metPuppi_jecUp_rawPt*sin(metPuppi_jecUp_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_jecUp_rawPt*cos(metPuppi_jecUp_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt',
      #                '((-metPuppi_jecDown_pt*sin(metPuppi_jecDown_phi)- zll_pt*sin(zll_phi))*zll_pt*cos(zll_phi) - (-metPuppi_jecDown_pt*cos(metPuppi_jecDown_phi) - zll_pt*cos(zll_phi))*zll_pt*sin(zll_phi))/zll_pt']
      #  uParaPuppiRaw = ['((-metPuppi_rawPt*sin(metPuppi_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_rawPt*cos(metPuppi_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt', 
      #                '((-metPuppi_jecUp_rawPt*sin(metPuppi_jecUp_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_jecUp_rawPt*cos(metPuppi_jecUp_rawPhi) - zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt + zll_pt',
      #                '((-metPuppi_jecDown_rawPt*sin(metPuppi_jecDown_rawPhi)- zll_pt*sin(zll_phi))*zll_pt*sin(zll_phi)+(-metPuppi_jecDown_rawPt*cos(metPuppi_jecDown_rawPhi)-zll_pt*cos(zll_phi))*zll_pt*cos(zll_phi))/zll_pt+ zll_pt']
    
    #variable = [uPara]
    variable = [uPerp, uPara]
    #variable = [uPerp, uPara, uPerpPuppi, uParaPuppi]
    #variablename = ['_uPara']
    variablename = ['_uPerp', '_uPara']
    #variablename = ['_uPerp', '_uPara', '_uPerpPuppi', '_uParaPuppi']
    #dependences = ['met_sumEt-zll_pt']
    #dependences = ['nVert']
    dependences = ['zll_pt']
    #dependences = ['zll_pt',   'nVert', 'met_sumEt-zll_pt']
    dependence = "zll_pt"
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
    
    vtxbins = [[0, 6],[6, 8],[8, 10], [10, 12],[12 ,14],[14, 16],[16, 20], [20,  40]]
    qtbins = [[0, 15], [15, 30],[30, 51],  [51, 65],[65, 80],[80,110], [110, 140], [140, 170], [170, 200],[200, 250] ,[250, 330]]
    sumetbins = [ [400, 500],[500, 600], [600, 700], [700, 800],[800, 900],  [900, 1000],[1000, 1100], [1100, 1300],  [1300, 1600], [1600, 2500]]

    cutsList = []
    binposition = []
    binerror = []
    if doBKGSubtraction:
        x = RooRealVar("x", "x", -90, 90)
    else: 
        x = RooRealVar("x", "x", -90, 90)

    g_w = RooRealVar("g_w", "width Gaus", 10.,0. , 100., "GeV") # sigma
    gamma_Z0 = RooRealVar("gamma_Z0_U", "Z0 width", 2.3, 0., 100., "GeV") # gamma
    v_m = RooRealVar("v_m", "v_m",0,-10.,10.)
    for dependence in dependences :                                                                                                                                                                  
        if dependence == 'nVert':
            bins = vtxbins        
        if dependence == 'zll_pt':
            bins = qtbins       
        if dependence == 'met_sumEt-zll_pt':
            bins = sumetbins           
        for dire in direction:
            f2 = TFile(dire+channel + "ZllResolution.root", "UPDATE");   
            upd = updown[direction.index(dire)]
            for vari in variable:
                var = vari[direction.index(dire)]
                for reg in regions:
                    cutsList = []
                    binposition = []
                    binerror = []
                    ratio = []
                    data_fwhms = []
                    data_fwhm = []
                    data_errors = []
                    data_error = []
                    mc_fwhm = []
                    mc_fwhms = []
                    mc_error = []
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
                        x = RooRealVar("x", "x", -100,100)
                        for tree in trees:
                            if doBKGSubtraction:
                                if tree.name == 'tt':
                                    tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var,65, -100., 100.,  i, '', variablename[variable.index(vari)]+i)    
                                    tt_Hist = RooDataHist("tt","tt"+i,RooArgList(x),tt_hist)                 
                                    m_tt = tt_hist.GetMean()
                                    um_tt = tt_hist.GetMean()-tt_hist.GetRMS()
                                    uM_tt = tt_hist.GetMean()+tt_hist.GetRMS()
                                if tree.name == 'da':
                                    data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,65, -100., 100.,  i, '', variablename[variable.index(vari)]+i)
                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()  
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            if tree.name == 'dy':
                                if doBKGSubtraction:
                                    dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,65, -100., 100.,  i, '', variablename[variable.index(vari)]+i)
                                else:
                                    if dependence == 'met_sumEt-zll_pt':
                                        if ((cutsList.index(i) == 2) or (cutsList.index(i) == 3)):
                                            dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,100, -60., 60.,  i, '', variablename[variable.index(vari)]+i)
                                        else:
                                            dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,100, -90., 90.,  i, '', variablename[variable.index(vari)]+i)
                                    else:
                                        dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var,100, -90., 90.,  i, '', variablename[variable.index(vari)]+i)
                                dy_Hist = RooDataHist("dy","dy"+i, RooArgList(x), dy_hist)
                                dy_BkgHist = RooDataHist()                                                          
                                m_dy = dy_hist.GetMean() 
                                um_dy = dy_hist.GetMean()-dy_hist.GetRMS()
                                uM_dy = dy_hist.GetMean()+dy_hist.GetRMS()
                        
                        if doBKGSubtraction:
                            data_fwhm, data_error = constructModel(data_Hist, tt_Hist, m_da, um_da, uM_da, True, channel, upd,  str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            data_fwhms.append(data_fwhm) 
                            data_errors.append(data_error)
                        else:
                            mc_fwhm, mc_error = constructModel(dy_Hist, dy_BkgHist, m_dy, um_dy, uM_dy, False, channel, upd, str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            mc_fwhms.append(mc_fwhm)
                            mc_errors.append(mc_error)
                    if doBKGSubtraction:
                        graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_fwhms), array("f", binerror), array("f", data_errors))
                        graph.Write (channel+"_met" +variablename[variable.index(vari)]+ "_vs_" + dependence);
                    else:
                        graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_fwhms), array("f", binerror), array("f", mc_errors))
                        graph_mc.Write (channel+"_met"+variablename[variable.index(vari)] + "_vs_" + dependence);
        f2.Close()







