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



def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm, upd,  cut, var, plot_name):
    f = 0
    efwhm = 0
    v_m.setVal(Hist.mean(x) );
    v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));
    voigt = RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w)
    xFrame = x.frame()
    Hist.plotOn(xFrame)
    if BKGSubtraction:
        plot = Canvas.Canvas("met/gammajets/resolution/fit/"+dire+var+dependence+cut , "png",0.6, 0.7, 0.8, 0.9)
        bkg_pdf =  RooHistPdf("bkg_pdf","bkg_pdf",RooArgSet(x),bkg_hist)
        lAbkgFrac = RooRealVar("AbkgFrac","AbkgFrac",0.5,0.,1.)
        sigbkgFrac = RooFormulaVar("bkgfrac","@0",RooArgList(lAbkgFrac))
        model = RooAddPdf("modelSB","modelSB",voigt,bkg_pdf,sigbkgFrac)
        result = model.fitTo (Hist, RooFit.Minimizer("Minuit","Migrad"),RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel (-1)) 
        model.plotOn(xFrame)
        model.plotOn(xFrame, RooFit.Components("bkg_pdf"), RooFit.LineColor(r.kRed)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kRed)  ,RooFit.DrawOption("F"))
        model.plotOn(xFrame, RooFit.Components("voigt")  , RooFit.LineColor(r.kGreen)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kGreen)  ,RooFit.DrawOption("L")) 
        xFrame.Draw()
        plot.save(0, 1, 0, xFrame.chiSquare())
    else:   
        plotmc = Canvas.Canvas("met/gammajets/resolution/mcfit/"+dire+var+dependence+cut , "png",0.6, 0.7, 0.8, 0.9)
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
        if (upd ==""):
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
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples2.dat', help='the samples file. default \'samples2.dat\'')
    (opts, args) = parser.parse_args()

    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'
    channel = 'gamma'
    if doBKGSubtraction: 
        bkgDatasets = ['TTGJets',  'WGJets', 'ZGJets','WJetsToLNu_HT100to200_ext', 'WJetsToLNu_HT200to400', 'WJetsToLNu_HT400to600', 'WJetsToLNu_HT600to800', 'WJetsToLNu_HT800to1200', 'WJetsToLNu_HT1200to2500', 'WJetsToLNu_HT2500toInf']
        gjetsDatasets = ['QCD_HT200to300_ext', 'QCD_HT300to500', 'QCD_HT500to700','QCD_HT700to1000', 'QCD_HT1000to1500', 'QCD_HT1500to2000', 'QCD_HT2000toInf','GJets_HT40to100', 'GJets_HT100to200', 'GJets_HT200to400', 'GJets_HT400to600', 'GJets_HT600toInf']
        daDatasets = ['SinglePhoton_Run2016B_PromptReco_v2NEW'] 
        plot_name = 'Data'
        destination = 'Not_BKG_Subtraction'
        treeBKG = Sample.Tree(helper.selectSamples(opts.sampleFile, bkgDatasets, 'bkg'), 'bkg'  , 0)
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'da'), 'da', 1)
    else:
        gjetsDatasets = [ 'QCD_HT200to300_ext', 'QCD_HT300to500', 'QCD_HT500to700','QCD_HT700to1000', 'QCD_HT1000to1500', 'QCD_HT1500to2000', 'QCD_HT2000toInf','GJets_HT40to100', 'GJets_HT100to200', 'GJets_HT200to400', 'GJets_HT400to600', 'GJets_HT600toInf']
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets'  , 0)
        plot_name =  'GJets'
        destination = 'Not_BKG_Subtraction'

    if doBKGSubtraction:
        trees = [treeBKG, treeGJETS, treeDA] 
        direction = [plot_name+'tgraphs']  
        updown = [""]  
        uPerp = ['((-met_pt*sin(met_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt'] 
        uPara = ['((-met_pt*sin(met_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt + gamma_pt']
        uPerpPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-metPuppi_pt*cos(metPuppi_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt'] 
        uParaPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_pt*cos(metPuppi_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt + gamma_pt']
    else:      
        trees = [treeGJETS]
        direction = [plot_name+'tgraphs',plot_name+'_up_tgraphs_jes_GJets', plot_name+'_down_tgraphs_jes_GJets']  
        updown = ["","_up", "_down"]  
        uPerp = ['((-met_pt*sin(met_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt', 
            '((-met_JetEnUp_Pt*sin(met_JetEnUp_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-met_JetEnUp_Pt*cos(met_JetEnUp_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt',
            '((-met_JetEnDown_Pt*sin(met_JetEnDown_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-met_JetEnDown_Pt*cos(met_JetEnDown_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt']
        uPara = ['((-met_pt*sin(met_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_pt*cos(met_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt + gamma_pt', 
            '((-met_JetEnUp_Pt*sin(met_JetEnUp_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_JetEnUp_Pt*cos(met_JetEnUp_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt + gamma_pt',
            '((-met_JetEnDown_Pt*sin(met_JetEnDown_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-met_JetEnDown_Pt*cos(met_JetEnDown_Phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt+ gamma_pt']
        uPerpPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-metPuppi_pt*cos(metPuppi_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt',                                    
            '((-metPuppi_JetEnUp_Pt*sin(metPuppi_JetEnUp_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-metPuppi_JetEnUp_Pt*cos(metPuppi_JetEnUp_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt',
            '((-metPuppi_JetEnDown_Pt*sin(metPuppi_JetEnDown_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*cos(gamma_phi) - (-metPuppi_JetEnDown_Pt*cos(metPuppi_JetEnDown_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*sin(gamma_phi))/gamma_pt']
        uParaPuppi = ['((-metPuppi_pt*sin(metPuppi_phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_pt*cos(metPuppi_phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt + gamma_pt', 
            '((-metPuppi_JetEnUp_Pt*sin(metPuppi_JetEnUp_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_JetEnUp_Pt*cos(metPuppi_JetEnUp_Phi) - gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt + gamma_pt',
            '((-metPuppi_JetEnDown_Pt*sin(metPuppi_JetEnDown_Phi)- gamma_pt*sin(gamma_phi))*gamma_pt*sin(gamma_phi)+(-metPuppi_JetEnDown_Pt*cos(metPuppi_JetEnDown_Phi)-gamma_pt*cos(gamma_phi))*gamma_pt*cos(gamma_phi))/gamma_pt+ gamma_pt']

        
    variable = [uPerp, uPara]
    #variable = [uPerp, uPara, uPerpPuppi, uParaPuppi]
    variablename = ['_uPerp', '_uPara']
    #variablename = ['_uPerp', '_uPara', '_uPerpPuppi', '_uParaPuppi']
    dependences = ['gamma_pt', 'met_sumEt-gamma_pt', 'nVert']
    dependence = "gamma_pt"
    print 'Trees successfully loaded...'
    lumi = 0.6237

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    
    GammaJets = Region.region('met' +variablename[0], 
            [cuts.gammas()],
            'met_pt',
            dependence,
            variablename[0], 
            [range(-100,300,50)],
            False)
    regions.append(GammaJets)                                 
    
    vtxbins = [[0, 5],[5, 7],[7, 9], [9, 11],[11 ,13],[13, 15],[15, 18], [18, 25], [25, 40]]
    qtbins = [[0, 20], [20, 40],[40, 60],  [60, 80],[80, 100], [100, 120], [120, 140], [140, 175], [175, 225], [225, 320],[320,600]]
    sumetbins = [[300, 400], [400, 500],[500, 600], [600, 700], [700, 800],[800, 900],  [900, 1000],[1000, 1100], [1100, 1200], [1200, 1350], [1350, 1600], [1600, 2500]]

    cutsList = []
    binposition = []
    binerror = []
    x = RooRealVar("x", "x", -90, 90)

    g_w = RooRealVar("g_w", "width Gaus", 10.,0. , 100., "GeV") # sigma
    gamma_Z0 = RooRealVar("gamma_Z0_U", "Z0 width", 2.3, 0., 100., "GeV") # gamma
    v_m = RooRealVar("v_m", "v_m",0,-10.,10.)
    for dependence in dependences :                                                                                                                                                                  
        if dependence == 'nVert':
            bins = vtxbins        
        if dependence == 'gamma_pt':
            bins = qtbins       
        if dependence == 'met_sumEt-gamma_pt':
            bins = sumetbins           
        
        for dire in direction:
            f2 = TFile(destination+dire+channel + "gjets80XForResolution.root", "UPDATE");   
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
                        if ((dependence == 'nVert') or (dependence == 'met_sumEt')):
                            cutsList.append(cuts.AddList(reg.cuts + [  dependence +  "<"+ str(maxi) + "&&" + dependence + " >" + str(mini) + " &&  gamma_pt > 50 "]))
                        else:
                            cutsList.append(cuts.AddList(reg.cuts + [  dependence +  "<"+ str(maxi) + "&&" + dependence + " >" + str(mini)]))
                    for i in cutsList:
                        print i
                        for tree in trees:
                            if doBKGSubtraction:
                                if tree.name == 'bkg':
                                    bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,150, -90., 90.,  i, '', variablename[variable.index(vari)]+i)    
                                    bkg_Hist = RooDataHist("bkg","bkg"+i,RooArgList(x),bkg_hist)                 
                                    m_bkg = bkg_hist.GetMean()
                                    um_bkg = bkg_hist.GetMean()-bkg_hist.GetRMS()
                                    uM_bkg = bkg_hist.GetMean()+bkg_hist.GetRMS()
                                if tree.name == 'da':
                                    data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,150, -90., 90.,  i, '', variablename[variable.index(vari)]+i)
                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()  
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            if tree.name == 'gjets':
                                gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var,150, -90., 90.,  i, '', variablename[variable.index(vari)]+i)
                                gjets_Hist = RooDataHist("gjets","gjets"+i, RooArgList(x), gjets_hist)
                                gjets_BkgHist = RooDataHist()                                                          
                                m_gjets = gjets_hist.GetMean() 
                                um_gjets = gjets_hist.GetMean()-gjets_hist.GetRMS()
                                uM_gjets = gjets_hist.GetMean()+gjets_hist.GetRMS()
                        
                        if doBKGSubtraction:
                            if dependence == 'gamma_pt':
                                if cutsList.index(i) < 2:
                                    data_fwhm = 0
                                    data_error = 0      
                                else:                             
                                    data_fwhm, data_error = constructModel(data_Hist, bkg_Hist, m_da, um_da, uM_da, True, channel, upd,  str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            else:
                                data_fwhm, data_error = constructModel(data_Hist, bkg_Hist, m_da, um_da, uM_da, True, channel, upd,  str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            data_fwhms.append(data_fwhm)
                            data_errors.append(data_error)
                        else: 
                            if dependence == 'gamma_pt':
                                if cutsList.index(i) < 2:
                                    mc_fwhm = 0.
                                    mc_error = 0.                   
                                else:                             
                                    mc_fwhm, mc_error = constructModel(gjets_Hist, gjets_BkgHist, m_gjets, um_gjets, uM_gjets, False, channel, upd, str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            else:
                                mc_fwhm, mc_error = constructModel(gjets_Hist, gjets_BkgHist, m_gjets, um_gjets, uM_gjets, False, channel, upd, str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            mc_fwhms.append(mc_fwhm)
                            mc_errors.append(mc_error)
                    if doBKGSubtraction:
                        graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_fwhms), array("f", binerror), array("f", data_errors))
                        graph.Write (channel+"_met" +variablename[variable.index(vari)]+ "_vs_" + dependence);
                    else:
                        graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_fwhms), array("f", binerror), array("f", mc_errors))
                        graph_mc.Write (channel+"_met"+variablename[variable.index(vari)] + "_vs_" + dependence);
        f2.Close()







