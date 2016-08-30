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

def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm, upd,  cut, var, plot_name):
    f = 0
    efwhm = 0
    v_m.setVal(Hist.mean(x) );
    v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));
    voigt = RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w)
    xFrame = x.frame()
    Hist.plotOn(xFrame)
    if BKGSubtraction:
        plotda = Canvas.Canvas("met/gammajets/resolution/fit/"+dire+var+dependence+cut , "png",0.6, 0.7, 0.8, 0.9)
        bkg_pdf =  RooHistPdf("bkg_pdf","bkg_pdf",RooArgSet(x),bkg_hist)
        lAbkgFrac = RooRealVar("AbkgFrac","AbkgFrac",0.5,0.,1.)
        sigbkgFrac = RooFormulaVar("bkgfrac","@0",RooArgList(lAbkgFrac))
        model = RooAddPdf("modelSB","modelSB",voigt,bkg_pdf,sigbkgFrac)
        result = model.fitTo (Hist, RooFit.Minimizer("Minuit","Migrad"),RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel (-1)) 
        model.plotOn(xFrame)
        model.plotOn(xFrame, RooFit.Components("bkg_pdf"), RooFit.LineColor(r.kRed)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kRed)  ,RooFit.DrawOption("F"))
        model.plotOn(xFrame, RooFit.Components("voigt")  , RooFit.LineColor(r.kGreen)  ,RooFit.LineStyle(r.kDashed),RooFit.FillColor(r.kGreen)  ,RooFit.DrawOption("L")) 
        xFrame.Draw()
    else:   
        plotmc = Canvas.Canvas("met/gammajets/resolution/mcfit/"+dire+var+dependence+cut , "png",0.6, 0.7, 0.8, 0.9)
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
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
    if BKGSubtraction:
        plotda.save(0, 1, 0, xFrame.chiSquare(), "FWHM = %.1f " %(f), "x", "")
    else:
        if (upd ==""):
            plotmc.save(0, 1, 0, xFrame.chiSquare(),  "FWHM = %.1f " %(f), "x", "") 
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
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()
    print 'Going to load DATA and MC trees...'
    bkgDatasets = ['QCD_HT200to300', 'QCD_HT500to700',  'QCD_HT1000to1500', 'QCD_HT1500to2000', 'QCD_HT2000toInf', 'TGJets', 'TTGJets', 'ZGJets','ZGJets40to130', 'WGToLNuG', 'WJetsToLNu_HT100to200', 'WJetsToLNu_HT200to400', 'WJetsToLNu_HT400to600', 'WJetsToLNu_HT600to800', 'WJetsToLNu_HT800to1200', 'WJetsToLNu_HT1200to2500', 'WJetsToLNu_HT2500toInf']
    gjetsDatasets = ['GJets_HT40to100', 'GJets_HT100to200', 'GJets_HT200to400',  'GJets_HT600toInf']
    daDatasets = ['SinglePhoton_Run2016B_PromptReco_v2', 'SinglePhoton_Run2016C_PromptReco_v2', 'SinglePhoton_Run2016D_PromptReco_v2']
    channel = 'gamma'
    if doBKGSubtraction: 
        treeBKG = Sample.Tree(helper.selectSamples(opts.sampleFile, bkgDatasets, 'bkg'), 'bkg'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'da'), 'da', 1)
        trees = [treeBKG, treeDA] 
        plot_name = 'Data'
    else:
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets'  , 0)
        trees = [treeGJETS]
        plot_name =  'GJets'
    direction = [plot_name]  
    updown = [""]  
    uPerp = [makeUperp('met_pt', 'met_phi', 'gamma_pt', 'gamma_phi')] 
    uPara = [makeUpara('met_pt', 'met_phi', 'gamma_pt', 'gamma_phi')]
    metx = ['met_pt*sin(met_phi)']
    mety = ['met_pt*cos(met_phi)']
    uPerpPuppi = [makeUperp('metPuppi_pt', 'metPuppi_phi', 'gamma_pt', 'gamma_phi')] 
    uParaPuppi = [makeUpara('metPuppi_pt', 'metPuppi_phi', 'gamma_pt', 'gamma_phi')] 

    variable = [uPerp, uPara, uPerpPuppi, uParaPuppi]
    variablename = ['_uPerp', '_uPara', '_uPerpPuppi', '_uParaPuppi']
    dependences = [  'gamma_pt']
    print 'Trees successfully loaded...'
    lumi = 12.9

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    
    GammaJets = Region.region('met' +variablename[0], 
            [cuts.gammas()],
            'met_pt',
            'gamma_pt',
            variablename[0], 
            [range(-100,300,50)],
            False)
    regions.append(GammaJets)                                 
    
    vtxbins = [[0, 6],[6, 8],[8, 10], [10, 12],[12 ,14],[14, 16],[16, 20], [20,  40]]
    qtbins = [[0, 15], [15, 30],[30, 51],  [51, 65],[65, 80],[80,110], [110, 140], [140, 170], [170, 200],[200, 250] ,[250, 330]]
    sumetbins = [ [400, 500],[500, 600], [600, 700], [700, 800],[800, 900],  [900, 1000],[1000, 1100], [1100, 1300],  [1300, 1600], [1600, 2500]]
    scalebins = [ [18,26],[26, 34],[34, 42],[42, 50], [50, 58], [58, 66], [66, 74], [74, 82],[82,90],[90,100], [100, 120], [120, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305], [305, 335], [335, 365], [365, 395], [395, 430],  [430, 500]]

    cutsList = []
    binposition = []
    binerror = []
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
            f2 = TFile(dire+ "gjetsResolution.root", "UPDATE");                  
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
                        binLow = -80.; binHigh = 80.; nbins = 70;
                        x = RooRealVar("x", "x", binLow, binHigh)
                        for tree in trees:
                            if doBKGSubtraction:
                                if tree.name == 'bkg':
                                    bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,nbins, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                    bkg_Hist = RooDataHist("bkg","bkg"+i,RooArgList(x),bkg_hist)                 
                                    m_bkg = bkg_hist.GetMean()
                                    um_bkg = bkg_hist.GetMean()-bkg_hist.GetRMS()
                                    uM_bkg = bkg_hist.GetMean()+bkg_hist.GetRMS()
                                if tree.name == 'da':
                                    data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,nbins, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)
                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()  
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            if tree.name == 'gjets':
                                gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var, nbins, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)
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







