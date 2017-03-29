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

def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm, upd, cut, var, plot_name):
    f = 0
    efwhm = 0
    v_m.setVal(Hist.mean(x) );
    v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));
    
    voigt = RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w)
    xFrame = x.frame()
    Hist.plotOn(xFrame)
    if BKGSubtraction:
        plotda = Canvas.Canvas("met/zjets/resolution/fit/"+dire+var+eemm +cut , "png",0.6, 0.7, 0.8, 0.9)
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
        plotmc = Canvas.Canvas("met/zjets/resolution/mcfit/"+dire+var+eemm +cut , "png",0.6, 0.7, 0.8, 0.9)
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

    doBKGSubtraction = False
    doee = False
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()
    
    print 'Going to load DATA and MC trees...'

    if doBKGSubtraction: 
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo2L2Q_ee']
            tt = ['TTJets_DiLepton_total_ee' , 'WZTo3LNu_ee', 'ZZTo2L2Nu_ee', 'WWTo2L2Nu_ee', 'TBarToLep_tch_ee', 'TBar_tWch_ee', 'T_tWch_ee', 'TToLep_sch_ee','TToLep_tch_ee', 'WWW_ee', 'WWZ_ee', 'WZZ_ee', 'ZZZ_ee']
            da = ['DoubleEG_Run2016B_ReReco', 'DoubleEG_Run2016C_ReReco', 'DoubleEG_Run2016D_ReReco', 'DoubleEG_Run2016E_ReReco', 'DoubleEG_Run2016F_ReReco', 'DoubleEG_Run2016G_ReReco', 'DoubleEG_Run2016H_PromptReco_v2', 'DoubleEG_Run2016H_PromptReco_v3']  
            channel = 'E'
            lumi =  36.4
        else:
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo2L2Q_mm']
            tt = ['TTJets_DiLepton_total_mm', 'WWTo2L2Nu_mm', 'WZTo3LNu_mm',  'ZZTo2L2Nu_mm', 'TBarToLep_tch_mm', 'TBar_tWch_mm', 'T_tWch_mm', 'TToLep_sch_mm','TToLep_tch_mm','WWW_mm', 'WWZ_mm', 'WZZ_mm', 'ZZZ_mm']
            #da = ['DoubleMuon_Run2016B_PromptReco_v2', 'DoubleMuon_Run2016C_PromptReco_v2', 'DoubleMuon_Run2016D_PromptReco_v2', 'DoubleMuon_Run2016E_PromptReco_v2', 'DoubleMuon_Run2016F_PromptReco_v2', 'DoubleMuon_Run2016G_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v3']  
            da = ['DoubleMuon_Run2016B_ReReco', 'DoubleMuon_Run2016C_ReReco', 'DoubleMuon_Run2016D_ReReco', 'DoubleMuon_Run2016E_ReReco', 'DoubleMuon_Run2016F_ReReco', 'DoubleMuon_Run2016G_ReReco', 'DoubleMuon_Run2016H_PromptReco_v2', 'DoubleMuon_Run2016H_PromptReco_v3']  
            channel = 'M'           
            lumi = 36.4
        plot_name = 'Data'
        destination = 'Not_BKG_Subtraction'
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)
        treeTT = Sample.Tree(helper.selectSamples(opts.sampleFile, tt, 'tt'), 'tt', 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, da, 'da'), 'da', 1)
    else:
        if doee:
            dy = ['DYJetsToLL_M50_ee', 'ZZTo4L_ee', 'ZZTo2L2Q_ee', 'WZTo2L2Q_ee']
            channel = 'E'
            lumi =  36.4
        else:
            dy = ['DYJetsToLL_M50_mm', 'ZZTo4L_mm', 'ZZTo2L2Q_mm', 'WZTo2L2Q_mm']
            channel = 'M'
            lumi = 36.4
        plot_name =  'DY'
        treeDY = Sample.Tree(helper.selectSamples(opts.sampleFile, dy, 'dy'), 'dy', 0)

    if doBKGSubtraction:
        trees = [treeTT, treeDA] 
        updown = [""]  
        direction = [plot_name] 
        uPerp = [makeUperp('met_pt', 'met_phi', 'zll_pt', 'zll_phi')] 
        uPara = [makeUpara('met_pt', 'met_phi', 'zll_pt', 'zll_phi')] 
        metx = ['met_pt*sin(met_phi)']                                                                                                      
        mety = ['met_pt*cos(met_phi)']                                                                                                      
        uPerpPuppi = [makeUperp('metPuppi_pt', 'metPuppi_phi', 'zll_pt', 'zll_phi')] 
        uParaPuppi = [makeUpara('metPuppi_pt', 'metPuppi_phi', 'zll_pt', 'zll_phi')] 
        uPerpPuppiRaw = [makeUperp('metPuppi_rawPt', 'metPuppi_rawPhi', 'zll_pt', 'zll_phi')] 
        uParaPuppiRaw = [makeUpara('metPuppi_rawPt', 'metPuppi_rawPhi', 'zll_pt', 'zll_phi')] 
        uPerpRaw = [makeUperp('met_rawPt', 'met_rawPhi', 'zll_pt', 'zll_phi')] 
        uParaRaw = [makeUpara('met_rawPt', 'met_rawPhi', 'zll_pt', 'zll_phi')] 
    else:
        trees = [treeDY]
        direction = [plot_name]  
        #direction = [plot_name,plot_name+'_up_jes_DY', plot_name+'_down_jes_DY', plot_name+'_up_jer_DY', plot_name+'_down_jer_DY',plot_name+'_up_uncl_DY', plot_name+'_down_uncl_DY']  
        updown = [""]  
        #updown = ["","_up", "_down",  "_jerUp", "_jerDown", "_unclUp", "_unclDown"]  
        uPara = [makeUpara('met_pt', 'met_phi', 'zll_pt', 'zll_phi'), makeUpara('met_jecUp_pt', 'met_jecUp_phi', 'zll_pt', 'zll_phi'), makeUpara('met_jecDown_pt', 'met_jecDown_phi', 'zll_pt', 'zll_phi'), makeUpara('met_shifted_JetResUp_pt', 'met_shifted_JetResUp_phi', 'zll_pt', 'zll_phi'), makeUpara('met_shifted_JetResDown_pt', 'met_shifted_JetResDown_phi', 'zll_pt', 'zll_phi'), makeUpara('met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUpara('met_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        uPerp = [makeUperp('met_pt', 'met_phi', 'zll_pt', 'zll_phi'), makeUperp('met_jecUp_pt', 'met_jecUp_phi', 'zll_pt', 'zll_phi'), makeUperp('met_jecDown_pt', 'met_jecDown_phi', 'zll_pt', 'zll_phi'), makeUperp('met_shifted_JetResUp_pt', 'met_shifted_JetResUp_phi', 'zll_pt', 'zll_phi'), makeUperp('met_shifted_JetResDown_pt', 'met_shifted_JetResDown_phi', 'zll_pt', 'zll_phi'), makeUperp('met_shifted_UnclusteredEnUp_pt', 'met_shifted_UnclusteredEnUp_phi', 'zll_pt', 'zll_phi'),  makeUperp('met_shifted_UnclusteredEnDown_pt', 'met_shifted_UnclusteredEnDown_phi', 'zll_pt', 'zll_phi')] 
        metx = ['met_pt*sin(met_phi)',  'met_jecUp_pt*sin(met_jecUp_phi)','met_jecDown_pt*sin(met_jecDown_phi)', 'met_shifted_UnclusteredEnUp_pt*sin(met_shifted_UnclusteredEnUp_phi)','met_shifted_UnclusteredEnDown_pt*sin(met_shifted_UnclusteredEnDown_phi)']
        mety = ['met_pt*cos(met_phi)',  'met_jecUp_pt*cos(met_jecUp_phi)','met_jecDown_pt*cos(met_jecDown_phi)', 'met_shifted_UnclusteredEnUp_pt*cos(met_shifted_UnclusteredEnUp_phi)','met_shifted_UnclusteredEnDown_pt*cos(met_shifted_UnclusteredEnDown_phi)']    

    #variable = [metx, mety]
    variable = [uPara , uPerp]
    #variable = [uPerp]
    #variable = [uPerp, uPara, uPerpPuppiRaw, uParaPuppiRaw ]
    #variable = [uPerpRaw, uParaRaw, uPerpPuppi,uParaPuppi]
    #variablename = ['_metx_', '_mety_']
    variablename = ['_uPara_', '_uPerp_' ]
    #variablename = ['_uPerp_']
    #variablename = ['_uPerp', '_uPara', '_uPerpPuppiRaw_', '_uParaPuppiRaw_' ]
    #variablename = ['_uPerpRaw', '_uParaRaw','_uPerpPuppi_', '_uParaPuppi_']
    dependences = ['nVert']
    #dependences = ['zll_pt']
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
    
    vtxbins = [[0.5, 6],[6, 8],[8, 10], [10, 12],[12,14], [14 ,16],[16, 18], [18, 20], [20,  22], [22, 24], [24, 26], [26, 28], [28, 32], [32, 40]]
    qtbins = [[0, 15], [15, 30],[30, 51],  [51, 65],[65, 80],[80,110], [110, 140], [140, 170], [170, 200],[200, 250] ,[250, 330]]
    sumetbins = [ [400, 500],[500, 600], [600, 700], [700, 800],[800, 900],  [900, 1000],[1000, 1100], [1100, 1300],  [1300, 1600], [1600, 2500]]
    scalebins = [ [0,24],[24, 30],[30, 38],[38, 46], [46, 52], [52, 60], [60, 68], [68, 76], [76, 84],[84,92],[92,100], [100, 115], [115, 130],  [130, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305], [305, 335], [335, 365], [365, 400], [400, 440],  [440, 500]]
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
        if dependence == 'met_sumEt-zll_pt':
            bins = sumetbins           
        for dire in direction:
            if doBKGSubtraction:
                fileD =  TFile("DataEZllScaleJan10ReReco.root");
            else: 
                fileD =  TFile("DYEZllScaleJan10ReReco.root");
            hScale = fileD.Get("E_met_uPara_over_zll_pt");
            a = "*( "
            for q in scalebins:   
                miniqt = float(min(scalebins[scalebins.index(q)]))
                maxiqt = float(max(scalebins[scalebins.index(q)]))
                midqt = float((miniqt+maxiqt)/2)
                x = hScale.GetY() 
                if scalebins.index(q)== 0:
                    a = a + " ("+ str(miniqt) +" < zll_pt && zll_pt < "+ str(maxiqt) +")*(1/ " + str(x[scalebins.index(q)]) +")" 
                else:
                    a = a + " + ("+ str(miniqt)+" < zll_pt && zll_pt <"+str(maxiqt)+")*(1/" + str(x[scalebins.index(q)]) +")"
            a = a+ ")"       
            print a
            fileD.Close()                                                                                                                
            f2 = TFile(dire+channel + "ZllResolutionJan10.root", "UPDATE");                                                           
            upd = updown[direction.index(dire)]
            for vari in variable:
                var = vari[direction.index(dire)]
                for reg in regions:
                    cutsList = [];binposition = [];binerror = [];ratio = [];data_fwhms = [];data_fwhm = [];data_errors = [];data_error = [];mc_fwhm = [];mc_fwhms = [];mc_error = [];mc_errors = []
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
                        binLow = -80.; binHigh = 80.; nbins = 250;
                        x = RooRealVar("x", "x", binLow, binHigh)
                        for tree in trees:
                            if doBKGSubtraction:
                                if tree.name == 'tt':
                                    tt_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'tt', var+a,nbins, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                    tt_Hist = RooDataHist("tt","tt"+i,RooArgList(x),tt_hist)                 
                                    m_tt = tt_hist.GetMean()
                                    um_tt = tt_hist.GetMean()-tt_hist.GetRMS()
                                    uM_tt = tt_hist.GetMean()+tt_hist.GetRMS()
                                if tree.name == 'da':
                                    data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var+a, nbins, binLow, binHigh,  i , 'noOF', variablename[variable.index(vari)]+i)
                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()  
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            else:
                                dy_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'dy', var+a, nbins, binLow, binHigh,  i, 'noOF', variablename[variable.index(vari)]+i)
                                dy_Hist = RooDataHist("dy","dy"+i, RooArgList(x), dy_hist)
                                dy_BkgHist = RooDataHist()                                                          
                                m_dy = dy_hist.GetMean() 
                                um_dy = dy_hist.GetMean()-dy_hist.GetRMS()
                                uM_dy = dy_hist.GetMean()+dy_hist.GetRMS()
                        
                        if doBKGSubtraction:
                            data_fwhm, data_error = constructModel(data_Hist, tt_Hist, m_da, um_da, uM_da, True, channel+dependence, upd,  str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            data_fwhms.append(data_fwhm) 
                            data_errors.append(data_error)
                        else:
                            mc_fwhm, mc_error = constructModel(dy_Hist, dy_BkgHist, m_dy, um_dy, uM_dy, False, channel+dependence, upd, str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            mc_fwhms.append(mc_fwhm)
                            mc_errors.append(mc_error)
                    if doBKGSubtraction:
                        graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_fwhms), array("f", binerror), array("f", data_errors))
                        graph.Write (channel+"_met" +variablename[variable.index(vari)]+ "_vs_" + dependence);
                    else:
                        graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_fwhms), array("f", binerror), array("f", mc_errors))
                        graph_mc.Write (channel+"_met"+variablename[variable.index(vari)] + "_vs_" + dependence);
        f2.Close()







