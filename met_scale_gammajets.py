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
    uPara = ( "(((( -"+ met + "*cos("+phi +" ) "+ " - " +boson + "*cos( " + boson_phi +"))* " + boson + "*cos( " +boson_phi +" )+(- "+ met + "*sin( "+ phi+ " )- " +boson+"*sin(" + boson_phi +" ))* " +boson +"*sin( " +boson_phi +" ))/" +boson + ") / " + boson +")" )
    return uPara

def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction,  channel, updo,  cut, var, plot_name):
    plotda = Canvas.Canvas("met/gammajets/scale/fit/"+var+channel+"_"+updo+cut , "png",0.6, 0.7, 0.8, 0.9)
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
        plotmc = Canvas.Canvas("met/gammajets/scale/mcfit/"+var+channel+"_"+ updo+cut , "png",0.6, 0.7, 0.8, 0.9)
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
    mean = v_m.getVal()  
    print 'after fit, mean is ', mean
    error = v_m.getError()  
    print 'error is ', error
    if BKGSubtraction:
        plotda.save(0, 1, 0, xFrame.chiSquare(), "mean = ", "u_{||}/q_{T}", "")
    else:
        if (upd ==""):
            plotmc.save(0, 1, 0, xFrame.chiSquare(),  "mean = ", "u_{||}/q_{T}", "") 
    return mean, error

if __name__ == "__main__":

    doBKGSubtraction = True
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args() 
    print 'Going to load DATA and MC trees...'
    bkgDatasets = ['QCD_HT200to300', 'QCD_HT300to500', 'QCD_HT500to700','QCD_HT700to1000',   'QCD_HT1000to1500', 'QCD_HT1500to2000', 'QCD_HT2000toInf', 'TGJets', 'TTGJets', 'ZGJets','ZGJets40to130', 'WGToLNuG', 'WJetsToLNu_HT100to200', 'WJetsToLNu_HT200to400', 'WJetsToLNu_HT400to600', 'WJetsToLNu_HT600to800', 'WJetsToLNu_HT800to1200', 'WJetsToLNu_HT1200to2500', 'WJetsToLNu_HT2500toInf']
    gjetsDatasets = ['GJets_HT40to100', 'GJets_HT100to200', 'GJets_HT200to400', 'GJets_HT400to600',  'GJets_HT600toInf', 'ZGTo2LG']
    daDatasets = ['SinglePhoton_Run2016B_PromptReco_v2', 'SinglePhoton_Run2016C_PromptReco_v2', 'SinglePhoton_Run2016D_PromptReco_v2']
    channel = 'gamma'
    if doBKGSubtraction:
        treeBKG = Sample.Tree(helper.selectSamples(opts.sampleFile, bkgDatasets, 'bkg'), 'bkg'  , 0)
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets'  , 0)
        treeDA = Sample.Tree(helper.selectSamples(opts.sampleFile, daDatasets, 'da'), 'da', 1)
        trees = [treeBKG, treeGJETS, treeDA] 
        plot_name = 'Data'
    else:
        treeGJETS = Sample.Tree(helper.selectSamples(opts.sampleFile, gjetsDatasets, 'gjets'), 'gjets', 0)
        trees = [treeGJETS] 
        plot_name = 'GJets'
    direction = [plot_name] 
    updown = [""] 
    uPara = [makeUpara('met_pt', 'met_phi', 'gamma_pt', 'gamma_phi')]
    metx = ['met_pt*sin(met_phi)']
    mety = ['met_pt*cos(met_phi)']
    uParaPuppi = [makeUpara('metPuppi_pt', 'metPuppi_phi', 'gamma_pt', 'gamma_phi')] 

    #variable = [metx, mety]
    variable = [uPara]
    variablename = [ '_uPara_']
    #variablename = [ '_metx_', '_mety_']
    dependences = ['gamma_pt' ]
    print 'Trees successfully loaded...'
    lumi = 12.9

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()

    regions = []
    
    region = Region.region('met' +variablename[0], 
            [cuts.gammas()],
            'met_uPerp_gamma',
            'gamma_pt',
            variablename[0], 
            [range(-100,300,50)],
            False)
    regions.append(region)
    qtbins = [ [18,26],[26, 34],[34, 42],[42, 50], [50, 58], [58, 66], [66, 74], [74, 82],[82,90],[90,100], [100, 120], [120, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 305], [305, 335], [335, 365], [365, 395], [395, 430],  [430, 500]]
    #vtxbins = [[2, 4],[4 ,6],[6, 8], [8,  10],[10, 12],[12, 14],[14 ,16],[16,  18], [18, 20], [20,  22],[22,  24], [24, 26],[26,  28], [28,  30], [30, 32],[32,  34], [34 ,36],[36, 38], [38,  40] ]
    vtxbins = [[1, 6],[6, 8],[8, 10], [10, 12],[12 ,14],[14, 16],[16, 20], [20,  40]]

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
        for dire in direction:
            f2 = TFile(dire + "gjetsScalewoOF.root", "UPDATE");   
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
                        cutsList.append(cuts.AddList(reg.cuts + [  dependence +  "<"+ str(maxi) + "&&" + dependence + " >" + str(mini)]))
                    
                    for i in cutsList:
                        print i
                        if (dependence == 'nVert'):
                            x = RooRealVar("x", "x", -100,100)
                        else:
                            x = RooRealVar("x", "x", -1.6,-0.4)
                        for tree in trees:
                            if doBKGSubtraction:
                                if tree.name == 'bkg' :
                                    if (dependence == 'nVert'):
                                        bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,100, -100, 100,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                    else:
                                        if ((cutsList.index(i) < 9) ):
                                            if ((cutsList.index(i) < 6) ):
                                                x = RooRealVar("x", "x", -2,0)
                                                bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,400, -2, 0,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                            else:
                                                bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,400, -1.6, -0.4,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                        else:
                                            bkg_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'bkg', var,200, -1.6, -0.4,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                    bkg_Hist = RooDataHist("bkg","bkg"+i,RooArgList(x),bkg_hist)                 
                                    m_bkg = bkg_hist.GetMean()
                                    um_bkg = bkg_hist.GetMean()-bkg_hist.GetRMS()
                                    uM_bkg = bkg_hist.GetMean()+bkg_hist.GetRMS()                                                                                                      
                                if tree.name == 'da':
                                    if (dependence == 'nVert'):
                                        data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,100, -100, 100,  i, 'noOF', variablename[variable.index(vari)]+i)
                                    else:
                                        if ((cutsList.index(i) < 9) ):
                                            if ((cutsList.index(i) < 6) ):
                                                x = RooRealVar("x", "x", -2,0)
                                                data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,400, -2, 0,  i, 'noOF', variablename[variable.index(vari)]+i)
                                            else:
                                                data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,400, -1.6, -0.4,  i, 'noOF', variablename[variable.index(vari)]+i)
                                        else:
                                            data_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'da', var,200, -1.6, -0.4,  i, 'noOF', variablename[variable.index(vari)]+i)
                                    data_Hist = RooDataHist("da","da"+i, RooArgList(x), data_hist)                  
                                    m_da = data_hist.GetMean()
                                    um_da = data_hist.GetMean()-data_hist.GetRMS()
                                    uM_da = data_hist.GetMean()+data_hist.GetRMS()                                                   
                            if tree.name == 'gjets':
                                if (dependence == 'nVert'):
                                    gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var,100, -100., 100.,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                else:
                                    if ((cutsList.index(i) < 9) ):
                                        if ((cutsList.index(i) < 6) ):
                                            x = RooRealVar("x", "x", -2, 0)
                                            gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var,400, -2, 0,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                        else:
                                            gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var,400, -1.6, -0.4,  i, 'noOF', variablename[variable.index(vari)]+i)    
                                    else: 
                                        gjets_hist = tree.getTH1F(lumi, variablename[variable.index(vari)]+i+'gjets', var,200, -1.6, -0.4,  i, 'noOF', variablename[variable.index(vari)]+i)
                                gjets_Hist = RooDataHist("gjets","gjets"+i, RooArgList(x), gjets_hist)
                                gjets_BkgHist = RooDataHist()
                                m_gjets = gjets_hist.GetMean() 
                                um_gjets = gjets_hist.GetMean()-gjets_hist.GetRMS()
                                uM_gjets = gjets_hist.GetMean()+gjets_hist.GetRMS()
                     
                        if doBKGSubtraction:
                            data_mean, data_error = constructModel(data_Hist, bkg_Hist, m_da, um_da, uM_da, True, channel,  upd,  str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            if (dependence == 'nVert'):
                                data_means.append(data_mean)
                            else:
                                data_means.append(-data_mean)
                            data_errors.append(data_error)
                        else:
                            mc_mean, mc_error = constructModel(gjets_Hist, gjets_BkgHist, m_gjets, um_gjets, uM_gjets, False, channel, upd, str(cutsList.index(i)), variablename[variable.index(vari)], plot_name)
                            if (dependence == 'nVert'):
                                mc_means.append(mc_mean)
                            else:
                                mc_means.append(-mc_mean)
                            mc_errors.append(mc_error)                                                                                                                                                       
                    if doBKGSubtraction:
                        graph = TGraphErrors(len(binposition), array("f", binposition), array("f", data_means), array("f", binerror), array("f", data_errors))
                        graph.Write ("met" +variablename[variable.index(vari)]+ "over_qt");
                    else:
                        graph_mc = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_means), array("f", binerror), array("f", mc_errors))
                        graph_mc.Write ("met"+variablename[variable.index(vari)] + "over_qt");
        f2.Close()







