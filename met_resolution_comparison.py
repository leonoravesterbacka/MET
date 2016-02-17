#####################################################################
#################### MET RESOLUTION##################################
#####################################################################
# this one plots scale of the different met types, in ee/mumu 


import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TGraph, TGraphErrors, TPad, TMultiGraph, TLine, TLegend, RooHistPdf, RooFormulaVar, RooPlot, RooAddPdf, RooFitResult,  RooDataHist, RooDataSet,  RooFit, RooArgList, RooRealVar, RooArgSet, RooLinkedList, RooVoigtian
import math, sys, optparse, copy, re, array
import numpy as np
from array import array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample

def constructModel(Hist, bkg_hist,  m, um,uM, BKGSubtraction, eemm, cut, var, plot_name):
    print 'here now'
    plot = Canvas.Canvas("met/resolution_plots/"+var+"_vs_"+str(reg.dependence)+"_"+eemm+"_"+plot_name+"_"+cut , "png",0.6, 0.7, 0.8, 0.9)
    f = 0
    efwhm = 0
    v_m.setVal(m)
    v_m.setRange(um,uM) 
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
        #plot.save(0, 1, 0, lumi)
    else:    
        result = voigt.fitTo(Hist, RooFit.Minimizer("Minuit","Migrad"), RooFit.Strategy(2), RooFit.SumW2Error(False), RooFit.Save(True), RooFit.PrintLevel(-1)) 
        voigt.plotOn(xFrame,RooFit.FillColor(r.kGray),RooFit.VisualizeError(result,1),RooFit.Components("voigt"))
        voigt.plotOn(xFrame,RooFit.LineColor(r.kGray))  
        xFrame.Draw()
        #plot.save(0, 1, 0, lumi)
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
    #doBKGSubtraction = False
    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    ## make the options globa.. also the lumi
    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'
    dy_ee = ['DYJetsToLL_ee']
    dy_mm = ['DYJetsToLL_mm']
    tt_ee = ['TTJets_ee']
    tt_mm = ['TTJets_mm']
    data_ee = ['Data_ee']
    data_mm = ['Data_mm']


    treeDY_ee = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_ee, 'dy_ee'), 'dy_ee', 0)
    treeDY_mm = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_mm, 'dy_mm'), 'dy_mm', 0)
    treeTT_ee = Sample.Tree(helper.selectSamples(opts.sampleFile, tt_ee, 'tt_ee'), 'tt_ee', 0)
    treeTT_mm = Sample.Tree(helper.selectSamples(opts.sampleFile, tt_mm, 'tt_mm'), 'tt_mm', 0)
    treeData_ee = Sample.Tree(helper.selectSamples(opts.sampleFile, data_ee, 'da_ee'), 'da_ee', 1)
    treeData_mm = Sample.Tree(helper.selectSamples(opts.sampleFile, data_mm, 'da_mm'), 'da_mm', 1)
    
    if doBKGSubtraction: 
        #trees = [ treeData_mm, treeData_ee, treeTT_mm, treeTT_ee]
        plot_name = 'data'
    else:    
        #trees = [ treeDY_mm, treeDY_ee]
        plot_name = 'MC'
    print 'Trees successfully loaded...'
    lumi = 2.260

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    regions = []
    vector = []
    region = Region.region('Type 1 ME_{T}', 
            [cuts.GoodLepton()],
            'met_uPerp_zll',
            'zll_pt',
            [range(-100,300,50)],
            False)
    regions.append(region)                 
    region = Region.region('Type 1 ME_{T}',
            [cuts.GoodLepton()],
            'met_uPara_zll+zll_pt',
            'zll_pt',
            [range(-100,300,50)],
            False)
    regions.append(region)                  
    region = Region.region('Raw ME_{T}', 
            [cuts.GoodLepton()],
            'met_raw_uPerp_zll',
            'zll_pt',
            [range(-100,300,50)],
            False)
    regions.append(region)                 
    region = Region.region('Raw ME_{T}',
            [cuts.GoodLepton()],
            'met_raw_uPara_zll+zll_pt',
            'zll_pt',
            [range(-100,300,50)],
            False)
    regions.append(region)  
    ############sumEt
#    region = Region.region('Type 1 ME_{T}',
#            [cuts.GoodLepton()],
#            'met_uPerp_zll',
#            'met_sumEt',
#            [range(-100,300,50)],
#            False)
#    regions.append(region)                 
#    region = Region.region('Type 1 ME_{T}',
#            [cuts.GoodLepton()],
#            'met_uPara_zll+zll_pt',
#            'met_sumEt',
#            [range(-100,300,50)],
#            False)
#    regions.append(region)                 
#    region = Region.region('Raw 1 ME_{T}', 
#            [cuts.GoodLepton()],
#            'met_uPerp_zll',
#            'met_sumEt',
#            [range(-100,300,50)],
#            False)
#    regions.append(region)                 
#    region = Region.region('Raw 1 ME_{T}',
#            [cuts.GoodLepton()],
#            'met_uPara_zll+zll_pt',
#            'met_sumEt',
#            [range(-100,300,50)],
#            False)
#    regions.append(region)      

    qtbins = [ [25, 50], [50, 75]]
    #qtbins = [ [0, 10], [10, 20], [20, 30], [30, 40], [40, 60], [60, 80], [80, 100], [100, 130], [130, 160], [160, 200],[200, 300], [300,450]]
    metbins = [[0, 5], [5, 10], [10, 15], [15, 20], [20, 25], [25, 30], [30, 35], [35, 40], [40, 45], [45, 50]]
    sumetbins = [[100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000],[1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 1900], [1900, 2000],[2000, 2100], [2100, 2200]]


    graph_ee = TGraphErrors()
    graph_mm = TGraphErrors()
    cutsList = []
    binposition = []
    binerror = []
    

    x = RooRealVar("x", "x", -200, 200)
    g_w = RooRealVar("g_w", "width Gaus", 10.,0. , 100., "GeV") # sigma
    gamma_Z0 = RooRealVar("gamma_Z0_U", "Z0 width", 2.3, 0., 100., "GeV") # gamma
    v_m = RooRealVar("v_m", "v_m",0,-10.,10.)
    trees = [treeData_mm, treeData_ee, treeTT_mm, treeTT_ee, treeDY_ee, treeDY_mm]

    for reg in regions:
        if reg.name == 'Raw ME_{T}':
            regname = 'met_raw'       
        if reg.name == 'Type 1 ME_{T}':
            regname = 'met_type1'       
        if reg.rvars == 'met_raw_uPara_zll+zll_pt' or reg.rvars == 'met_uPara_zll+zll_pt':
            label = 'u_{|| }'
            var = 'uPara'
        if reg.rvars == 'met_raw_uPerp_zll' or reg.rvars == 'met_uPerp_zll'  :    
            var = 'uPerp'
            label = 'u_{#perp}'
        cutsList = []
        binposition = []
        binerror = []
        ratio_ee = []
        ratio_mm = []
        data_fwhms_ee = []
        data_fwhms_mm = []
        data_error_ee = []
        data_error_mm = []
        mc_fwhms_ee = []
        mc_fwhms_mm = []
        mc_error_mm = []
        mc_error_ee = []

        for i in qtbins:
            mini = float(min(qtbins[qtbins.index(i)]))
            maxi = float(max(qtbins[qtbins.index(i)]))
            mid = float((mini+maxi)/2)
            binposition.append(mid)
            binerror.append(0.0)
            cutsList.append(cuts.AddList([  reg.dependence +  "<"+ str(maxi) + "&&" + reg.dependence + " >" + str(mini)]))
        for i in cutsList:
            print i
            for tree in trees:
                print 'trees', tree.name
                if tree.name == 'dy_ee':                   
                    dy_hist_ee = tree.getTH1F(lumi, reg.name+i, reg.rvars,40, -100., 100.,  i, '', reg.name)
                    dy_Hist_ee = RooDataHist("dy_ee","dy_ee"+i, RooArgList(x), dy_hist_ee)
                    dy_BkgHist_ee = RooDataHist()                                                          
                    m_dy_ee = dy_hist_ee.GetMean()                            
                    um_dy_ee = dy_hist_ee.GetMean()-dy_hist_ee.GetRMS()
                    uM_dy_ee = dy_hist_ee.GetMean()+dy_hist_ee.GetRMS()
                elif tree.name == 'dy_mm':    
                    dy_hist_mm = tree.getTH1F(lumi, reg.name+i, reg.rvars,40, -100., 100.,  i, '', reg.name)
                    dy_Hist_mm = RooDataHist("dy_mm","dy_mm"+i, RooArgList(x), dy_hist_mm)
                    dy_BkgHist_mm = RooDataHist()                                                          
                    m_dy_mm = dy_hist_mm.GetMean()                            
                    um_dy_mm = dy_hist_mm.GetMean()-dy_hist_mm.GetRMS()
                    uM_dy_mm = dy_hist_mm.GetMean()+dy_hist_mm.GetRMS()
                elif tree.name == 'tt_ee':  
                    tt_hist_ee = tree.getTH1F(lumi, reg.name+i, reg.rvars,50, -150., 150.,  i, '', reg.name)
                    tt_Hist_ee = RooDataHist("tt_ee","tt_ee"+i,RooArgList(x),tt_hist_ee)               
                    m_tt_ee = tt_hist_ee.GetMean()                       
                    um_tt_ee = tt_hist_ee.GetMean()-tt_hist_ee.GetRMS()
                    uM_tt_ee = tt_hist_ee.GetMean()+tt_hist_ee.GetRMS()
                elif tree.name == 'tt_mm': 
                    tt_hist_mm = tree.getTH1F(lumi, reg.name+i, reg.rvars,50, -150., 150.,  i, '', reg.name)    
                    tt_Hist_mm = RooDataHist("tt_mm","tt_mm"+i,RooArgList(x),tt_hist_mm)                 
                    m_tt_mm = tt_hist_mm.GetMean()                       
                    um_tt_mm = tt_hist_mm.GetMean()-tt_hist_mm.GetRMS()
                    uM_tt_mm = tt_hist_mm.GetMean()+tt_hist_mm.GetRMS()
                elif tree.name== 'da_ee':
                    data_hist_ee = tree.getTH1F(lumi, reg.name+i, reg.rvars,50, -150., 150.,  i, '', reg.name)
                    data_Hist_ee = RooDataHist("da_ee","da_ee"+i, RooArgList(x), data_hist_ee)                   
                    m_da_ee = data_hist_ee.GetMean()                       
                    um_da_ee = data_hist_ee.GetMean()-data_hist_ee.GetRMS()
                    uM_da_ee = data_hist_ee.GetMean()+data_hist_ee.GetRMS()
                elif tree.name== 'da_mm':
                    data_hist_mm = tree.getTH1F(lumi, reg.name+i, reg.rvars,50, -150., 150.,  i, '', reg.name)
                    data_Hist_mm = RooDataHist("da_mm","da_mm"+i, RooArgList(x), data_hist_mm)                  
                    m_da_mm = data_hist_mm.GetMean()                               
                    um_da_mm = data_hist_mm.GetMean()-data_hist_mm.GetRMS()
                    uM_da_mm = data_hist_mm.GetMean()+data_hist_mm.GetRMS()   
                                                                                                                                                                              
            data_f_ee, data_e_ee = constructModel(data_Hist_ee, tt_Hist_ee, m_da_ee, um_da_ee, uM_da_ee, True, 'ee',  str(cutsList.index(i)), var, 'data')
            data_f_mm, data_e_mm = constructModel(data_Hist_mm, tt_Hist_mm, m_da_mm, um_da_mm, uM_da_mm, True, 'mm',  str(cutsList.index(i)), var, 'data')
            data_fwhms_ee.append(data_f_ee)
            data_error_ee.append(data_e_ee)
            data_fwhms_mm.append(data_f_mm)
            data_error_mm.append(data_e_mm)
            mc_f_ee, mc_e_ee = constructModel(dy_Hist_ee, dy_BkgHist_ee, m_dy_ee, um_dy_ee, uM_dy_ee, False, 'ee', str(cutsList.index(i)), var, 'MC')
            mc_f_mm, mc_e_mm = constructModel(dy_Hist_mm, dy_BkgHist_mm, m_dy_mm, um_dy_mm, uM_dy_mm, False, 'mm', str(cutsList.index(i)), var, 'MC')
            mc_fwhms_ee.append(mc_f_ee)
            mc_fwhms_mm.append(mc_f_mm)
            mc_error_ee.append(mc_e_ee)
            mc_error_mm.append(mc_e_mm)
            ratio_ee.append(data_f_ee/mc_f_ee)
            ratio_mm.append(data_f_mm/mc_f_mm)
        print 'mc_fwhms_ee', mc_fwhms_ee
        print 'mc_fwhms_ee', mc_error_ee
        print binposition
        print binerror
        plot_da = Canvas.Canvas("met/resolution_plots/"+regname+"_"+var+"_vs_"+str(reg.dependence)+"_fit", "png",0.6, 0.65, 0.85, 0.82)        
        graph_ee = TGraphErrors(len(binposition), array("f", binposition), array("f", data_fwhms_ee), array("f", binerror), array("f", data_error_ee))
        graph_mm = TGraphErrors(len(binposition), array("f", binposition), array("f", data_fwhms_mm), array("f", binerror), array("f", data_error_mm))
        graph_ee.SetMarkerColor(4)
        graph_mm.SetMarkerColor(8)
        graph_ee.SetMarkerStyle(22)
        graph_mm.SetMarkerStyle(23)
        graph_ee.SetTitle("Z->ee")
        graph_mm.SetTitle("Z->mm")
        mg = TMultiGraph() 
        mg.Add(graph_ee)
        mg.Add(graph_mm)
        mg.Draw("AP")
        #mg.GetXaxis().SetTitle('Z q_{T} [GeV]')
        mg.GetYaxis().SetTitle("#sigma (" + label + " ) [GeV]")
        if var == 'uPerp':
            mg.SetMaximum(23.)
        if var == 'uPara':
            mg.SetMaximum(40.)
        leg = TLegend(0.65,0.7,0.85,0.9)
        #leg.SetFillColor(r.kWhite)
        #leg.SetTextFont(42)
        #leg.SetTextSize(0.04)
        #leg.SetLineWidth(0)
        #leg.SetBorderSize(0)
        leg.Draw()
        leg.AddEntry(graph_ee, "Z->ee", "P")
        leg.AddEntry(graph_mm, "Z->#mu#mu", "P")
        plot_da.saveRatio(0, 1, 0, lumi, ratio_ee, binposition)
        #plot_da.save(0, 1, 0, lumi)
        del  graph_ee, graph_mm,  mg, leg
        plot_mc = Canvas.Canvas("met/resolution_plots/"+regname+"_"+var+"_vs_"+str(reg.dependence)+"_MC", "png",0.5, 0.5, 0.8, 0.8)       
        graph_ee = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_fwhms_ee), array("f", binerror), array("f", mc_error_ee))
        graph_mm = TGraphErrors(len(binposition), array("f", binposition), array("f", mc_fwhms_mm), array("f", binerror), array("f", mc_error_mm))
        graph_ee.SetMarkerColor(4)
        graph_mm.SetMarkerColor(8)
        graph_ee.SetMarkerStyle(22)
        graph_mm.SetMarkerStyle(23)
        mg = TMultiGraph() 
        mg.Add(graph_ee)
        mg.Add(graph_mm)
        mg.Draw("AP")
        mg.GetXaxis().SetTitle('Z q_{T} [GeV]')
        mg.GetYaxis().SetTitle("#sigma (" +label + " "+ ")" +" [GeV]")
        if var == 'uPerp':
            mg.SetMaximum(23.)
        if var == 'uPara':
            mg.SetMaximum(40.)
        leg = TLegend(0.17,0.7,0.37,0.9)
        leg.SetFillColor(r.kWhite)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.SetLineWidth(0)
        leg.SetBorderSize(0)
        leg.Draw()
        leg.AddEntry(graph_ee, "Z->ee", "P")
        leg.AddEntry(graph_mm, "Z->#mu#mu", "P")
        plot_mc.addLatex(0.43, 0.84, reg.name)
        plot_mc.save(0, 0, 0, lumi)
        #del leg, binposition, binerror, graph_ee,  mg,error_ee 


