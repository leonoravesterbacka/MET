#####################################################################
######                                                              #
###### 8=========D         OO                        8==========D   #  
###### OO                  ||                                ,88   #
###### ||                  ||                              ,88'    #  
###### ||O---\     ,adPPYb,|| ,adPPYb,d8  ,adPPYba,      ,88'      #
###### ||O---/    a8'    `Y||a8'    `Y88 a8P_____88    ,88'        #
###### ||         8b       ||8b       88 8PP'''''''  ,88'          #
###### \/         '8a,   ,d||'8a,   ,d88 '8b,   ,aa 88'            #
###### 8=========D `'8bbdP'\/  `'YbbdP'Y8  `'Ybbd8'' 8==========D   #
######                       aa,    ,88                             #
######                         'Y8bbdP'                             #
######                                                              #
#####################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TF1, TPaveStats, TGraph, TGraphErrors, TMultiGraph, TLine, TLegend
import math, sys, optparse, copy, re, array
import numpy as np
from array import array


import include.helper     as helper
import include.Region     as Region
import include.Canvas     as Canvas
import include.CutManager as CutManager
import include.Sample     as Sample


if __name__ == "__main__":



    parser = optparse.OptionParser(usage="usage: %prog [opts] FilenameWithSamples", version="%prog 1.0")
    parser.add_option("-t", "--test", action="store_true", dest="test", default=False, help="just do a testrun. takes one variable in one eta for one region")
    parser.add_option('-s', '--samples', action='store', type=str, dest='sampleFile', default='samples.dat', help='the samples file. default \'samples.dat\'')
    (opts, args) = parser.parse_args()

    ## make the options globa.. also the lumi
    (opts, args) = parser.parse_args()

    print 'Going to load DATA and MC trees...'
    dy_mm = ['DYJetsToLL_mm']
    dy_ee = ['DYJetsToLL_ee']

    treeDY_ee = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_ee, 'ee'), 'ee', 0)
    treeDY_mm = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_mm, 'mm'), 'mm', 0)
    trees = [ treeDY_mm, treeDY_ee]
    print 'Trees successfully loaded...'
    lumi = 2.26

   
    gROOT.ProcessLine('.L include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle() 
    cuts = CutManager.CutManager()


    regions = []
    vector = []
    region = Region.region('Type 1 ME_{T}', 
            [cuts.GoodLepton()],
            'met_uPara_zll/zll_pt',
            'met_Type1',
            [range(-100,300,50)],
            False)
    regions.append(region)                 
    region = Region.region('Raw ME_{T}',
            [cuts.GoodLepton()],
            'met_raw_uPara_zll/zll_pt',
            'met_raw',
            [range(-100,300,50)],
            False)
    regions.append(region)                 
   # region = Region.region('Puppi ME_{T}',
   #         [cuts.GoodLepton()],
   #         'metPuppi_uPara_zll/zll_pt',
   #         'metPuppi',
   #         [range(-100,300,50)],
   #         False)
   # regions.append(region)                 
    
    
    
    
    
    
    
    #met = Region.region('ME_{T} [GeV]',
    #        [cuts.GoodLepton()],
    #        'met_uPara_zll/zll_pt',
    #        'met_pt',
    #        [range(-100,300,50)],  
    #        False)
    #regions.append(met)  
    #sumet = Region.region('Sum E_{T} [GeV]',
    #        [cuts.GoodLepton()],
    #        'met_uPara_zll/zll_pt',
    #        'met_sumEt',
    #        [range(-100,300,50)],  
    #        False)
    #regions.append(sumet)
    #nvtx = Region.region('nVtx',
    #        [cuts.GoodLepton()],
    #        'met_uPara_zll/zll_pt',
    #        'nVert',
    #        [range(-100,100,20)],
    #        False)
    #regions.append(nvtx)  
  



    #vtxbins = [[0, 2],[2, 4], [4,6],[6, 8]]
    vtxbins = [[0, 4],[4, 8], [8,12],[12, 16],[16, 20],[20, 24], [24,28],[28, 32],[32, 36],[36, 40]]
    #bins = [ [25, 50], [50, 75], [75, 100], [100, 125]]
    qtbins = [ [25, 50], [50, 75], [75, 100], [100, 125], [125, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 300], [300, 325], [325, 350], [350, 375], [375, 400],[400, 425],  [425, 450], [450, 475], [475, 500]]
    metbins = [[0, 5], [5, 10], [10, 15], [15, 20], [20, 25], [25, 30], [30, 35], [35, 40], [40, 45], [45, 50]]
    sumetbins = [[100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000],[1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 1900], [1900, 2000],[2000, 2100], [2100, 2200]]


    graph_ee = TGraphErrors()
    graph_mm = TGraphErrors()
    leg = TLegend(0.68,0.7,0.88,0.9)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetLineWidth(0)
    leg.SetBorderSize(0)      
    cutsList = []
    mean_ee = []
    error_ee = []
    mean_mm = []
    error_mm = []
    binposition = []
    binerror = []
    for reg in regions:
        graph_ee = TGraphErrors()
        graph_mm = TGraphErrors()
        leg = TLegend(0.68,0.7,0.88,0.9)
        leg.SetFillColor(r.kWhite)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.SetLineWidth(0)
        leg.SetBorderSize(0)
        leg.Draw()
        cutsList = []
        mean_ee = []
        error_ee = []
        mean_mm = []
        error_mm = []
        binposition = []
        binerror = []
        plot_var = Canvas.Canvas("met/scale_plots/"+str(reg.dependence)+"_scale" , "png",0.6, 0.7, 0.8, 0.9)
        #if reg.dependence == "nVert":
        #    for i in vtxbins:
        #        mini = float(min(vtxbins[vtxbins.index(i)]))
        #        maxi = float(max(vtxbins[vtxbins.index(i)]))
        #        mid = float((mini+maxi)/2)
        #        binposition.append(mid)
        #        cutsList.append(cuts.AddList([ reg.dependence+  "<"+ str(maxi) + "&&"+reg.dependence +" >" + str(mini)]))
        #elif reg.dependence == "met_pt": 
        #    for i in metbins:
        #        mini = float(min(metbins[metbins.index(i)]))
        #        maxi = float(max(metbins[metbins.index(i)]))
        #        mid = float((mini+maxi)/2)
        #        binposition.append(mid)
        #        cutsList.append(cuts.AddList([ reg.dependence+  "<"+ str(maxi) + "&&"+reg.dependence +" >" + str(mini)]))
        #elif reg.dependence == "met_sumEt": 
        #    for i in sumetbins:
        #        mini = float(min(sumetbins[sumetbins.index(i)]))
        #        maxi = float(max(sumetbins[sumetbins.index(i)]))
        #        mid = float((mini+maxi)/2)
        #        binposition.append(mid)
        #        cutsList.append(cuts.AddList([ reg.dependence+  "<"+ str(maxi) + "&&"+reg.dependence +" >" + str(mini)]))
        #else: 
        for i in qtbins:
            mini = float(min(qtbins[qtbins.index(i)]))
            maxi = float(max(qtbins[qtbins.index(i)]))
            mid = float((mini+maxi)/2)
            binposition.append(mid)
            cutsList.append(cuts.AddList([   "zll_pt<"+ str(maxi) + "&&zll_pt >" + str(mini)]))

        for i in cutsList:
            print i
            for tree in trees:    
                if tree == treeDY_mm:
                    mean_mm.append(-tree.getTH1F(lumi, reg.name+"_mm_"+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                    error_mm.append(tree.getTH1F(lumi, reg.name+"_mm_"+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
                if tree == treeDY_ee:
                    mean_ee.append(-tree.getTH1F(lumi, reg.name+"_ee_"+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                    error_ee.append(tree.getTH1F(lumi, reg.name+"_ee_"+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
            binerror.append(0.0)
        print 'mean_mm', mean_mm
        print 'mean_ee', mean_ee
        print 'error_ee', error_ee
        print 'error_mm', error_mm
        print binposition
        print 'binerror', binerror

        graph_ee = TGraphErrors(len(binposition), array("f", binposition), array("f", mean_ee), array("f", binerror), array("f", error_ee))
        graph_mm = TGraphErrors(len(binposition), array('f', binposition), array("f", mean_mm), array("f", binerror), array("f", error_mm))
        graph_ee.SetMarkerColor(8)
        graph_ee.SetMarkerStyle(22)
        graph_mm.SetMarkerColor(46)
        graph_mm.SetMarkerStyle(23)
        graph_mm.SetTitle("Z->#mu#mu")
        graph_ee.SetTitle("Z->ee")
        mg = TMultiGraph() 
        mg.Add(graph_mm)
        mg.Add(graph_ee)
        mg.Draw("AP")
        mg.GetXaxis().SetTitle('Z q_{T} [GeV]')
        mg.GetYaxis().SetTitle(reg.name +"-<u_{||}/q_{T}>")
        mg.SetMinimum(0.6)
        mg.SetMaximum(1.3)
        leg.Draw()
        leg.AddEntry(graph_mm, "Z#rightarrow #mu#mu", "P")
        leg.AddEntry(graph_ee, "Z#rightarrow ee", "P")
        plot_var.save(0, 0, 0, lumi)
        del leg, binposition, binerror, graph_ee, graph_mm, mg, plot_var, mean_ee, mean_mm, error_ee, error_mm
