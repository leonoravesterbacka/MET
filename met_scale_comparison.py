#####################################################################
############################ MET ####################################
#####################################################################
# this one compares different met types in the same plot, for ee and mumu separatly 
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
    dy_mm_old = ['DYJetsToLL_mm_old']
    dy_ee_old = ['DYJetsToLL_ee_old']

    ee = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_ee, 'ee'), 'ee', 0)
    mm = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_mm, 'mm'), 'mm', 0)
    ee_old = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_ee_old, 'ee'), 'ee', 1)
    mm_old = Sample.Tree(helper.selectSamples(opts.sampleFile, dy_mm_old, 'mm'), 'mm', 1)
    #trees = [ mm_old, ee_old]
    trees = [ mm]
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
            'met',
            [range(-100,300,50)],
            False)
    regions.append(region)                 
    region = Region.region('Raw ME_{T}',
            [cuts.GoodLepton()],
            'met_raw_uPara_zll/zll_pt',
            'metRaw',
            [range(-100,300,50)],
            False)
    regions.append(region)                 
    #region = Region.region('Raw ME_{T} 74X',
    #        [cuts.GoodLepton()],
    #        'met_raw_uPara_zll/zll_pt',
    #        'met_raw_74X',
    #        [range(-100,300,50)],
    #        False)
    #regions.append(region)                 
#    region = Region.region('Puppi ME_{T} Raw',
#            [cuts.GoodLepton()],
#            '(((-metPuppi_rawPt*sin(metPuppi_rawPhi))-(zll_pt*sin(zll_phi)))*zll_pt*(sin(zll_phi)) + ((-metPuppi_rawPt*cos(metPuppi_rawPhi)-(zll_pt*cos(zll_phi)))*zll_pt*cos(zll_phi)))/(zll_pt*zll_pt)',
#            'metPuppiRaw',
#            [range(-100,300,50)],                                                                                                                                                                   
#            False)
#    regions.append(region)   
    
    
    
    #qtbins = [ [25, 50], [50, 75], [75, 100]]
    qtbins = [ [25, 50], [50, 75], [75, 100], [100, 125], [125, 150], [150, 175], [175, 200], [200, 225], [225, 250], [250, 275],[275, 300], [300, 325], [325, 350], [350, 375], [375, 400],[400, 425],  [425, 450], [450, 475], [475, 500]]


    graph_type1 = TGraphErrors()
    graph_raw = TGraphErrors()
    graph_puppi = TGraphErrors()
    graph_puppi_raw = TGraphErrors()
    leg = TLegend(0.68,0.7,0.88,0.9)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetLineWidth(0)
    leg.SetBorderSize(0)      
    cutsList = []
    mean_type = []
    error_type = []
    mean_raw = []
    error_raw = []
    mean_puppi = []
    error_puppi = []
    mean_puppi_raw = []
    error_puppi_raw = []
    mean_raw_74X = []
    error_raw_74X = []
    binposition = []
    mean_type_74X = []
    error_type_74X = []
    binerror = []
    for i in qtbins:
        mini = float(min(qtbins[qtbins.index(i)]))
        maxi = float(max(qtbins[qtbins.index(i)]))
        mid = float((mini+maxi)/2)
        binposition.append(mid)
        binerror.append(0.0)
        cutsList.append(cuts.AddList([   "zll_pt<"+ str(maxi) + "&&zll_pt >" + str(mini)]))

    for tree in trees:
        treename = tree.name.lower()
        print 'doing tree ', treename, tree.isData
        leg = TLegend(0.55,0.7,0.75,0.9)
        leg.SetFillColor(r.kWhite)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.SetLineWidth(0)
        leg.SetBorderSize(0)      
        
        mean_type = []
        mean_raw = []
        mean_puppi = []
        mean_puppi_raw = []
        mean_raw_74X = []
        mean_raw_74X = []
        error_type = []
        error_raw = []
        error_puppi = []
        error_puppi_raw = []
        error_raw_74X = []
        error_raw_74X = []
        for reg in regions:
            print 'doing reg, ', reg.name
            for i in cutsList:
                print i   
                if reg.dependence == 'met':
                    mean_type.append(-tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                    error_type.append(tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
                if reg.dependence == 'metPuppi':
                    mean_puppi.append(-tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                    error_puppi.append(tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
                if reg.dependence == 'metRaw':
                    mean_raw.append(-tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                    error_raw.append(tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
                if reg.dependence == 'metPuppiRaw' :
                    mean_puppi_raw.append(-tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                    error_puppi_raw.append(tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
                #if reg.dependence == 'metRaw':
                #    print 'hereeeeeee'
                #    mean_raw_74X.append(-tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMean())
                #    error_raw_74X.append(tree.getTH1F(lumi, reg.name+str(cutsList.index(i)), reg.rvars, 100, -4, 4,  i, '', reg.name).GetMeanError())
       
       
        print 'mean_type', treename,  mean_type
        print 'error_type', treename, error_type
        #print 'mean_raw', treename, mean_raw
        #print 'error_raw', treename, error_raw
       # print 'mean_puppiRaw', mean_puppi_raw
       # print 'binposition', binposition
       # print 'binerror', binerror
    mean_raw_74X_mm = [0.92519138875378593, 0.94887751607291504, 0.95671773579949815, 0.96110654839151133, 0.96282807528355452, 0.96435372180434109, 0.96026212355661478, 0.9604774540216231, 0.96065045476710398, 0.95784148397976387, 0.96090156960901552, 0.95206042343364761, 0.95180844615073745, 0.94904832949048323, 0.95857749040210127, 0.95134981124011664, 0.94753400573702218, 0.95250031332247143, 0.94869486948694859]
    error_raw_74X_mm = [0.00061478655264394935, 0.00057036260452820856, 0.00064338751440089497, 0.00076242702769798419, 0.00092067786316567854, 0.00111876085897293, 0.0013383936550002485, 0.0015859686672233008, 0.0018453122938648036, 0.0024106783606974444, 0.0027100334246656678, 0.0031244847373714762, 0.0034746042904681257, 0.0043017647080849939, 0.0045717545124740189, 0.0050404989082165142, 0.0063215525155510925, 0.0067983089333386202, 0.0099133802804301455]
    mean_raw_74X_ee = [0.92673039638346355, 0.94915183552666271, 0.95836585659489415, 0.96290730134751601, 0.96566412650256916, 0.9652168910276363, 0.9678388152605516, 0.96558895862996441, 0.96430348917244668, 0.96157525864946036, 0.96110029412983133, 0.96197040198719519, 0.9634973547606015, 0.96359858208043037, 0.96428309497616427, 0.95270550310845026, 0.95788778877887792, 0.96204620462046209, 0.96883793642522142]
    error_raw_74X_ee = [0.00073299988574178701, 0.00067923965820228054, 0.00076393582821369151, 0.00090059133185762391, 0.0010852710276287481, 0.0013037342372973061, 0.0015397545284986667, 0.0017895295631191019, 0.0021726386412524996, 0.0025415682177892065, 0.0031276605106428756, 0.0034922708641963317, 0.0041405147945084457, 0.004685304001334681, 0.0053839724596755241, 0.0064847565977753244, 0.0070982030566530754, 0.0064132144204404745, 0.0078615601897505816]
    mean_type_74X_mm = [1.0115215742324328, 1.0346102352942947, 1.039619567457134, 1.042818916228005, 1.04378881458345, 1.0452758843723566, 1.0411122265351223, 1.0411342191443509, 1.0431233701528611, 1.0406224433573206, 1.043561890435619, 1.0342425546902516, 1.0357648009698928, 1.036574890365749, 1.0420286926651849, 1.0322672555025287, 1.0308133617100028, 1.040731921293395, 1.0297029702970295]
    error_type_74X_mm = [0.00064320180657504908, 0.00059587588537972998, 0.00067312868284676714, 0.00079902596203431004, 0.00096273077023077322, 0.0011775970604070611, 0.0014093830953413806, 0.0016838528350346031, 0.0019594543651900347, 0.0025518632016700874, 0.0028791366109178975, 0.0034212248066974016, 0.0037222692008727983, 0.004613026747157919, 0.0050616493948883188, 0.0052591479237458025, 0.0066103058079465277, 0.0077491303253399389, 0.010999281596053729]
    mean_type_74X_ee =  [1.0099155969199087, 1.0281505416188899, 1.0343471040962799, 1.0373271696025286, 1.0397392151457927, 1.0396520349320613, 1.0437654830718412, 1.0432874885095424, 1.043906743615538, 1.0424296362220491, 1.0410373254898715, 1.0444436669815391, 1.043899364810853, 1.0426109277594424, 1.051411807847451, 1.0389131936449458, 1.0450165016501651, 1.0401540154015398, 1.0505471599791558]
    error_type_74X_ee = [0.00077431033290634328, 0.00071651871525676481, 0.0008055237963884407, 0.00095413389338820859, 0.0011513141033270336, 0.0013836510386203574, 0.0016434751412119271, 0.0019209594323548842, 0.0023142857387496644, 0.0027495983396691897, 0.0033592557962893667, 0.0038541293318729145, 0.0044556867293860805, 0.0050367500617548534, 0.0057014981516678035, 0.0072576842669976171, 0.0069131235210610998, 0.0066119507293136461, 0.0082715902135763134]





    plot_var_type = Canvas.Canvas("met/scale_plots/met_type1_76vs74"+str(treename) , "png",0.6, 0.7, 0.8, 0.9)                                              
    graph_type = TGraphErrors(len(binposition), array("f", binposition), array("f", mean_type), array("f", binerror), array("f", error_type))
    if treename == 'ee':
        graph_type_74X = TGraphErrors(len(binposition), array('f', binposition), array("f", mean_type_74X_ee), array("f", binerror), array("f", error_type_74X_ee))
    if treename == 'mm':
        graph_type_74X = TGraphErrors(len(binposition), array('f', binposition), array("f", mean_type_74X_mm), array("f", binerror), array("f", error_type_74X_mm))
    
    graph_type.SetMarkerColor(46)
    graph_type_74X.SetMarkerColor(9)
    graph_type.SetMarkerStyle(23)
    graph_type_74X.SetMarkerStyle(22)
    graph_type.SetTitle("ME_{T} Type 1 76X")
    graph_type_74X.SetTitle("ME_{T} Type 1 74X")
    mg = TMultiGraph() 
    mg.Add(graph_type_74X)
    mg.Add(graph_type)
    mg.Draw("AP")
    print treename
    if treename == 'mm':
        mg.GetXaxis().SetTitle('Z#rightarrow #mu#mu q_{T} [GeV]')
    else:
        mg.GetXaxis().SetTitle('Z#rightarrow ee q_{T} [GeV]')
    
    mg.GetYaxis().SetTitle("-<u_{||}/q_{T}>")
    mg.SetMinimum(0.6)
    mg.SetMaximum(1.3)
    leg.Draw()
    leg.AddEntry(graph_type, "ME_{T} Type 1 76X", "P")
    leg.AddEntry(graph_type_74X, "ME_{T} Type 1 74X" , "P")
    
    plot_var_type.save(0, 0, 0, lumi)
    del leg, mg, plot_var_type

    plot_var_raw = Canvas.Canvas("met/scale_plots/met_raw_76vs74"+str(treename) , "png",0.6, 0.7, 0.8, 0.9)                                              
    graph_raw = TGraphErrors(len(binposition), array('f', binposition), array("f", mean_raw), array("f", binerror), array("f", error_raw))
    if treename == 'ee':
        graph_raw_74X = TGraphErrors(len(binposition), array('f', binposition), array("f", mean_raw_74X_ee), array("f", binerror), array("f", error_raw_74X_ee))
    if treename == 'mm':
        graph_raw_74X = TGraphErrors(len(binposition), array('f', binposition), array("f", mean_raw_74X_mm), array("f", binerror), array("f", error_raw_74X_mm))
    
    graph_raw.SetMarkerColor(46)
    graph_raw_74X.SetMarkerColor(9)
    graph_raw.SetMarkerStyle(23)
    graph_raw_74X.SetMarkerStyle(22)
    graph_raw.SetTitle("ME_{T} Raw 76X")
    graph_raw_74X.SetTitle("ME_{T} Raw 74X")
    mg = TMultiGraph() 
    mg.Add(graph_raw_74X)
    mg.Add(graph_raw)
    mg.Draw("AP")
    print treename
    if treename == 'mm':
        mg.GetXaxis().SetTitle('Z#rightarrow #mu#mu q_{T} [GeV]')
    else:
        mg.GetXaxis().SetTitle('Z#rightarrow ee q_{T} [GeV]')
    
    mg.GetYaxis().SetTitle("-<u_{||}/q_{T}>")
    mg.SetMinimum(0.6)
    mg.SetMaximum(1.3)
    leg = TLegend(0.6,0.7,0.8,0.9)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetLineWidth(0)
    leg.SetBorderSize(0)      
    
    leg.Draw()
    leg.AddEntry(graph_raw, "ME_{T} Raw 76X", "P")
    leg.AddEntry(graph_raw_74X, "ME_{T} Raw 74X", "P")
    
    plot_var_raw.save(0, 0, 0, lumi)
