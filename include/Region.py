from ROOT import TGraphErrors
import sys


class collection:
    def __init__(self, bins, vname):
        self.bins  = bins
        self.vname = vname
        self.cen_mc, self.cen_mc_gr = 0, 0
        self.fwd_mc, self.fwd_mc_gr = 0, 0
        self.cen_da, self.cen_da_gr = 0, 0
        self.fwd_da, self.fwd_da_gr = 0, 0

    def getHisto(self, dataMC, eta):
        if   (dataMC, eta) == ('MC'  , 'central'): 
            return self.cen_mc
        elif (dataMC, eta) == ('MC'  , 'forward'): 
            return self.fwd_mc
        elif (dataMC, eta) == ('DATA', 'central'): 
            return self.cen_da
        elif (dataMC, eta) == ('DATA', 'forward'): 
            return self.fwd_da

    def getGraph(self, dataMC, eta):
        if   (dataMC, eta) == ('MC'  , 'central'): 
            return self.cen_mc_gr
        elif (dataMC, eta) == ('MC'  , 'forward'): 
            return self.fwd_mc_gr
        elif (dataMC, eta) == ('DATA', 'central'): 
            return self.cen_da_gr
        elif (dataMC, eta) == ('DATA', 'forward'): 
            return self.fwd_da_gr

    def setHisto(self, histo, dataMC, eta):
        if   (dataMC, eta) == ('MC'  , 'central'): 
            self.cen_mc = histo
            self.cen_mc_gr = TGraphErrors(histo)
        elif (dataMC, eta) == ('MC'  , 'forward'): 
            self.fwd_mc = histo
            self.fwd_mc_gr = TGraphErrors(histo)
        elif (dataMC, eta) == ('DATA', 'central'): 
            self.cen_da = histo
            self.cen_da_gr = TGraphErrors(histo)
        elif (dataMC, eta) == ('DATA', 'forward'): 
            self.fwd_da = histo
            self.fwd_da_gr = TGraphErrors(histo)
        else:
            print 'you are not calling setHisto correctly'

    def printValues(self):
        for _histo in [self.cen_mc, self.cen_da, self.fwd_mc, self.fwd_da]:
            if _histo == 0: continue
            if _histo in [self.cen_mc, self.cen_da]:
                eta = 'central'
            else:
                eta = 'forward'
            for _bin in range(1,_histo.GetNbinsX()+1):
                print '%-6s in [%.0f, %.0f] in %s: %.3f +- %.3f' %( self.vname,
                      _histo.GetXaxis().GetBinLowEdge(_bin),
                      _histo.GetXaxis().GetBinUpEdge (_bin),
                      eta,
                      _histo.GetBinContent(_bin),
                      _histo.GetBinError(_bin))

    def saveInFile(self, pattern, systErr, findBin = 0):
        print 'writing calculated values into file...'
        filename = 'ingredients.dat'
        f = open(filename, 'r')
        lines = f.readlines()
        newlines = []
        for line in lines:
            appended = False
            for t in ['MC', 'DATA'] if self.cen_da else ['MC']:
                if all(s in line for s in pattern+[t]):
                    newlines.append('%-6s \t %-15s %-6s \t %.4f \t %.4f \t %-6s \t %.4f \t %.4f \t %-6s\n' %(
                            line.split()[0], line.split()[1], t,
                            self.getHisto(t, 'central').GetBinContent(1 if not findBin else self.getHisto(t, 'central').FindBin(findBin)),
                            self.getHisto(t, 'central').GetBinError  (1 if not findBin else self.getHisto(t, 'central').FindBin(findBin)),
                            str(systErr),
                            self.getHisto(t, 'forward').GetBinContent(1 if not findBin else self.getHisto(t, 'forward').FindBin(findBin)),
                            self.getHisto(t, 'forward').GetBinError  (1 if not findBin else self.getHisto(t, 'forward').FindBin(findBin)),
                            str(systErr) ))
                    appended = True
    
            if not appended:
                newlines.append(line)
        f.close()
        g = open(filename, 'w')
        g.writelines(newlines)

class region():
    def __init__(self, name, cuts, rvars, dependence, varname, bins,  doData):
        self.name      = name
        self.cuts      = cuts
        self.rvars     = rvars
        self.dependence = dependence
        self.varname   = varname
        self.bins      = bins
        self.doData    = doData
        #if len(rvars) is not len(bins):
        #    print 'length of variables and bins has to be equal'
        #    sys.exit('exiting!')
        self.setVariables()

    def setVariables(self):
        for v in self.rvars:

            if v == 'met_pt':
                self.met_pt   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pt_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pt_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_pt_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_pt_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.met_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Up_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Up_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_Up_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_Up_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.met_Down   = collection(self.bins[self.rvars.index(v)], v)          
                self.met_Down_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Down_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_Down_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_Down_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.met_pt   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pt_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pt_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_pt_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_pt_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_pt_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_pt_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_pt_da = collection(self.bins[self.rvars.index(v)], v)         
                self.met_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Up_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Up_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Up_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Up_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Up_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Up_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Up_da = collection(self.bins[self.rvars.index(v)], v)         
                self.met_Down   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Down_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Down_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_Down_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Down_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Down_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Down_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_Down_da = collection(self.bins[self.rvars.index(v)], v)         
            if v == 'zll_pt':
                self.zll_pt   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_pt_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_pt_da = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_pt_Up_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_pt_Up_da = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_Up_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_Up_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_Down   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_pt_Down_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_pt_Down_da = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_Down_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_pt_Down_ewk = collection(self.bins[self.rvars.index(v)], v)    
            if v == 'zll_mass':
                self.zll_mass   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_mass_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_mass_da = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_mass_Up_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_mass_Up_da = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_Up_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_Up_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_Down   = collection(self.bins[self.rvars.index(v)], v)       
                self.zll_mass_Down_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.zll_mass_Down_da = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_Down_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.zll_mass_Down_ewk = collection(self.bins[self.rvars.index(v)], v)    
            if v == 'nVert':
                self.nVert   = collection(self.bins[self.rvars.index(v)], v)
                self.nVert_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.nVert_da = collection(self.bins[self.rvars.index(v)], v)    
                self.nVert_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.nVert_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.nVert_wg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVert_zg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVert_wj = collection(self.bins[self.rvars.index(v)], v) 
                self.nVert_ttg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVert_qcd = collection(self.bins[self.rvars.index(v)], v) 
                self.nVert_gjets = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertUp   = collection(self.bins[self.rvars.index(v)], v)
                self.nVertUp_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.nVertUp_da = collection(self.bins[self.rvars.index(v)], v)    
                self.nVertUp_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.nVertUp_ewk = collection(self.bins[self.rvars.index(v)], v)   
                self.nVertUp_wg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertUp_zg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertUp_wj = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertUp_ttg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertUp_qcd = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertUp_gjets = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertDown   = collection(self.bins[self.rvars.index(v)], v)
                self.nVertDown_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.nVertDown_da = collection(self.bins[self.rvars.index(v)], v)    
                self.nVertDown_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.nVertDown_ewk = collection(self.bins[self.rvars.index(v)], v)   
                self.nVertDown_wg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertDown_zg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertDown_wj = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertDown_ttg = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertDown_qcd = collection(self.bins[self.rvars.index(v)], v) 
                self.nVertDown_gjets = collection(self.bins[self.rvars.index(v)], v) 
            if v == 'met_uPerp_zll':
                self.met_uPerp_zll   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_zll_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_zll_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_zll_Up_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_zll_Up_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_Up_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_Up_ewk = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_Down   = collection(self.bins[self.rvars.index(v)], v)                
                self.met_uPerp_zll_Down_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_zll_Down_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_Down_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPerp_zll_Down_ewk = collection(self.bins[self.rvars.index(v)], v)    
            if v == 'met_uPara_zll':
                self.met_uPara_zll     = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_zll_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_zll_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPara_zll_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPara_zll_ewk = collection(self.bins[self.rvars.index(v)], v)         
                self.met_uPara_zll_Up     = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_zll_Up_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_zll_Up_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPara_zll_Up_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPara_zll_Up_ewk = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPara_zll_Down     = collection(self.bins[self.rvars.index(v)], v)        
                self.met_uPara_zll_Down_dy   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_zll_Down_da = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPara_zll_Down_tt = collection(self.bins[self.rvars.index(v)], v)    
                self.met_uPara_zll_Down_ewk = collection(self.bins[self.rvars.index(v)], v)       
            if v == 'met_uPara_gamma':
                self.met_uPara_gamma   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_da = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPara_gamma_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_Up_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_Up_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_Up_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Up_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Up_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Up_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Up_da = collection(self.bins[self.rvars.index(v)], v)      
                self.met_uPara_gamma_Down   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_Down_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_Down_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPara_gamma_Down_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Down_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Down_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Down_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPara_gamma_Down_da = collection(self.bins[self.rvars.index(v)], v)      
            if v == 'met_uPerp_gamma':
                self.met_uPerp_gamma   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_gamma_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_gamma_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_da = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPerp_gamma_qcd = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPerp_gamma_Up   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_gamma_Up_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_gamma_Up_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Up_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Up_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Up_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Up_da = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPerp_gamma_Up_qcd = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPerp_gamma_Down   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_gamma_Down_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.met_uPerp_gamma_Down_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Down_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Down_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Down_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.met_uPerp_gamma_Down_da = collection(self.bins[self.rvars.index(v)], v)       
                self.met_uPerp_gamma_Down_qcd = collection(self.bins[self.rvars.index(v)], v)       
            if v == 'gamma_pt':
                self.gamma_pt   = collection(self.bins[self.rvars.index(v)], v)
                self.gamma_pt_qcd   = collection(self.bins[self.rvars.index(v)], v)
                self.gamma_pt_gjets   = collection(self.bins[self.rvars.index(v)], v)
                self.gamma_pt_wg = collection(self.bins[self.rvars.index(v)], v)  
                self.gamma_pt_zg = collection(self.bins[self.rvars.index(v)], v)  
                self.gamma_pt_wj = collection(self.bins[self.rvars.index(v)], v)  
                self.gamma_pt_ttg = collection(self.bins[self.rvars.index(v)], v)  
                self.gamma_pt_da = collection(self.bins[self.rvars.index(v)], v)      

