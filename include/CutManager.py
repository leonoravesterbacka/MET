

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.trigMMc = "(HLT_mu17mu8_dz > 0 || HLT_mu27tkmu8 > 0 || HLT_mu17tkmu8_dz > 0)"
      self.twoLeptons = "nLeptons > 1"
      self.bla = "met_pt >= 0"
      self.trigEEc = "(HLT_el17el12_dz > 0 || HLT_ele33ele33 > 0)"
      self.trigEMc = "(HLT_mu8el17 > 0 || HLT_mu17el12 > 0 || HLT_mu30ele30 > 0)"
      self.leptonPt = "t.Lep1_pt_Edge > 20. && t.Lep2_pt_Edge > 20."
      self.leptonDR = "t.lepsDR_Edge > 0.3"       
      self.ECALCrack = "abs(abs(Lep1_eta_Edge) - 1.5) > 0.1 && abs(abs(Lep2_eta_Edge) - 1.5) > 0.1"
      self.leptonsMll = "t.lepsMll_Edge > 20"
      self.goodLepton = self.twoLeptons + "&&" + self.leptonPt + "&&" + self.leptonDR + "&&" + self.ECALCrack + "&&" + self.leptonsMll
      self.goodLeptonMET = self.twoLeptons + "&&" + self.leptonPt + "&&" + "t.lepsMll_Edge > 81 && t.lepsMll_Edge < 101"  
      self.SSee = "(t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 121) "
      self.SSmm = "(t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 169)"
      self.SSOF = "(t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 143)"
      self.SSSF =  "(t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 121 ||  t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 169)"
      self.SS =  "(t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 121 ||  t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 169 ||t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == 143) "
      self.OS =  "(t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == -121 ||  t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == -169 ||t.Lep1_pdgId_Edge * t.Lep2_pdgId_Edge == -143) "
      self.fakes = "((t.Lep1_mcMatchId_Edge * t.Lep2_mcMatchId_Edge)== 0 ) "
      self.fake1 =  "(t.Lep1_mcMatchId_Edge==0) && (t.Lep2_mcMatchId_Edge!=0 )"
      self.fake2 =  "(t.Lep2_mcMatchId_Edge==0) && (t.Lep1_mcMatchId_Edge!=0) " 
      self.twofakes = "(t.Lep1_mcMatchId_Edge == 0) && (t.Lep2_mcMatchId_Edge ==0)"
      self.nophoton = "ngamma == 0 "       
      self.genweights = "genWeight > 0"
      self.notau = "((t.Lep1_mcMatchTau_Edge*t.Lep2_mcMatchTau_Edge) ==0 )"
      self.ee = "t.Lep1_pdgId_Edge*t.Lep2_pdgId_Edge == -121"
      self.mm = "t.Lep1_pdgId_Edge*t.Lep2_pdgId_Edge == -169"
      self.OF = "t.Lep1_pdgId_Edge*t.Lep2_pdgId_Edge == -143"
      self.SF = "(" + self.ee + " || " +  self.mm + ")"
      self.AF = "(" + self.SF + " || " +  self.OF + ")"
      self.nj2 = "(t.nJetSel_Edge >= 2)"
      self.InclusiveCR = "(t.nJetSel_Edge >= 1 && t.lepsMll_Edge > 60 && t.lepsMll_Edge < 120)"
      self.nj0 = "(t.nJetSel_Edge >= 0)"
      self.METJetsSignalRegion = "((met_pt > 150 && t.nJetSel_Edge > 1) || (met_pt > 100 && t.nJetSel_Edge > 2))"
      self.JetsSignalRegion = "(t.nJetSel_Edge >= 0)"
      self.METJetsControlRegion = "(met_pt > 100 && met_pt < 150 && t.nJetSel_Edge == 2)"
      self.RSFOFControlAlternative = "(met_pt > 50 && met_pt < 150 && t.nJetSel_Edge == 2 && t.nBJetLoose35_Edge >=1 )"
      self.METBJetsControlRegion = "(met_pt > 100 && met_pt < 150 && t.nJetSel_Edge == 2 && t.nBJetLoose35_Edge >=2 )"
      self.DYControlRegion = "(met_pt < 50 && t.nJetSel_Edge >= 2)"
      self.blinded = "!"+self.METJetsSignalRegion
      self.DYmet = "(met_pt > 50)"
      self.TestMET = "(met_pt < 170 && met_pt > 160)"
      self.DYmass = "t.lepsMll_Edge > 60 && t.lepsMll_Edge < 120"
      self.ZmassVeto = "(t.lepsMll_Edge < 81 || t.lepsMll_Edge > 101)"
      self.lowmass = "t.lepsMll_Edge > 20 && t.lepsMll_Edge < 70"
      self.Zmass = "t.lepsMll_Edge > 81 && t.lepsMll_Edge < 101"
      self.highmass = "t.lepsMll_Edge > 120"
      self.twolight = "(t.Lep1_mcMatchAny_Edge ==1 && t.Lep2_mcMatchAny_Edge ==1)"
      self.unmatched = "(t.Lep1_mcMatchAny_Edge ==0 && t.Lep2_mcMatchAny_Edge ==0)"
      self.twoheavy = "(t.Lep1_mcMatchAny_Edge >= 4 && t.Lep2_mcMatchAny_Edge >= 4)"
      self.mix =      "((t.Lep1_mcMatchAny_Edge >= 4 && t.Lep2_mcMatchAny_Edge < 4 ) || (t.Lep2_mcMatchAny_Edge >= 4 && t.Lep1_mcMatchAny_Edge <4 ))  "
      self.charm =  "t.Lep2_mcMatchId_Edge == 0 && t.Lep2_mcMatchAny_Edge == 4"
      self.bottom = "t.Lep2_mcMatchId_Edge == 0 && t.Lep2_mcMatchAny_Edge == 5"
      self.prompt = "t.Lep2_mcMatchId_Edge !=0"
      self.fake = "t.Lep2_mcMatchId_Edge ==0 "
      self.central = "(abs(t.Lep1_eta_Edge)<1.4 && abs(t.Lep2_eta_Edge)<1.4)"
      self.forward = "(abs(t.Lep1_eta_Edge)>1.4 || abs(t.Lep2_eta_Edge)>1.4)"
      self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      self.noOStrigger = "((" + self.trigMMc +" && " + self.SSmm + ") || (" + self.trigEEc  +" && " + self.SSee + ") || (" + self.trigEMc +" && " + self.SSOF + "))"
      self.noSignTrigger = "((" + self.trigMMc +") || (" + self.trigEEc  + ") || (" + self.trigEMc + "))"
      self.HT = "(t.htJet35j_Edge > 200)"
      self.triggerHT = "(HLT_pfht200 > 0 || HLT_pfht250 > 0 || HLT_pfht300 > 0 || HLT_pfht300 > 0 || HLT_pfht400>0 || HLT_pfht475>0 || HLT_pfht600>0 || HLT_pfht800>0 || HLT_at51>0 || HLT_at52 >0 || HLT_at53 > 0 || HLT_at55 > 0)"

   def brackets(self, cut):
      return '('+cut+')'

   def AddList(self, cutlist):
      returncut = ''
      for cut in cutlist:
          returncut += cut
          if not cutlist.index(cut) == len(cutlist)-1:
            returncut += ' && '
      return self.brackets(returncut)
  
   def Add(self, cut1, cut2):

      return self.brackets(cut1 + " && " + cut2 )
  
   def MaxRun(self, run):
      
      return self.brackets("run <= %d"%(run))

   def Central(self):
      
      return self.brackets(self.central)

   def Charm(self):

       return self.brackets(self.charm)

   def Light(self):
       
       return self.brackets(self.light)

   def Bottom(self):       
       
       return self.brackets(self.bottom)

   def trigMM(self):

      return self.brackets(self.trigMMc)
 
   def trigEE(self):

      return self.brackets(self.trigEEc)
 
   def trigEM(self):

      return self.brackets(self.trigEMc)
 
   def Forward(self):
      
      return self.brackets(self.forward)

   def DYMass(self):
 
      return self.brackets(self.DYmass)
 
   def GoodLeptonNoTriggerSF(self):

      return self.brackets(self.goodLepton + " && " + self.SF)

   def GoodLeptonNoTriggerOF(self):

      return self.brackets(self.goodLepton + " && " + self.OF) 

   def GoodLeptonNoTriggeree(self):

      return self.brackets(self.goodLepton + " && " + self.ee)

   def GoodLeptonNoTriggermm(self):

      return self.brackets(self.goodLepton + " && " + self.mm) 

   def GoodLeptonNoSignTrigger(self):

      return self.brackets(self.goodLepton + " && " + self.noSignTrigger)

   def GoodLepton(self):

      return self.brackets(self.goodLeptonMET)

   def GoodLeptonNoOSTrigger(self):

      return self.brackets(self.goodLepton + " && " + self.noOStrigger)

   def GoodLeptonAF(self):

      return self.brackets(self.goodLepton + " && " + self.AF + " && " + self.trigger)

   def GoodLeptonSS(self):

      return self.brackets( self.SS )

   def GoodLeptonOS(self):

      return self.brackets( self.OS )

   def GoodLeptonSF(self):

      return self.brackets(self.goodLepton + " && " + self.SF + " && " + self.trigger)

   def GoodLeptonOF(self):

      return self.brackets(self.goodLepton + " && " + self.OF + " && " + self.tigger)

   def GoodLeptonee(self):

      return self.brackets(self.goodLepton + " && " + self.ee + " && " + self.trigger)

   def GoodLeptonmm(self):

      return self.brackets(self.goodLepton + " && " + self.mm + " && " + self.trigger)

   def GoodLeptonSSee(self):

      return self.brackets(self.goodLepton + " && " + self.SSee + " && " + self.trigger)

   def GoodLeptonSSmm(self):                                                               
  
        return self.brackets(self.goodLepton + " && " + self.SSmm + " && " + self.trigger)
  
   def GoodLeptonSSem(self):

      return self.brackets(self.goodLepton + " && " + self.SSem + " && " + self.trigger)   

   def SignalNoMassLeptonSF(self):

      return self.brackets(self.GoodLeptonSF() + " && " + self.METJetsSignalRegion)
   
   def SignalNoMassLeptonOF(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.METJetsSignalRegion)
   
   def SignalNoMassLeptonee(self):

      return self.brackets(self.GoodLeptonee() + " && " + self.METJetsSignalRegion)

   def SignalNoMassLepton(self):

      return self.brackets(self.GoodLeptonNoSignTrigger() + " && " + self.METJetsSignalRegion)

   def ControlNoMassLeptonSF(self):

      return self.brackets(self.GoodLeptonSF() + " && " + self.METJetsControlRegion  + " && " + self.trigger)
   
   def ControlNoMassLeptonOF(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.METJetsControlRegion  + " && " + self.trigger )
   
   def METBJetsOFControlRegion(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.METBJetsControlRegion  + " && " + self.trigger )

   def ControlNoMassLeptonee(self):

      return self.brackets(self.GoodLeptonee() + " && " + self.METJetsControlRegion)

   def ControlNoMassLeptonmm(self):

      return self.brackets(self.GoodLeptonmm() + " && " + self.METJetsControlRegion)

   def DYControlNoMassLeptonSF(self):

      return self.brackets(self.GoodLeptonSF() + " && " + self.DYControlRegion  + " && " + self.trigger)
   
   def DYControlNoMassLeptonOF(self):

      return self.brackets(self.GoodLeptonOF() + " && " + self.DYControlRegion  + " && " + self.trigger)

   def DYControlNoMassLeptonee(self):

      return self.brackets(self.GoodLeptonee() + " && " + self.DYControlRegion)

   def DYControlNoMassLeptonmm(self):

      return self.brackets(self.GoodLeptonmm() + " && " + self.DYControlRegion)

   def SignalLowMassSF(self):

      return self.brackets(self.SignalNoMassLeptonSF() + " && " + self.lowmass)
   
   def SignalLowMassOF(self):

      return self.brackets(self.SignalNoMassLeptonOF() + " && " + self.lowmass)
   
   def SignalLowMassee(self):

      return self.brackets(self.SignalNoMassLeptonee() + " && " + self.lowmass)
   
   def SignalLowMass(self):

      return self.brackets(self.SignalNoMassLepton() + " && " + self.lowmass)

   def SignalZMassSF(self):

      return self.brackets(self.SignalNoMassLeptonSF() + " && " + self.Zmass)
   
   def SignalZMassOF(self):

      return self.brackets(self.SignalNoMassLeptonOF() + " && " + self.Zmass)
   
   def SignalZMassee(self):

      return self.brackets(self.SignalNoMassLeptonee() + " && " + self.Zmass)
   
   def SignalZMass(self):

      return self.brackets(self.SignalNoMassLepton() + " && " + self.Zmass)

   def SignalHighMassSF(self):

      return self.brackets(self.SignalNoMassLeptonSF() + " && " + self.highmass)
   
   def SignalHighMassOF(self):

      return self.brackets(self.SignalNoMassLeptonOF() + " && " + self.highmass)
   
   def SignalHighMassee(self):

      return self.brackets(self.SignalNoMassLeptonee() + " && " + self.highmass)
   
   def SignalHighMass(self):

      return self.brackets(self.SignalNoMassLepton() + " && " + self.highmass)
 
   def ControlLowMassSF(self):

      return self.brackets(self.ControlNoMassLeptonSF() + " && " + self.lowmass)
   
   def ControlLowMassOF(self):

      return self.brackets(self.ControlNoMassLeptonOF() + " && " + self.lowmass)
   
   def ControlLowMassee(self):

      return self.brackets(self.ControlNoMassLeptonee() + " && " + self.lowmass)
   
   def ControlLowMassmm(self):

      return self.brackets(self.ControlNoMassLeptonmm() + " && " + self.lowmass)

   def ControlZMassSF(self):

      return self.brackets(self.ControlZMassLeptonSF() + " && " + self.Zmass)
   
   def ControlZMassOF(self):

      return self.brackets(self.ControlZMassLeptonOF() + " && " + self.Zmass)
   
   def ControlZMassee(self):

      return self.brackets(self.ControlZMassLeptonee() + " && " + self.Zmass)
   
   def ControlZMassmm(self):

      return self.brackets(self.ControlNoMassLeptonmm() + " && " + self.Zmass)

   def ControlHighMassSF(self):

      return self.brackets(self.ControlNoMassLeptonSF() + " && " + self.highmass)
   
   def ControlHighMassOF(self):

      return self.brackets(self.ControlNoMassLeptonOF() + " && " + self.highmass)
   
   def ControlHighMassee(self):

      return self.brackets(self.ControlNoMassLeptonee() + " && " + self.highmass)
   
   def Control2JetsSF(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonSF()  + " && " + self.trigger)
   
   def Control2JetsOF(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonOF()  + " && " + self.trigger)
   
   def Control2JetsAF(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonAF()  + " && " + self.trigger)
   
   def RSFOFControlRegion(self):

      return self.brackets(self.RSFOFControlAlternative + " && " + self.trigger)
   
   def Control2Jetsee(self):

      return self.brackets(self.nj2 + " && " + self.GoodLeptonee())
   
   def Control2Jetsmm(self):
  
	  return self.brackets(self.nj2 + " && " + self.GoodLeptonmm())
      
   def Control2JetsMETSF(self):

      return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonSF())
   
   def Control2JetsMETOF(self):

      return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonOF())
   
   def Control2JetsMETee(self):

      return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonee())
   
   def Control2JetsMETmm(self):

	  return self.brackets(self.nj2 + " && " + self.DYmet + " && " + self.GoodLeptonmm())
      
   def InclusiveCROF(self):

      return self.brackets(self.InclusiveCR + " && "  + self.GoodLeptonOF())
   
   def InclusiveCRSF(self):

	  return self.brackets(self.InclusiveCR + " && "  + self.GoodLeptonSF())
  
   def InclusiveCRee(self):

      return self.brackets(self.InclusiveCR + " && " + self.GoodLeptonee())
    
   def InclusiveCRmm(self):

      return self.brackets(self.InclusiveCR + " && " + self.GoodLeptonmm())

   def DYFake(self):
                                                                                         
      return self.brackets(self.DYControlRegion + " && " + self.GoodLeptonNoSignTrigger() + " && " + self.fakes )
   
   def DYSameSign(self):

      return self.brackets(self.DYControlRegion + " && " + self.GoodLeptonNoOSTrigger())

   def DYSameSignEE(self): 
  
      return self.brackets(self.DYControlRegion + " && " + self.GoodLeptonNoOSTrigger() + " && " + self.SSee)
  
   def DYSameSignOF(self):
   
      return self.brackets(self.DYControlRegion + " && " + self.GoodLeptonNoOSTrigger() + " && " + self.SSOF)
  
   def DYSameSignMM(self):

      return self.brackets(self.DYControlRegion + " && " + self.GoodLeptonNoOSTrigger() + " && " + self.SSmm)

   def SRSameSignEE(self):

      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSee)
    
   def SRSameSignMM(self): 
  
      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSmm)
  
   def SRSameSignOF(self):
     
      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSOF)
  
   def SRSameSignSF(self):
     
      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSSF)
  
   def SRSameSign(self):
         
     #return self.brackets(self.METJetsSignalRegion + " && " + self.trigger)
     return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoOSTrigger() + " && "  + self.SS)
   def SRFake(self):
                                                                                         
      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLepton() + " && " + self.fakes )

   def SRFakeNoMETCut(self):
                                                                                         
      return self.brackets(self.JetsSignalRegion + " && " + self.GoodLepton()+ " && "+ self.twofakes )

   def SRFakeNoMETCut_SS(self):

      return self.brackets(self.JetsSignalRegion + " && " + self.GoodLeptonNoOSTrigger() + " && " + self.twofakes  )

   def SRFakeMM(self):

      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoSignTrigger() + " && " + self.fakes+ " && " + self.mm  )

   def SRFakeNoSignTrigger(self):

      return self.brackets(self.METJetsSignalRegion + " && " + self.GoodLeptonNoSignTrigger() + " && " + self.fakes  )

   def TTSameSignEE(self):
    
      return self.brackets(self.METJetsControlRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSee)

   def TTSameSignMM(self): 
  
      return self.brackets(self.METJetsControlRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSmm)
  
   def TTSameSignOF(self):

      return self.brackets(self.METJetsControlRegion + " && " + self.GoodLeptonNoOSTrigger()+ " && " + self.SSOF)
  
   def TTSameSign(self): 
  
      return self.brackets(self.METJetsControlRegion + " && " + self.GoodLeptonNoOSTrigger()) 
  
   def TTFake(self):
                                                                                         
      return self.brackets(self.METJetsControlRegion + " && " + self.GoodLeptonNoSignTrigger() + " && " + self.fakes )

   def charm(self):
    
      return self.brackets(self.c)  
