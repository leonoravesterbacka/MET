

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.filters = "(Flag_badMuonFilter > 0 && Flag_badChargedHadronFilter > 0 && Flag_HBHENoiseIsoFilter > 0 && Flag_HBHENoiseFilter > 0 &&   Flag_goodVertices > 0 && Flag_eeBadScFilter > 0 && Flag_EcalDeadCellTriggerPrimitiveFilter > 0 && Flag_CSCTightHalo2015Filter > 0)"
      self.trigMMc = "(HLT_mu17mu8_dz > 0 || HLT_mu27tkmu8 > 0 || HLT_mu17tkmu8_dz > 0)"
      self.twoLeptons = "nLeptons > 1"
      #self.twoGoodLeptons = "(nLepGood20 ==2 )"
      self.twoGoodLeptons = "((nLepGood20 ==2 ) && (abs(lep_eta < 2.4))  && abs(abs(lep_eta) - 1.5) > 0.1) "
      self.noLeptons = "nLepGood10 ==0"
      self.bla = "met_pt >= 0"
      self.trigEEc = "(HLT_el17el12_dz > 0 || HLT_ele33ele33 > 0)"
      self.trigEMc = "(HLT_mu8el17 > 0 || HLT_mu17el12 > 0 || HLT_mu30ele30 > 0)"
      self.leptonPt = "genLep_pt > 20. "
      self.leptonDR = "t.lepsDR_Edge > 0.3"       
      self.leptonEta = "(abs(genLep_pdgId) == 11 ) && (abs(abs(genLep_eta) - 1.5) > 0.1 && abs(abs(genLep_eta) - 1.5) > 0.1)  || ((abs(genLep_pdgId) == 13 ) && (abs(genLep_eta < 2.1)))"
      self.leptonsMll = "t.lepsMll_Edge > 20"
      self.goodLepton =   self.leptonPt + "&&" + self.leptonDR + "&&" + self.leptonEta + "&&" + self.leptonsMll
      self.nj2 = "(t.nJetSel_Edge >= 2)"
      self.jet1_pt = "(jet1_pt >= 0) && (jet2_pt < 0)"
      self.jet2_pt = "(jet2_pt >= 0) && (jet1_pt>= 0)"
      self.InclusiveCR = "(t.nJetSel_Edge >= 1 && t.lepsMll_Edge > 60 && t.lepsMll_Edge < 120)"
      self.nj0 = "(t.nJetSel_Edge >= 0)"
      self.puWeight = "puWeight"
      self.gamma = "(gamma_pt > 50)" + " && "  +  "(ngamma == 1)"+ "&&" +  "(gamma_r9 < 1.0)" + " && " + "(gamma_r9 > 0.9) " + " && "+ "gamma_hOverE < 0.05" + " && " + "gamma_sigmaIetaIeta < 0.01 " 
      self.noLeptons = "nLepGood10 ==0"
      self.Zmasswindow = "zll_mass > 81 && zll_mass < 101"
      self.Zpt = "zll_pt > 50 "
      self.zcut = "zll_pt < 0.5 "
      self.centralPhoton = "(abs(gamma_eta)<1.4) "
      self.central = "(abs(t.Lep1_eta_Edge)<1.4 && abs(t.Lep2_eta_Edge)<1.4)"
      self.forward = "(abs(t.Lep1_eta_Edge)>1.4 || abs(t.Lep2_eta_Edge)>1.4)"
      #self.trigger = "((" + self.trigMMc + " && " + self.mm + ") || (" + self.trigEEc + " && " + self.ee + ") || (" + self.trigEMc + " && " + self.OF + "))"
      #self.noOStrigger = "((" + self.trigMMc +" && " + self.SSmm + ") || (" + self.trigEEc  +" && " + self.SSee + ") || (" + self.trigEMc +" && " + self.SSOF + "))"
      #self.noSignTrigger = "((" + self.trigMMc +") || (" + self.trigEEc  + ") || (" + self.trigEMc + "))"
      #self.triggerHT = "(HLT_pfht200 > 0 || HLT_pfht250 > 0 || HLT_pfht300 > 0 || HLT_pfht300 > 0 || HLT_pfht400>0 || HLT_pfht475>0 || HLT_pfht600>0 || HLT_pfht800>0 || HLT_at51>0 || HLT_at52 >0 || HLT_at53 > 0 || HLT_at55 > 0)"
      self.trigger12090 = "HLT_Photon120>0 && HLT_Photon90>0"
      self.trigger120 = "(HLT_Photon120>0  && gamma_pt > 85 && gamma_pt < 120  )* (HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale)"                     
      self.trigger9075 = "HLT_Photon90>0 && HLT_Photon75>0"
      self.trigger90 = "(HLT_Photon90>0  && gamma_pt > 75 && gamma_pt < 90 )* (HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale)"                     
      self.trigger7550 = "HLT_Photon75>0 && HLT_Photon50>0"
      self.trigger75 = "(HLT_Photon75>0  && gamma_pt > 65 && gamma_pt < 80 )* (HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale)"                     
      self.trigger5030 = "HLT_Photon50>0 && HLT_Photon30>0"
      self.trigger50 = "(HLT_Photon50>0 && gamma_pt > 35 && gamma_pt < 70)* (HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale) "                     
      self.trigger30 = "(HLT_Photon30>0 && gamma_pt > 30 )* (HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v_Prescale)"                     
      self.triggerPhoton = "(HLT_Photon120 == 1 || HLT_Photon90==1 ||HLT_Photon75 == 1 || HLT_Photon50 ==1 || HLT_Photon30 ==1 )"
      self.triggerLepton = ""
      #self.triggerLepton = "( HLT_DoubleMu > 0 || HLT_DoubleEG > 0)"

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
  
   def gammas(self):
                                                                                                                                                     
      return self.brackets(self.centralPhoton + " && " + self.gamma  +" && " + self.noLeptons  + "&&" + self.filters )
  
   def gammasData(self):

       return self.brackets(self.centralPhoton + " && " + self.gamma + " && " + self.triggerPhoton  +" && " + self.noLeptons+ "&&" + self.filters )
   
   def leps(self):

          return self.brackets(self.filters + " && " +self.Zmasswindow  + " && "+ self.twoGoodLeptons ) 

      


