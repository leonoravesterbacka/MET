import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, THStack, TCanvas


class Sample:
   'Common base class for all Samples'

   def __init__(self, name, location, xsection, isdata, doee, dokfactorweight, iszjets):
      self.name = name
      self.location = location
      self.xSection = xsection
      self.isData = isdata
      self.isee = doee
      self.dokfactorWeight = dokfactorweight
      self.isZjets = iszjets
      #self.tfile = TFile(self.location+self.name+'/METtree.root')
      self.tfile = TFile(self.location+self.name+'/jes_ALLMETtree.root')
      #self.tfile = TFile(self.location+self.name+'/jes_jes_BCDMETtree.root')
      self.ttree = self.tfile.Get('METtree')
      print self.name      
      self.puWeight  = "1.0"
      if not self.isData:
          gw = 0.
          for i in self.ttree:
              gw = abs(i.genWeight)
              if gw: break
          
          self.count = self.tfile.Get('SumGenWeights').GetBinContent(1)/abs(gw)
      else:
          self.count = self.tfile.Get('Count').GetBinContent(1)
      self.lumWeight =  1.0
      if not self.isData:
        self.lumWeight = self.xSection / self.count
        print 'name ',self.name
        print 'xsec ',self.xSection
        print 'count ',self.count
        print 'lumweight ' ,self.lumWeight

   def printSample(self):
      print "#################################"
      print "Sample Name: ", self.name
      print "Sample Location: ", self.location
      print "Sample XSection: ", self.xSection
      print "Sample IsData: ", self.isData
      print "Sample LumWeight: ", self.lumWeight
      print "#################################"


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, doNPV):
 
      if(xmin == xmax):
        h = TH1F(name, "", len(nbin)-1, array('d', nbin))
        ylabel = "# events"
      else:
        h = TH1F(name, "", nbin, xmin, xmax)
        bw = int((xmax-xmin)/nbin)
        ylabel = "Events / " + str(bw) + " GeV"
      h.Sumw2()
      h.GetXaxis().SetTitle(xlabel)
      h.GetYaxis().SetTitle(ylabel)

      addCut = "1."
      addPrescale = "1." 
      if self.isData:
          addDataFilters = "&&(  (Flag_eeBadScFilter == 1  ))"
          if self.isZjets:
              #addTriggers = "&&( ((HLT_DoubleMu == 1 || HLT_SingleMu == 1 || HLT_HighPTMuNonIso == 1) && ((lep_tightId[0] ==1 && lep_mediumMuonId[0] ==1 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_tightId[1] ==1 && lep_mediumMuonId[1] ==1 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 )  )) || ((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && ( (lep_tightId[0] ==2 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0 && lep_tightId[1] ==2 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0) )))"
              addTriggers = "&&( ((HLT_DoubleMu == 1 || HLT_SingleMu == 1 || HLT_HighPTMuNonIso == 1) && ((lep_tightId[0] ==1 && lep_mediumMuonId[0] ==1 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_tightId[1] ==1 && lep_mediumMuonId[1] ==1 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 )  )) || ((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && ( (lep_tightId[0] ==3 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0 && lep_tightId[1] ==3 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0) )))"
          else:   
              addTriggers = "&&( HLT_Photon120 ==1 || HLT_Photon165 == 1 || HLT_Photon90 == 1 || HLT_Photon75 == 1 || HLT_Photon50 == 1 || HLT_Photon30 == 1   ) " 
              
              #addPrescale = "(( tr165 == 1  && tr90 ==0 && tr75==0 && tr50==0 && tr30==0  ) + (ps120)*( tr120 == 1 &&  tr165 == 0 && tr75 == 0 && tr50==0 &&  tr30==0  ) + (ps90)*( tr90 == 1 && tr120 == 0 &&  tr75 == 1  &&tr50 == 1  && tr30==1  ) +   (ps75)*( tr75 == 1 && tr120 == 0  &&  tr50 == 0 && tr30 == 0 ) +   (ps50)*( tr50 == 1 && tr120 == 0 && tr90 == 0&& tr75 == 0 && tr30 == 0) ) " 
              #addPrescale = "(( HLT_Photon165 == 1  && HLT_Photon90==0 && HLT_Photon50==0 && HLT_Photon30==0  ) + (HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon120 == 1 &&  HLT_Photon165 == 0 && HLT_Photon50==0 &&  HLT_Photon30==0  ) + (HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon90 == 1 && HLT_Photon120 == 0 &&  HLT_Photon50 == 1  && HLT_Photon30==1  ) +   (HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon50 == 1 && HLT_Photon120 == 0 && HLT_Photon90 == 0&& HLT_Photon30 == 0) ) " 
              #addPrescale = "(( HLT_Photon165 == 1  && HLT_Photon90==0 && HLT_Photon75==0 && HLT_Photon50==0 && HLT_Photon30==0  ) + (HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon120 == 1 &&  HLT_Photon165 == 0 && HLT_Photon75 == 0 && HLT_Photon50==0 &&  HLT_Photon30==0  ) + (HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon90 == 1 && HLT_Photon120 == 0 &&  HLT_Photon75 == 1  &&  HLT_Photon50 == 1  && HLT_Photon30==1  ) +   (HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale)*(0.95)*( HLT_Photon75 == 1 && HLT_Photon120 == 0  &&  HLT_Photon50 == 0 && HLT_Photon30 == 0 ) +   (HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon50 == 1 && HLT_Photon120 == 0 && HLT_Photon90 == 0&& HLT_Photon75 == 0 && HLT_Photon30 == 0) ) " 
              addPrescale = "(( HLT_Photon165 == 1  && HLT_Photon90==0 && HLT_Photon75==0 && HLT_Photon50==0 && HLT_Photon30==0  ) + (HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon120 == 1 &&  HLT_Photon165 == 0 && HLT_Photon75 == 0 && HLT_Photon50==0 &&  HLT_Photon30==0  ) + (HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon90 == 1 && HLT_Photon120 == 0 &&  HLT_Photon75 == 1  &&  HLT_Photon50 == 1  && HLT_Photon30==1  ) +   (HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon75 == 1 && HLT_Photon120 == 0  &&  HLT_Photon50 == 0 && HLT_Photon30 == 0 ) +   (HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale)*( HLT_Photon50 == 1 && HLT_Photon120 == 0 && HLT_Photon90 == 0&& HLT_Photon75 == 0 && HLT_Photon30 == 0) ) " 
              #addPrescale = "(( HLT_Photon165 == 1  && HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v==0 && HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v==0 && HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v==0 && HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v==0  ) + (HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale)*( HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v == 1 &&  HLT_Photon165 == 0 && HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v == 0 && HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v==0 &&  HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v==0  ) + (HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale)*( HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v == 1 && HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v == 0 &&  HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v == 1  &&  HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v == 1  && HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v==1  ) +   (HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale)*( HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v == 1 && HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v == 0  &&  HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v == 0 && HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v == 0 ) +   (HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale)*( HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v == 1 && HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v == 0 && HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v == 0&& HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v == 0 && HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v == 0) ) " 
          cut = "("+ cut + addTriggers + addDataFilters+ ")" + "* (" + addPrescale +")" 

      if(self.isData == 0):
          addCut = "1."
          if doNPV: 
              #addCut = "1"
              addCut = "puWeightNPV"
          else:
              addCut = "puWeightNTI"
          if self.dokfactorWeight:
              #cut =  cut  + "* ( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )" 
              cut =  cut  + "*( puWeightGJets )"+  "*( photonEff )"  +  "*( kfactorWeight )" + "* ( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )" 
          else:
              if self.isee:
                  #cut =  cut  + "* ( " + addCut + " )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && (lep_tightId[0] ==1)*lep1eff )"+ "*((lep_tightId[1] ==1) *lep2eff )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "*((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && (lep_tightId[0] ==3 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0))"+ "*((lep_tightId[1] ==3 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0))" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && (lep_tightId[0] ==3 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0) )"+ "*((lep_tightId[1] ==3 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0) )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && (lep_tightId[0] ==3 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0))"+ "*((lep_tightId[1] ==3 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0) )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleEG  == 1 || HLT_SingleEle == 1 || HLT_HighPTEleNonIso ==1 ) && (lep_tightId[0] ==3 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0)*lep1eff )"+ "*((lep_tightId[1] ==3 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0) *lep2eff )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
              else:
                  #cut =  cut  + "*((HLT_DoubleMu == 1 || HLT_SingleMu == 1 || HLT_HighPTMuNonIso == 1) && (lep_tightId[0] ==1 && lep_mediumMuonId[0] ==1 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0 && (abs(lep_eta[0]) < 2.4)))"+ "*((lep_tightId[1] ==1 && lep_mediumMuonId[1] ==1 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0 && (abs(lep_eta[1]) < 2.4)) )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleMu == 1 || HLT_SingleMu == 1 || HLT_HighPTMuNonIso == 1) && (lep_tightId[0] ==1 && lep_mediumMuonId[0] ==1 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0 && (abs(lep_eta[0]) < 2.4)))"+ "*((lep_tightId[1] ==1 && lep_mediumMuonId[1] ==1 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0 && (abs(lep_eta[1]) < 2.4)) )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleMu == 1 || HLT_SingleMu == 1 || HLT_HighPTMuNonIso == 1) && (lep_tightId[0] ==1 && lep_mediumMuonId[0] ==1 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0 && (abs(lep_eta[0]) < 2.4)))"+ "*((lep_tightId[1] ==1 && lep_mediumMuonId[1] ==1 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0 && (abs(lep_eta[1]) < 2.4)) )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  cut =  cut  + "* ( " + addCut + " )" + "*((HLT_DoubleMu == 1 || HLT_SingleMu == 1 || HLT_HighPTMuNonIso == 1) && (lep_tightId[0] ==1 && lep_mediumMuonId[0] ==1 && lep_dxy[0] < 0.05 && lep_dz[0]< 0.1 && lep_convVeto[0] ==1 && lep_lostHits[0] == 0 && (abs(lep_eta[0]) < 2.4))*lep1eff )"+ "*((lep_tightId[1] ==1 && lep_mediumMuonId[1] ==1 && lep_dxy[1] < 0.05 && lep_dz[1]< 0.1 && lep_convVeto[1] ==1 && lep_lostHits[1] == 0 && (abs(lep_eta[1]) < 2.4)) *lep2eff )" + "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
      self.ttree.Project(name, var, cut, options)
      #print cut
      return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight/abs(genWeight) * " + self.puWeight + " )" 
     
     self.ttree.Project(name, var, cut, options) 
     return h

class Block:
   'Common base class for all Sample Blocks'

   def __init__(self, name, label, color, isdata):
      self.name  = name
      self.color = color
      self.isData = isdata
      self.label = label
      self.samples = []

   def printBlock(self):

      print "####################"
      print "Block Name: ", self.name
      print "Block Color: ", self.color
      print "Block IsData: ", self.isData
      print "####################"
      print "This block contains the following Samples"

      for l in self.samples:
        l.printSample()
     

   def addSample(self, s):
      self.samples.append(s)

   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, rereco):

     if(xmin == xmax):
       h = TH1F(name, "", len(nbin)-1, array('d', nbin))
       ylabel = "# events"
     else:
       h = TH1F(name, "", nbin, xmin, xmax)
       bw = int((xmax-xmin)/nbin)
       ylabel = "Events / " + str(bw) + " GeV"
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)

     for s in self.samples:
       AuxName = "auxT1_sample" + s.name
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, rereco)
       h.Add(haux)
       del haux


     h.SetLineColor(self.color)
     h.SetMarkerColor(self.color)
     h.SetTitle(self.label)

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for s in self.samples:
     
       AuxName = "auxT2_block" + s.name
       haux = s.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   

       

class Tree:
   'Common base class for a physics meaningful tree'

   def __init__(self, fileName, name, isdata):
      self.name  = name
      self.isData = isdata
      self.blocks = []
      self.parseFileName(fileName)

   def parseFileName(self, fileName):

      f = open(fileName)

      for l in f.readlines():

        if (l[0] == "#" or len(l) < 2):
          continue

        splitedLine = str.split(l)
        block       = splitedLine[0]
        theColor    = splitedLine[1]
        name        = splitedLine[2]
        label       = splitedLine[3]
        location    = splitedLine[4]
        xsection    = float(splitedLine[5])
        isdata =  int(splitedLine[6])
        doboson =  int(splitedLine[7])
        dokfactor =  int(splitedLine[8])
        iszjets =  int(splitedLine[9])

        color = 0
        plusposition = theColor.find("+")
        if(plusposition == -1):
          color = eval(theColor)
        else:
          color = eval(theColor[0:plusposition])
          color = color + int(theColor[plusposition+1:len(theColor)])

        sample = Sample(name, location,  xsection, isdata, doboson, dokfactor, iszjets)
        coincidentBlock = [l for l in self.blocks if l.name == block]

        if(coincidentBlock == []):

          newBlock = Block(block, label, color, isdata)
          newBlock.addSample(sample)
          self.addBlock(newBlock)

        else:

          coincidentBlock[0].addSample(sample)





   def printTree(self):

      print "######"
      print "Tree Name: ", self.name
      print "Tree IsData: ", self.isData
      print "######"
      print "This Tree contains the following Blocks"

      for l in self.blocks:
        l.printBlock()
     

   def addBlock(self, b):
      self.blocks.append(b)



   def getYields(self, lumi, var, xmin, xmax, cut):
  
      h = self.getTH1F(lumi, "yields", var, 1, xmin, xmax, cut, "", "", rereco)
      nbinmin = h.FindBin(xmin)
      nbinmax = h.FindBin(xmax)
      error = r.Double()
      value = h.IntegralAndError(nbinmin, nbinmax, error)
      y = [value, error]
      
      del h
      return y

   def getStack(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
   
     hs = THStack(name, "")
     for b in self.blocks:
     
       AuxName = "auxStack_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, rereco)
       haux.SetFillColor(b.color)
       hs.Add(haux)
       del haux


     can_aux = TCanvas("can_aux")
     can_aux.cd()
     hs.Draw()
     del can_aux

     hs.GetXaxis().SetTitle(xlabel)
     b = int((xmax-xmin)/nbin)
     ylabel = "Events / " + str(b) + " GeV"
     hs.GetYaxis().SetTitle(ylabel)

     return hs   


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel, rereco):
     
     if(xmin == xmax):
       _nbins = len(nbin)-1
       _arr = array('d', nbin)
       h = TH1F(name, "", _nbins, _arr)
       if len(nbin) > 50 :  # this option is for when you want underflow bins, make sure to modify the nbins if you want underflow bins for a plot with fewer bins 
           _newarr = array('d', [ 2*_arr[0]-_arr[1] ]) +_arr +  array('d', [ 2*_arr[-1]-_arr[-2] ]) 
           h_of = TH1F(name+'_of', "", _nbins+2, _newarr)
       else:
           _newarr = _arr + array('d', [ 2*_arr[-1]-_arr[-2] ]) 
           h_of = TH1F(name+'_of', "", _nbins+1, _newarr)
       ylabel = "# events"
     else:
         h = TH1F(name, "", nbin, xmin, xmax)
         bw = int((xmax-xmin)/nbin)
         ylabel = "Events / " + str(bw) + " GeV"
         h_of = TH1F(name+'_of', '', nbin+1, xmin, xmax+bw)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
       AuxName = "auxh1_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, rereco)
       h.Add(haux)
       del haux

     if options == 'noOF':
         for _bin in range(0, h_of.GetNbinsX()):
            h_of.SetBinContent(_bin, h.GetBinContent(_bin))
            h_of.SetBinError(_bin, h.GetBinError(_bin))        
     else:
         if ((xmin == xmax) and len(nbin) > 50):    
             for _bin in range(0, h_of.GetNbinsX()+1):
                 h_of.SetBinContent(_bin+1, h.GetBinContent(_bin))
                 h_of.SetBinError  (_bin+1, h.GetBinError  (_bin)) 
         else:
             for _bin in range(1, h_of.GetNbinsX()+1):
                 h_of.SetBinContent(_bin, h.GetBinContent(_bin))
                 h_of.SetBinError(_bin, h.GetBinError(_bin))        


     return h_of
     del h_of
     del h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   

