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
      self.tfile = TFile(self.location+self.name+'/jes_jes_METtree.root')
      self.ttree = self.tfile.Get('METtree')
      
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


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
 
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
      addCut2 = ""
      
      if self.isData:
          if self.isZjets:
              addTriggers = "&&( HLT_DoubleMu == 1 || HLT_DoubleEG  == 1 || HLT_SingleMu == 1 || HLT_SingleEle == 1 )"
          else:   
              addCut = "prescaleWeight" 
              addTriggers = "" 
          cut = cut + addTriggers + "* (" + addCut +") "

      if(self.isData == 0):
          #addCut = "1"
          addCut = "puWeight"
          addCut2 = "1."
          if self.dokfactorWeight:
              cut =  cut  + "* ( " + addCut + " )" +  "*( kfactorWeight )" + "* ( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )" 
          else:
              if self.isee:
                  cut =  cut  + "* ( " + addCut + " )" +  "*0.83320698315*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
                  #cut =  cut  + "* ( " + addCut + " )" +  "*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )"
              else:
                  cut =  cut  + "* ( " + addCut + " )" +  "* 0.86905641268*( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )" 
                  #cut =  cut  + "* ( " + addCut + " )" +  "* ( " + str(self.lumWeight*lumi)  +  " )" + "* ( " + "genWeight/abs(genWeight) " +  " )" 
      #print "cuts ", cut
      self.ttree.Project(name, var, cut, options)
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

   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):

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
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
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
  
      h = self.getTH1F(lumi, "yields", var, 1, xmin, xmax, cut, "", "")
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
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
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


   def getTH1F(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
   
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
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
       h.Add(haux)
       del haux

     if ((xmin == xmax) and len(nbin) > 41):    
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










