from ROOT import TCanvas, TLegend, TPad, TLine,  TGraphAsymmErrors,  TLatex, TH1F, THStack, TGraphErrors, TLine, TPaveStats, TGraph, TArrow
import ROOT as r
import os, copy, math, array

class Canvas:
   'Common base class for all Samples'

   def __init__(self, name, format, x1, y1, x2, y2):
      self.name = name
      self.format = format
      self.plotNames = [name + "." + i for i in format.split(',')]
      self.myCanvas = TCanvas(name, name)
      self.ToDraw = []
      self.orderForLegend = []
      self.histos = []
      self.lines = []
      self.arrows= []
      self.latexs= []
      self.bands = []
      self.options = []
      self.labels = []      
      self.labelsOption = []
      self.myLegend = TLegend(x1, y1, x2, y2)
      self.myLegend.SetFillColor(r.kWhite)
      self.myLegend.SetTextFont(42)
      self.myLegend.SetTextSize(0.05)
      self.myLegend.SetLineWidth(0)
      self.myLegend.SetBorderSize(0)


   def banner(self, isData, lumi, isLog, event, isnVert, isSig1jet, run):
    
      latex = TLatex()                                
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.08);
      if isnVert:
          latex.DrawLatex(0.295, 0.93, "#bf{CMS}")
      else:
          latex.DrawLatex(0.23, 0.93, "#bf{CMS}")

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.06);
 
      if(isData):
          if isnVert:
              latexb.DrawLatex(0.485, 0.93, "#it{Preliminary}")
          else:
              latexb.DrawLatex(0.42, 0.93, "#it{Preliminary}")
      else:
        latexb.DrawLatex(0.42, 0.93, "#it{Simulation}")

      text_lumi = str(lumi)+" fb^{-1} (13 TeV, "+ run + ")"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.06);
      latexc.DrawLatex(0.91, 0.93, text_lumi)          

      latexd = TLatex()
      latexd.SetNDC();
      latexd.SetTextAngle(90);
      latexd.SetTextColor(r.kBlack);
      latexd.SetTextFont(42);
      latexd.SetTextAlign(31);
      latexd.SetTextSize(0.06);
      if event == 0:
          latexd.DrawLatex(0.059, 0.93, "Events")          
      else:
          latexd.DrawLatex(0.059, 0.93, "Events / " +str(event) + " GeV")          

      latexe = TLatex()
      latexe.SetNDC();
      latexe.SetTextColor(r.kBlack);
      latexe.SetTextFont(42);
      latexe.SetTextAlign(31);
      latexe.SetTextSize(0.045);         
      if isSig1jet:
          latexe.DrawLatex(0.86, 0.55, "p_{T}^{jet1} > 50 GeV")   


   def banner2(self, lumi, chisquare, mean, title, event, isNVert, isSig1jet):
    
      latex = TLatex()
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.06);
      if isNVert:
          latex.DrawLatex(0.29, 0.83, "#bf{CMS}")
      else:
          latex.DrawLatex(0.26, 0.93, "#bf{CMS}")

      latexd = TLatex()
      latexd.SetNDC();
      latexd.SetTextAngle(90);
      latexd.SetTextColor(r.kBlack);
      latexd.SetTextFont(42);
      latexd.SetTextAlign(31);
      latexd.SetTextSize(0.06);
      if event == 0:
          latexd.DrawLatex(0.047, 0.93, "Events")          
      if event == '':
          latexd.DrawLatex(0.047, 0.93, "")          
      else:
          latexd.DrawLatex(0.059, 0.93, "Events / " +str(event) + " GeV")        
      


      latexe = TLatex()
      latexe.SetNDC();
      latexe.SetTextColor(r.kBlack);
      latexe.SetTextFont(42);
      latexe.SetTextAlign(31);
      latexe.SetTextSize(0.045);         
      if isSig1jet:
          latexe.DrawLatex(0.86, 0.55, "p_{T}^{jet1} > 50 GeV")   
      
      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.05);
 
      if isNVert:
          latexb.DrawLatex(0.52, 0.83, "#it{Preliminary}")
      else:
          latexb.DrawLatex(0.49, 0.93, "#it{Preliminary}")

      text_lumi = str(lumi)+" fb^{-1} (Run B-H)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      latexc.DrawLatex(0.91, 0.93, text_lumi)
      latexc.DrawLatex(0.46, 0.8,  mean)
      if mean == "":
          latexc.DrawLatex(0.90, 0.04,  mean)
      else:
          latexc.DrawLatex(0.39, 0.7,  "#chi^{2} = %1.f " %(chisquare))

      latexd = TLatex()                                                  
      latexd.SetNDC();
      latexd.SetTextAngle(0);
      latexd.SetTextColor(r.kBlack);
      latexd.SetTextFont(42);
      latexd.SetTextAlign(31);
      latexd.SetTextSize(0.05);
      latexd.DrawLatex(0.89, 0.01,  title)




   def addBand(self, x1, y1, x2, y2, color, opacity):

      grshade = TGraph(4)
      grshade.SetPoint(0,x1,y1)
      grshade.SetPoint(1,x2,y1)
      grshade.SetPoint(2,x2,y2)
      grshade.SetPoint(3,x1,y2)
      #grshade.SetFillStyle(3001)
      grshade.SetFillColorAlpha(color, opacity)
      self.bands.append(grshade)

   def addLine(self, x1, y1, x2, y2, color, thickness = 0.):
      line = TLine(x1,y1,x2,y2)
      line.SetLineColor(color)
      if thickness:
          line.SetLineWidth(thickness)
      self.lines.append(line)

   def addArrow(self, x1, y1, x2, y2, color, option, thickness = 0.):
      arrow = TArrow(x1,y1,x2,y2, 0.05, option)
      arrow.SetLineColor(color)
      if thickness:
          arrow.SetLineWidth(thickness)
      self.arrows.append(arrow)

   def addLatex(self, x1, y1, text, font=42, size = 0.02):
      lat = [x1, y1, text, font, size]
      self.latexs.append(lat)

 
   def addHisto(self, h, option, label, labelOption, color, ToDraw, orderForLegend):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)

   def addGraph(self, h, option, label, labelOption, color, ToDraw, orderForLegend):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)


   def addStack(self, h, option, ToDraw, orderForLegend):

      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)

      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for h_c in h.GetHists():
          self.addHisto(h_c, "H", h_c.GetTitle(), "F", "", 0, legendCounter)
          legendCounter = legendCounter + 1
       

 
   def makeLegend(self, log):
      #if log:  
     for i in range(0, len(self.histos)):
         for j in range(0, len(self.orderForLegend)):
             if(self.orderForLegend[j] != -1 and self.orderForLegend[j] == i):
                 self.myLegend.AddEntry(self.histos[j], self.labels[j], self.labelsOption[j])
      #else:
      #    self.myLegend.AddEntry(self.histos[4], self.labels[4], self.labelsOption[4])
      #    self.myLegend.AddEntry(self.histos[3], self.labels[3], self.labelsOption[3])


          

   def ensurePath(self, _path):
      d = os.path.dirname(_path)
      if not os.path.exists(d):
         os.makedirs(d)

   def saveRatio(self, legend, isData, log, lumi, hdata, hMC, hjecUp, hjecDown ,  hjerUp, hjerDown ,hunclUp, hunclDown, title,  option, run = '2016',  r_ymin=0, r_ymax=2):
      
      events = 5  
      setUpAxis = 0 
      setLowAxis = 1 
      log = 1
      doAllErrors = 1
      isNVert = 0
      puppi = 0
      statOnly = 0
      tightRatio = 0
      isSig1jet = 0
      lowAxis = 0
      if option == "test": 
          doAllErrors = 0
          statOnly = 1
          isNVert = 1
          tightRatio = 1   
      if option == "test2": 
          doAllErrors = 0
          statOnly = 1
          isNVert = 1
          tightRatio = 1   
          lowAxis = 1
          setLowAxis = 0
          setUpAxis = 0        
      if option == 'nvert':
          log =0
          events = 0
          doAllErrors = 0
          statOnly = 1
          isNVert = 1
          tightRatio = 1    
      if option == 'sig0':
          log =1
          events = 2
          doAllErrors = 0  
          statOnly = 1     
      if option == 'sig1':
          log =1
          events = 2
          doAllErrors = 0  
          statOnly = 1  
          isSig1jet = 1      
      if option == 'chi':
          log =0
          events = 0
          doAllErrors = 0  
          statOnly = 1
      if option == 'qt':
          setLowAxis = 1
          setUpAxis = 1
          doAllErrors = 0  
          statOnly = 1
          tightRatio = 1
      if option == 'mass':
          fixAxis = 0
          doAllErrors = 0
          events = 2          
          statOnly = 1
          tightRatio = 1
      if option == 'puppi':
          puppi = 0
          doAllErrors = 1
          setLowAxis = 1          
      
      self.myCanvas.cd()
      pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0) 
      pad1.SetBottomMargin(0.01)
      pad1.Draw()
      pad2 = TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
      pad2.SetTopMargin(0.1);
      pad2.SetBottomMargin(0.3);
      pad2.Draw();

      pad1.cd()
      if(log):
          pad1.SetLogy(1)

      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):
              if lowAxis: 
                  self.histos[i].SetMinimum(0.00001)
              if setLowAxis:
                  self.histos[i].SetMinimum(10.)
              if setUpAxis:    
                  self.histos[i].SetMaximum(hMC.Integral()*10)
              self.histos[i].Draw(self.options[i])

      if(legend):
          self.makeLegend(log)
          self.myLegend.Draw()

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextSize(latex[-1])
          lat.SetTextFont(latex[-2])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      
      ratio = copy.deepcopy(hdata.Clone("ratio"))
      ratio.Divide(hMC)

      if tightRatio:
          ratio.GetYaxis().SetRangeUser(0.5, 1.5);
      else:
          ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
      if option == "test":
          ratio.GetYaxis().SetTitle("postfix/prefix")
      else:
          ratio.GetYaxis().SetTitle("Data / MC")
      ratio.GetYaxis().CenterTitle();
      ratio.GetYaxis().SetLabelSize(0.12);
      ratio.GetXaxis().SetLabelSize(0.12);
      ratio.GetXaxis().SetTitleOffset(0.91);
      ratio.GetYaxis().SetNdivisions(4);
      ratio.GetYaxis().SetTitleSize(0.13);
      ratio.GetXaxis().SetTitleSize(0.135);
      ratio.GetYaxis().SetTitleOffset(0.41);
      ratio.SetMarkerSize(0.6*ratio.GetMarkerSize());
      ratio.GetXaxis().SetTitle(title)
      
      den1tottot = copy.deepcopy(hMC.Clone("bkgden1"))
      den2tottot = copy.deepcopy(hMC.Clone("bkgden2"))
      	                                                                                                                                                                                 
      nvar = hMC.GetNbinsX()                                                                                                                                                           
      x = array.array('f', range(0, hMC.GetNbinsX()))
      y = array.array('f', range(0, hMC.GetNbinsX()))
      exl = array.array('f', range(0, hMC.GetNbinsX()))
      eyltottot = array.array('f', range(0, hMC.GetNbinsX()))
      exh = array.array('f', range(0, hMC.GetNbinsX()))
      eyhtottot = array.array('f', range(0, hMC.GetNbinsX()))
      ratiouptottot = copy.deepcopy(hMC.Clone("ratiouptottot"))
      ratiodowntottot = copy.deepcopy(hMC.Clone("ratiodowntottot"))      
      ymax = 2.                                                                                                                                                              
      #make tot tot errors 
      for km in range(0, hMC.GetNbinsX()):
          conte1tottot =  math.sqrt(hMC.GetBinError (km) * hMC.GetBinError (km) + (hjecUp.GetBinContent (km) - hMC.GetBinContent   (km))*(hjecUp.GetBinContent (km) -  hMC.GetBinContent (km)) + (hunclUp.GetBinContent (km) - hMC.GetBinContent   (km))*(hunclUp.GetBinContent (km) -  hMC.GetBinContent (km)) + (hjerUp.GetBinContent (km) - hMC.GetBinContent   (km))*(hjerUp.GetBinContent (km) -  hMC.GetBinContent (km))    );       
          conte2tottot =  math.sqrt(hMC.GetBinError (km) * hMC.GetBinError (km) + (hMC.GetBinContent (km) - hjecDown.GetBinContent (km))*(hMC.GetBinContent (km) -  hjecDown.GetBinContent (km)) + (hunclDown.GetBinContent (km) - hMC.GetBinContent   (km))*(hunclDown.GetBinContent (km) -  hMC.GetBinContent (km))  + (hjerDown.GetBinContent (km) - hMC.GetBinContent   (km))*(hjerDown.GetBinContent (km) -  hMC.GetBinContent (km)));   
          if conte1tottot > conte2tottot:
              den1tottot.SetBinContent (km, hMC.GetBinContent (km) + conte1tottot);                                                                                                                           
              den2tottot.SetBinContent (km, hMC.GetBinContent (km) - conte1tottot);                                                                                                                           
          else:
              den1tottot.SetBinContent (km, hMC.GetBinContent (km) + conte2tottot);  
              den2tottot.SetBinContent (km, hMC.GetBinContent (km) - conte2tottot);
          ymax = hMC.GetBinContent(km) + conte1tottot;           
          exl[km] = hMC.GetBinWidth (km) / 2;                                                                                                                                 
          exh[km] = hMC.GetBinWidth (km) / 2;                                                                                                                                 
          eyltottot[km] = conte2tottot;                                                                                                                                                           
          eyhtottot[km] = conte1tottot;                                                                                                                         
      
      ratiouptottot.Divide(den1tottot);                                                                                                                                                         
      ratiodowntottot.Divide(den2tottot);                   
      ratiodata = copy.deepcopy(hdata.Clone("ratiodata"))
      ratiodata.Divide (hMC);                                                                                                                                                    
                                                                                                                                                                                         
      for km in range(0, ratiodata.GetNbinsX()):
          if (ratiodata.GetBinContent (km) > ymax):
              ymax = ratiodata.GetBinContent (km) + ratiodata.GetBinError (km);                                                                                                              
          x[km] = ratiodata.GetBinCenter (km);                                                                                                                                          
          y[km] = 1;	                                                                                                                                                             
          exl[km] = ratiodata.GetBinWidth (km) / 2;                                                                                                                                     
          exh[km] = ratiodata.GetBinWidth (km) / 2;                                                                                                                                     
                                                                                                                                                                                         
          if (ratiouptottot.GetBinContent (km) != 0):
              eyhtottot[km] = (1. / ratiouptottot.GetBinContent (km) - 1)*ratiodata.GetBinContent (km);                                                                    
          else:                                                                                                                                                                       
              eyhtottot[km] = 0;                                                                                                                                                                 
                                                                                                                                                                                         
          if (ratiodowntottot.GetBinContent (km) != 0):
              eyltottot[km] = (1 - 1. / ratiodowntottot.GetBinContent (km))*ratiodata.GetBinContent (km);                                                              
          else:                                                                                                                                                                       
              eyltottot[km] = 0.                                                                                                                        
      errtottot = TGraphAsymmErrors(nvar, x, y, exl, exh, eyltottot, eyhtottot);                 



      #make tot errors
      den1tot = copy.deepcopy(hMC.Clone("bkgden1"))
      den2tot = copy.deepcopy(hMC.Clone("bkgden2"))
      	                                                           
      nvar = hMC.GetNbinsX()                                       
      x = array.array('f', range(0, hMC.GetNbinsX()))
      y = array.array('f', range(0, hMC.GetNbinsX()))
      exl = array.array('f', range(0, hMC.GetNbinsX()))
      eyltot = array.array('f', range(0, hMC.GetNbinsX()))
      exh = array.array('f', range(0, hMC.GetNbinsX()))
      eyhtot = array.array('f', range(0, hMC.GetNbinsX()))
      ratiouptot = copy.deepcopy(hMC.Clone("ratiouptot"))
      ratiodowntot = copy.deepcopy(hMC.Clone("ratiodowntot"))      
      
      
      for km in range(0, hMC.GetNbinsX()):
          conte1tot =  math.sqrt(hMC.GetBinError (km) * hMC.GetBinError (km) + (hjecUp.GetBinContent (km) - hMC.GetBinContent   (km))*(hjecUp.GetBinContent (km) -  hMC.GetBinContent (km)) + (hjerUp.GetBinContent (km) - hMC.GetBinContent   (km))*(hjerUp.GetBinContent (km) -  hMC.GetBinContent (km)));       
          conte2tot =  math.sqrt(hMC.GetBinError (km) * hMC.GetBinError (km) + (hMC.GetBinContent (km) - hjecDown.GetBinContent (km))*(hMC.GetBinContent (km) -  hjecDown.GetBinContent (km)) + (hjerDown.GetBinContent (km) - hMC.GetBinContent   (km))*(hjerDown.GetBinContent (km) -  hMC.GetBinContent (km)));   
          if conte1tot > conte2tot:
              den1tot.SetBinContent (km, hMC.GetBinContent (km) + conte1tot);                                                                                           
              den2tot.SetBinContent (km, hMC.GetBinContent (km) - conte1tot);                                                                                           
          else:
              den1tot.SetBinContent (km, hMC.GetBinContent (km) + conte2tot);  
              den2tot.SetBinContent (km, hMC.GetBinContent (km) - conte2tot);
          ymax = hMC.GetBinContent(km) + conte1tot;           
          exl[km] = hMC.GetBinWidth (km) / 2;                                                                                                                                 
          exh[km] = hMC.GetBinWidth (km) / 2;                                                                                                                                 
          eyltot[km] = conte2tot;                                                                                                                                                           
          eyhtot[km] = conte1tot;                                                                                   
      
      ratiouptot.Divide(den1tot);                                                                                                                                                         
      ratiodowntot.Divide(den2tot);                   
      ratiodata = copy.deepcopy(hdata.Clone("ratiodata"))
      ratiodata.Divide (hMC);                                                                                                                                                    
                                                                                                                                                                                         
      for km in range(0, ratiodata.GetNbinsX()):
          if (ratiodata.GetBinContent (km) > ymax):
              ymax = ratiodata.GetBinContent (km) + ratiodata.GetBinError (km);                                                                                                              
          x[km] = ratiodata.GetBinCenter (km);                                                                                                                                          
          y[km] = 1;	                                                                                                                                                             
          exl[km] = ratiodata.GetBinWidth (km) / 2;                                                                                                                                     
          exh[km] = ratiodata.GetBinWidth (km) / 2;                                                                                                                                     
                                                                                                                                                                                         
          if (ratiouptot.GetBinContent (km) != 0):
              eyhtot[km] = (1. / ratiouptot.GetBinContent (km) - 1)*ratiodata.GetBinContent (km);                                                                                                   
          else:                                                                                                                                                                       
              eyhtot[km] = 0;                                                                                                                                                                 
                                                                                                                                                                                         
          if (ratiodowntot.GetBinContent (km) != 0):
              eyltot[km] = (1 - 1. / ratiodowntot.GetBinContent (km))*ratiodata.GetBinContent (km);                                                                                                   
          else:                                                                                                                                                                       
              eyltot[km] = 0.                                                                                                                        
      errtot = TGraphAsymmErrors(nvar, x, y, exl, exh, eyltot, eyhtot);                                                                                 
      
      #make some jec + stat errors     
      den1 = copy.deepcopy(hMC.Clone("bkgden1"))
      den2 = copy.deepcopy(hMC.Clone("bkgden2"))
      	                                                                                                                                                                                 
      nvar = hMC.GetNbinsX()                                                                                                                                                           
      eyl = array.array('f', range(0, hMC.GetNbinsX()))
      eyh = array.array('f', range(0, hMC.GetNbinsX()))
      ratioup = copy.deepcopy(hMC.Clone("ratioup"))
      ratiodown = copy.deepcopy(hMC.Clone("ratiodown"))
      
      for km in range(0, hMC.GetNbinsX()):
          conte1 =  math.sqrt(hMC.GetBinError (km) * hMC.GetBinError (km) + (hjecUp.GetBinContent (km) - hMC.GetBinContent   (km))*(hjecUp.GetBinContent (km) -  hMC.GetBinContent (km)));       
          conte2 =  math.sqrt(hMC.GetBinError (km) * hMC.GetBinError (km) + (hMC.GetBinContent (km) - hjecDown.GetBinContent (km))*(hMC.GetBinContent (km) -  hjecDown.GetBinContent (km)));   
          if conte1 > conte2:
              den1.SetBinContent (km, hMC.GetBinContent (km) + conte1);                                                                                                                           
              den2.SetBinContent (km, hMC.GetBinContent (km) - conte1);                                                                                                                           
          else:
              den1.SetBinContent (km, hMC.GetBinContent (km) + conte2);  
              den2.SetBinContent (km, hMC.GetBinContent (km) - conte2);
          ymax = hMC.GetBinContent(km) + conte1;           
          eyl[km] = conte2;                                                                                                                                                           
          eyh[km] = conte1;                                                                                                                                                                        
      
      ratioup.Divide(den1);                                                                                                                                                         
      ratiodown.Divide(den2);                   
                                                                                                                                                                                         
      for km in range(0, ratiodata.GetNbinsX()):
          if (ratioup.GetBinContent (km) != 0):
              eyh[km] = (1. / ratioup.GetBinContent (km) - 1)*ratiodata.GetBinContent (km);                                                                                                   
          else:                                                                                                                                                                       
              eyh[km] = 0;                                                                                                                                                                 
                                                                                                                                                                                         
          if (ratiodown.GetBinContent (km) != 0):
              eyl[km] = (1 - 1. / ratiodown.GetBinContent (km))*ratiodata.GetBinContent (km);                                                                                                   
          else:                                                                                                                                                                       
              eyl[km] = 0.                                                                                                                        
      err = TGraphAsymmErrors(nvar, x, y, exl, exh, eyl, eyh);                                                                                                                                         

      #make some stat errors
      dens1 = copy.deepcopy(hMC.Clone("bkgdens1"))
      dens2 = copy.deepcopy(hMC.Clone("bkgdens2"))
      eyls = array.array('f', range(0, hMC.GetNbinsX()))
      eyhs = array.array('f', range(0, hMC.GetNbinsX()))
      ratioups = copy.deepcopy(hMC.Clone("ratioups"))
      ratiodowns = copy.deepcopy(hMC.Clone("ratiodowns"))
      ymaxs = 2.                                                     
      
      for km in range(0, hMC.GetNbinsX()):
          contes1 =  hMC.GetBinError(km);      
          contes2 =  hMC.GetBinError(km);      
          dens1.SetBinContent (km, hMC.GetBinContent (km) + contes1);                                                                                                                           
          dens2.SetBinContent (km, hMC.GetBinContent (km) - contes2);                                                                                                                           
          ymaxs = hMC.GetBinContent(km) + contes1;                                                                                                                         
          eyls[km] = conte2;                                                                                                                                                           
          eyhs[km] = conte1;

      ratioups.Divide(dens1);                                           
      ratiodowns.Divide(dens2);                                             
      
      for km in range(0, ratiodata.GetNbinsX()):                                                       
          if (ratioups.GetBinContent (km) != 0):
              eyhs[km] = (1. / ratioups.GetBinContent (km) - 1)*ratiodata.GetBinContent (km);            
          else:                                                                                       
              eyhs[km] = 0;                                                                            
                                                                                                      
          if (ratiodowns.GetBinContent (km) != 0):
              eyls[km] = (1 - 1. / ratiodowns.GetBinContent (km))*ratiodata.GetBinContent (km);            
          else:                                                                                       
              eyls[km] = 0.                                                                            
     
      staterr = TGraphAsymmErrors(nvar, x, y, exl, exh, eyls, eyhs);
      
      errtottot.SetFillColor (r.kCyan+1);
      errtottot.SetFillStyle (3002);   
      errtot.SetFillColor (r.kBlue+2);
      errtot.SetFillStyle (3002);   
      err.SetFillColor (r.kGreen);
      err.SetFillStyle (3002);   
      staterr.SetFillColor (r.kRed+1);
      staterr.SetFillStyle (3002);   
      pad2.cd()

      line = TLine(ratio.GetBinLowEdge(1), 1, ratio.GetBinLowEdge(ratio.GetNbinsX()+1), 1)
      line.SetLineColor(r.kRed)
      ratio.Draw()
      if doAllErrors:
          legratio = TLegend(0.14,0.31,0.75,0.45);
          legratio.SetFillColor(0);
          legratio.SetBorderSize(0);
          legratio.SetNColumns(4);
          legratio.AddEntry(errtottot, "Uncl + JER + JES +  Stat","f");
          legratio.AddEntry(errtot, " JER + JES + Stat","f");
          legratio.AddEntry(err, "JES + Stat","f");
          legratio.AddEntry(staterr, "Stat","f");                                         
          errtottot.Draw("2 same")
          errtot.Draw("2 same")
          err.Draw("2 same")
          staterr.Draw("2 same")
          legratio.Draw("same")
          line.Draw("same")
          ratio.Draw("same");                                                      
      if puppi: 
          legratio = TLegend(0.14,0.31,0.4,0.45);                                   
          legratio.SetFillColor(0);
          legratio.SetBorderSize(0);
          legratio.SetNColumns(2);
          legratio.AddEntry(err, "JES + Stat","f");
          legratio.AddEntry(staterr, "Stat","f");                                  
          err.Draw("2")
          staterr.Draw("2 same")
          legratio.Draw("same")
          line.Draw("same")
          ratio.Draw("same");                                                      
      if statOnly:
          legratio = TLegend(0.14,0.32,0.43,0.5);
          legratio.SetFillColor(0);
          legratio.SetBorderSize(0);
          legratio.AddEntry(staterr, " ","f");                                  
          staterr.Draw("2 same")
         # legratio.Draw("same")
          line.Draw("same")
          ratio.Draw("same");                                                      
         

      pad1.cd()
      self.banner(isData, lumi, log, events, isNVert, isSig1jet, run)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)

   def save(self, legend, lumi, log,chisquare, value, title, option):

      events = 5  
      setUpAxis = 0 
      setLowAxis = 1 
      log = 0
      isNVert = 0
      puppi = 0
      isSig1jet = 0
      logx = 0
      if option == '':  
          events  = ''
          lumi  = 36.2
      if option == 'nvert':
          log = 0
          events = 0
          isNVert = 1
      if option == 'sig0':
          log =1
          events = 2
      if option == 'sig1':
          log =1
          events = 2
          isSig1jet = 1      
      if option == 'qt':
          setLowAxis = 1
          setUpAxis = 1
          log = 1           
      if option == 'rereco':      
          setLowAxis = 1
          setUpAxis = 1
          log = 1           
          logx = 1           
      self.myCanvas.cd()
      if(log):
          self.myCanvas.GetPad(0).SetLogy(1)
      if(logx):
          self.myCanvas.GetPad(0).SetLogx(1)

      for i in range(0, len(self.histos)):
          if setLowAxis:
              self.histos[i].SetMinimum(0.000001)
              self.histos[i].SetMaximum(1.)
              if isNVert:
                  self.histos[i].SetMaximum(10000000.)
          if(self.ToDraw[i] != 0):        
              self.histos[i].Draw(self.options[i])

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextSize(latex[-1])
          lat.SetTextFont(latex[-2])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      if(legend):
          self.makeLegend(log)
          self.myLegend.Draw()

      self.banner2(lumi, chisquare, value , title, events, isNVert, isSig1jet)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)


