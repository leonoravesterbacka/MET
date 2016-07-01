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
      self.myLegend.SetTextSize(0.04)
      self.myLegend.SetLineWidth(0)
      self.myLegend.SetBorderSize(0)


   def banner(self, isData, lumi):
    
      latex = TLatex()                                
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.05);
      latex.DrawLatex(0.21, 0.93, "CMS")

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.03);
 
      if(isData):
        latexb.DrawLatex(0.33, 0.93, "Preliminary")
      else:
        latexb.DrawLatex(0.33, 0.93, "Simulation")

      text_lumi = "0.8 fb^{-1} (13 TeV, 2016)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      latexc.DrawLatex(0.90, 0.93, text_lumi)          

      latexd = TLatex()
      latexd.SetNDC();
      latexd.SetTextAngle(90);
      latexd.SetTextColor(r.kBlack);
      latexd.SetTextFont(42);
      latexd.SetTextAlign(31);
      latexd.SetTextSize(0.04);
      latexd.DrawLatex(0.035, 0.93, "Events / 5 GeV")          



   def banner2(self, isData, chisquare):
    
      latex = TLatex()
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(31);
      latex.SetTextSize(0.05);
      latex.DrawLatex(0.23, 0.93, "CMS")

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.03);
 
      if(isData):
        latexb.DrawLatex(0.37, 0.93, "Preliminary")
      else:
        latexb.DrawLatex(0.37, 0.93, "Simulation")

      text_lumi = "0.8 fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      latexc.DrawLatex(0.90, 0.93, text_lumi)
      latexc.DrawLatex(0.85, 0.8,  "#chi^{2} = %.2f " %(chisquare))

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

   def addLatex(self, x1, y1, text, font=42, size = 0.04):
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
       

 
   def makeLegend(self):

      for i in range(0, len(self.histos)):
          for j in range(0, len(self.orderForLegend)):
              if(self.orderForLegend[j] != -1 and self.orderForLegend[j] == i):
                  self.myLegend.AddEntry(self.histos[j], self.labels[j], self.labelsOption[j])
          

   def ensurePath(self, _path):
      d = os.path.dirname(_path)
      if not os.path.exists(d):
         os.makedirs(d)

   def saveRatio(self, legend, isData, log, lumi, hdata, hMC, hjecUp, hjecDown , r_ymin=0, r_ymax=2):
      self.myCanvas.cd()
      pad1 = TPad("pad1", "pad1", 0, 0.2, 1, 1.0) 
      pad1.SetBottomMargin(0.12)
      pad1.Draw()
      pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.2)
      pad2.SetTopMargin(0.1);
      pad2.SetBottomMargin(0.3);
      pad2.Draw();

      pad1.cd()
      if(log):
          pad1.SetLogy(1)

      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):
              self.histos[i].Draw(self.options[i])

      if(legend):
          self.makeLegend()
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

      ratio.SetTitle("")
      if hjecUp.Integral() == hjecDown.Integral():
          r_ymin = 0.5
          r_ymax = 1.5
      ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
      ratio.GetYaxis().SetTitle("Data / MC")
      ratio.GetYaxis().CenterTitle();
      ratio.GetYaxis().SetLabelSize(0.22);
      ratio.GetXaxis().SetLabelSize(0.22);
      ratio.GetYaxis().SetTitleOffset(0.3);
      ratio.GetYaxis().SetNdivisions(4);
      ratio.GetYaxis().SetTitleSize(0.22);
      ratio.GetXaxis().SetTitleSize(0.22);
      ratio.SetMarkerSize(0.6*ratio.GetMarkerSize());
      ratio.GetXaxis().SetTitle('');
                                                                                                                                                                                           
      #make some jec errors      
      den1 = copy.deepcopy(hMC.Clone("bkgden1"))
      den2 = copy.deepcopy(hMC.Clone("bkgden2"))
      	                                                                                                                                                                                 
      nvar = hMC.GetNbinsX()                                                                                                                                                           
      x = array.array('f', range(0, hMC.GetNbinsX()))
      y = array.array('f', range(0, hMC.GetNbinsX()))
      exl = array.array('f', range(0, hMC.GetNbinsX()))
      eyl = array.array('f', range(0, hMC.GetNbinsX()))
      exh = array.array('f', range(0, hMC.GetNbinsX()))
      eyh = array.array('f', range(0, hMC.GetNbinsX()))
      ratioup = copy.deepcopy(hMC.Clone("ratioup"))
      ratiodown = copy.deepcopy(hMC.Clone("ratiodown"))
      ymax = 2.                                                                                                                                                              
      
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
          exl[km] = hMC.GetBinWidth (km) / 2;                                                                                                                                 
          exh[km] = hMC.GetBinWidth (km) / 2;                                                                                                                                 
          eyl[km] = conte2;                                                                                                                                                           
          eyh[km] = conte1;                                                                                                                                                                        
      
      ratioup.Divide(den1);                                                                                                                                                         
      ratiodown.Divide(den2);                   
      ratiodata = copy.deepcopy(hdata.Clone("ratiodata"))
      ratiodata.Divide (hMC);                                                                                                                                                    
                                                                                                                                                                                         
      for km in range(0, ratiodata.GetNbinsX()):
          if (ratiodata.GetBinContent (km) > ymax):
              ymax = ratiodata.GetBinContent (km) + ratiodata.GetBinError (km);                                                                                                              
          x[km] = ratiodata.GetBinCenter (km);                                                                                                                                          
          y[km] = 1;	                                                                                                                                                             
          exl[km] = ratiodata.GetBinWidth (km) / 2;                                                                                                                                     
          exh[km] = ratiodata.GetBinWidth (km) / 2;                                                                                                                                     
                                                                                                                                                                                         
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
      ratiodatas = copy.deepcopy(hdata.Clone("ratiodatas"))              
      ratiodatas.Divide (hMC);                                          
      
      for km in range(0, ratiodatas.GetNbinsX()):                                                       
          if (ratioups.GetBinContent (km) != 0):
              eyhs[km] = (1. / ratioups.GetBinContent (km) - 1)*ratiodatas.GetBinContent (km);            
          else:                                                                                       
              eyhs[km] = 0;                                                                            
                                                                                                      
          if (ratiodowns.GetBinContent (km) != 0):
              eyls[km] = (1 - 1. / ratiodowns.GetBinContent (km))*ratiodatas.GetBinContent (km);            
          else:                                                                                       
              eyls[km] = 0.                                                                            
     
      staterr = TGraphAsymmErrors(nvar, x, y, exl, exh, eyls, eyhs);
      
      err.SetFillColor (r.kRed);
      err.SetFillStyle (3002);   
      staterr.SetFillColor (r.kGreen);
      staterr.SetFillStyle (3002);   

      pad2.cd();  
      line = TLine(ratio.GetBinLowEdge(1), 1, ratio.GetBinLowEdge(ratio.GetNbinsX()+1), 1)
      line.SetLineColor(r.kRed)
      ratio.Draw()
      if hjecUp.Integral() == hjecDown.Integral(): # this case, where the integrals of the up and down are the same, should only have stat errors 
          print "doing only stat errors"
          staterr.Draw("2 same");
      else:
          print "doing JEC + stat errors"
          err.Draw("2 same")
          staterr.Draw("2 same")
      line.Draw()
      ratio.Draw("same");

      pad1.cd()
      self.banner(isData, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)

   def save(self, legend, isData, log,chisquare):

      self.myCanvas.cd()
      if(log):
          self.myCanvas.GetPad(0).SetLogy(1)
     
      for i in range(0, len(self.histos)):
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
          self.makeLegend()
          self.myLegend.Draw()

      self.banner2(isData, chisquare )
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)



