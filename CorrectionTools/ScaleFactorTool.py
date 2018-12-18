#! /bin/usr/env bash
# Author: Izaak Neutelings (November 2018)
import os, re
from ROOT import TFile #, TH2F, TGraphAsymmErrors, Double()


def ensureTFile(filename,option='READ'):
  """Open TFile, checking if the file in the given path exists."""
  if not os.path.isfile(filename):
    print '>>> ERROR! ScaleFactorTool::ensureTFile: File in path "%s" does not exist!!'%(filename)
    exit(1)
  file = TFile(filename,option)
  if not file or file.IsZombie():
    print '>>> ERROR! ScaleFactorTool::ensureTFile Could not open file by name "%s"'%(filename)
    exit(1)
  return file
  
#def getSF_ptvseta(sftool, pt, eta):
#    """Get SF for a given pT vs. eta."""
#    #abseta = abs(eta)
#    xbin   = sftool.hist.GetXaxis().FindBin(eta)
#    ybin   = sftool.hist.GetYaxis().FindBin(pt)
#    sf     = sftool.hist.GetBinContent(xbin, ybin)
#    #print "ScaleFactor::getSF: %s, pt = %6.2f, eta = %6.3f, data = %6.3f, mc = %6.3f, sf = %6.3f"%(self.name,pt,eta,data,mc,sf)
#    #print "ScaleFactor::getSF_ptvseta: %s, pt = %6.2f, eta = %6.3f, sf = %6.3f"%(sftool.name,pt,eta,sf)
#    return sf
#    
#def getSF_etavspt(sftool, pt, eta):
#    """Get SF for a given pT vs. eta."""
#    #abseta = abs(eta)
#    xbin   = sftool.hist.GetXaxis().FindBin(pt)
#    ybin   = sftool.hist.GetYaxis().FindBin(eta)
#    sf     = sftool.hist.GetBinContent(xbin, ybin)
#    #print "ScaleFactor::getSF: %s, pt = %6.2f, eta = %6.3f, data = %6.3f, mc = %6.3f, sf = %6.3f"%(self.name,pt,eta,data,mc,sf)
#    #print "ScaleFactor::getSF_ptvseta: %s, pt = %6.2f, eta = %6.3f, sf = %6.3f"%(sftool.name,pt,eta,sf)
#    return sf
    


class ScaleFactor:
    
    def __init__(self, filename, histname, name="<noname>", ptvseta=True):
        #print '>>> ScaleFactor::init("%s","%s",name="%s",ptvseta=%r)'%(filename,histname,name,ptvseta)
        self.name     = name
        self.filename = filename
        self.file     = ensureTFile(filename)
        self.hist     = self.file.Get(histname)
        self.hist.SetDirectory(0)
        self.file.Close()
        
        #xtitle = self.hist.GetXaxis().GetTitle().lower().replace('_','').replace('{','').replace('}','')
        #ytitle = self.hist.GetYaxis().GetTitle().lower().replace('_','').replace('{','').replace('}','')
        #if not ptvseta and 'pt' in ytitle and 'eta' in xtitle:
        #  print '>>> Warning! ScaleFactor::init: ptvseta=False, but xtitle="%s" and ytitle="%s" for "%s" in "%s"'%(
        #  self.hist.GetXaxis().GetTitle(),self.hist.GetYaxis().GetTitle(),self.hist.GetName(),self.filename)
        #elif   ptvseta and 'pt' in xtitle and 'eta' in ytitle:
        #  print '>>> Warning! ScaleFactor::init: ptvseta=True, but xtitle="%s" and ytitle="%s" for "%s" in "%s"'%(
        #  self.hist.GetXaxis().GetTitle(),self.hist.GetYaxis().GetTitle(),self.hist.GetName(),self.filename)
        #if ptvseta:
        #  self.getSFMethod = getSF_ptvseta
        #else:
        #  self.getSFMethod = getSF_etavspt
        
    #def getSF(self, pt, eta):
    #    """Get SF for a given pT, eta."""
    #    return self.getSFMethod(self,pt,eta)
        
    def getSF(self, pt, eta):
        """Get SF for a given pT, eta."""
        #abseta = abs(eta)
        xbin   = self.hist.GetXaxis().FindBin(eta)
        ybin   = self.hist.GetYaxis().FindBin(pt)
        sf     = self.hist.GetBinContent(xbin,ybin)
        #print "ScaleFactor::getSF: %s, pt = %6.2f, eta = %6.3f, data = %6.3f, mc = %6.3f, sf = %6.3f"%(self.name,pt,eta,data,mc,sf)
        #print "ScaleFactor::getSF: %s, pt = %6.2f, eta = %6.3f, sf = %6.3f"%(self.name,pt,eta,sf)
        return sf
    


class ScaleFactorHTT:
    
    def __init__(self, filename, graphname="ZMass", name="<noname>"):
        #print '>>> ScaleFactor::init("%s","%s",name="%s")'%(filename,graphname,name)
        self.name      = name
        self.filename  = filename
        self.file      = ensureTFile(filename)
        self.hist_eta  = self.file.Get('etaBinsH')
        self.hist_eta.SetDirectory(0)
        self.effs_data = { }
        self.effs_mc   = { }
        for ieta in range(1,self.hist_eta.GetXaxis().GetNbins()+1):
          etalabel = self.hist_eta.GetXaxis().GetBinLabel(ieta)
          self.effs_data[etalabel] = self.file.Get(graphname+etalabel+"_Data")
          self.effs_mc[etalabel]   = self.file.Get(graphname+etalabel+"_MC")
        self.file.Close()
        
    def getSF(self, pt, eta):
        """Get SF for a given pT, eta."""
        abseta = abs(eta)
        etabin = self.hist_eta.GetXaxis().GetBinLabel(min(self.hist_eta.GetXaxis().GetNbins(),self.hist_eta.GetXaxis().FindBin(abseta)))
        data   = self.effs_data[etabin].Eval(pt)
        mc     = self.effs_mc[etabin].Eval(pt)
        if mc==0:
          sf   = 1.0
        else:
          sf   = data/mc
        #print "ScaleFactorHTT::getSF: %s, pt = %6.2f, eta = %6.3f, data = %6.3f, mc = %6.3f, sf = %6.3f"%(self.name,pt,eta,data,mc,sf)
        return sf
    


##etaLt = re.compile(r"EtaLt(\dp\d+)")
##etaTo = re.compile(r"Eta(\dp\d+to\dp\d+)")
##etaGt = re.compile(r"EtaGt(\dp\d+)")
#def getEtaRangeFromString(eta):
#    """Get eta range from string."""
#    etainf = 7.0
#    match  = re.match(r"EtaLt(\dp\d+)",eta)
#    if match:
#      return (0,float(match.group(1).replace('p','.')))
#    match = re.match(r"Eta(\dp\d+)to(\dp\d+)",eta)
#    if match:
#      return (float(match.group(1).replace('p','.')),float(match.group(2).replace('p','.')))
#    match = re.match(r"EtaGt(\dp\d+)",eta)
#    if match:
#      return (float(match.group(1).replace('p','.')),etainf)
#    print "ERROR! getEtaRange: Could not find a eta range pattern for the string '%s'"%(eta)
#    return None

#def getBinsFromTGraph(graph):
#    """Get xbins from TGraph."""
#    x, y  = Double(), Double()
#    xlast  = None
#    xbins = [ ]
#    for i in range(0,graph.GetN()):
#      graph.GetPoint(i,x,y)
#      xlow = float(x) - graph.GetErrorXlow(i)
#      xup  = float(x) + graph.GetErrorXhigh(i)
#      if xlow>=xup:
#        print 'Warning! getBinsFromTGraph: Point i=%d of graph "%s": lower x value %.1f >= upper x value %.1f.'%(i,graph.GetName(),xlow,xup)
#      if xlast!=None and abs(xlast-xlow)>1e-5:
#        print 'Warning! getBinsFromTGraph: Point i=%d of graph "%s": lower x value %.1f does not conincide with upper x value of last point, %.1f.'%(i,graph.GetName(),xlow,xlast)
#      xbins.append(xlow)
#      xlast = xup
#    xbins.append(xlast)
#    return xbins
