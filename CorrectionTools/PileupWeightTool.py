#! /bin/usr/env python
# Author: Izaak Neutelings (November 2018)
from ROOT import TFile
from ScaleFactorTool import ensureTFile

path = 'CorrectionTools/pileup/'

class PileupWeightTool:
    
    def __init__( self, year=2017 ):
        """Load data and MC pilup profiles."""
        if year==2017:
          self.datafile = ensureTFile( path+'Data_PileUp_2017_69p2.root', 'READ')
          self.mcfile   = ensureTFile( path+'MC_PileUp_Winter17_PU25ns_V2_fromMC.root', 'READ')
        else:
          self.datafile = ensureTFile( path+'Data_PileUp_2018_69p2.root', 'READ')
          self.mcfile   = ensureTFile( path+'MC_PileUp_2018_Autumn18.root', 'READ')
        self.datahist = self.datafile.Get('pileup')
        self.mchist   = self.mcfile.Get('pileup')
        self.datahist.SetDirectory(0)
        self.mchist.SetDirectory(0)
        self.datahist.Scale(1./self.datahist.Integral())
        self.mchist.Scale(1./self.mchist.Integral())
        self.datafile.Close()
        self.mcfile.Close()
        
    
    def getWeight(self,npu):
        """Get pileup weight for a given number of pileup interactions."""
        data = self.datahist.GetBinContent(self.datahist.GetXaxis().FindBin(npu))
        mc   = self.mchist.GetBinContent(self.mchist.GetXaxis().FindBin(npu))
        if mc>0.:
          return data/mc
        print ">>> Warning! PileupWeightTools::getWeight: Could not make pileup weight for npu=%s data=%s, mc=%s"%(npu,data,mc)  
        return 1.
    
