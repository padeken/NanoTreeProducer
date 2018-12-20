#! /bin/usr/env python
# Author: Izaak Neutelings (November 2018)
from ROOT import TFile

path = 'CorrectionTools/pileup/'

class PileupWeightTool:
    
    def __init__( self ):
        """Load data and MC pilup profiles."""
        self.datafile = TFile( path+'Data_PileUp_2017_69p2.root', 'READ')
        self.mcfile   = TFile( path+'MC_PileUp_Winter17_PU25ns_V2_fromMC.root', 'READ')
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
        
        data = self.datahist.GetBinContent(self.datahist.FindBin(npu))
        mc   = self.mchist.GetBinContent(self.mchist.FindBin(npu))
        
        if mc>0.:
          return data/mc
        
        print ">>> PileupWeightTools::getWeight: could not make pileup weight: data=%s, mc=%s"%(data,mc)  
        return 1.
    
