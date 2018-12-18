import ROOT
import math

#/shome/ineuteli/analysis/LQ_2017/NanoTreeProducer/leptonSF
#/shome/ytakahas/work/Leptoquark/CMSSW_9_4_4/src/PhysicsTools/NanoAODTools/NanoTreeProducer/leptonSF
path = 'pileup/'

class PileupWeightTool:
    
    def __init__( self ):
        """Load data and MC files."""
        datfile       = ROOT.TFile( path+'Data_PileUp_2017_69p2.root', 'r')
        mcfile        = ROOT.TFile( path+'MC_PileUp_Winter17_PU25ns_V2_fromMC.root', 'r')
        self.datahist = datfile.Get('pileup')
        self.mchist   = mcfile.Get('pileup')
        self.datahist.SetDirectory(0)
        self.mchist.SetDirectory(0)
        datafile.Close()
        mcfile.Close()
    
    
    def getLeptonTauFakeSF(npu):
        """Get pileup weight for a given number of pileup interactions."""
        
        data = self.datahist.GetBinContent(datahist.FindBin(npu))
        mc   = self.mchist.GetBinContent(mchist.FindBin(npu))
        
        if mc>0.:
          return data/mc
        
        print "PileupWeightTools: could not make pileup weight: data=%s, mc=%s"%(data,mc)  
        return 1

