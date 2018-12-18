from ROOT import TFile
from ScaleFactorTool import ScaleFactor, ScaleFactorHTT

#/shome/ineuteli/analysis/LQ_2017/NanoTreeProducer/leptonSF
#/shome/ytakahas/work/Leptoquark/CMSSW_9_4_4/src/PhysicsTools/NanoAODTools/NanoTreeProducer/leptonSF
path    = 'CorrectionTools/leptonEfficiencies/'
pathHTT = 'CorrectionTools/leptonEfficiencies/HTT/Electron/Run2017/'


class ElectronSFs:
    
    def __init__( self ):
        """Load histograms from files."""
        
        # TRIGGER (HTT)
        self.sftool_trig = ScaleFactorHTT(pathHTT+"Electron_Ele32orEle35.root","ZMass",'ele_trig')
        
        # RECO, IDISO (EGamme POG)
        self.sftool_reco  = ScaleFactor(path+"egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","EGamma_SF2D",'ele_reco',ptvseta=True)
        self.sftool_idiso = ScaleFactor(path+"gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp80iso.root","EGamma_SF2D",'ele_idiso',ptvseta=True)
        #self.sftool_idiso = ScaleFactorHTT(pathHTT+"Electron_IdIso_IsoLt0.15_IsoID_eff.root","ZMass",'ele_idiso')
        
    def getTriggerSF( self, pt, eta ):
        """Get SF for single electron trigger."""
        return self.sftool_trig.getSF(pt,eta)
        
    def getIdIsoSF( self, pt, eta ):
        """Get SF for electron identification + isolation."""
        sf_reco  = self.sftool_reco.getSF(pt,eta)
        sf_idiso = self.sftool_idiso.getSF(pt,eta)
        return sf_reco*sf_idiso
    
    def getLeptonTauFakeSF(genmatch,eta):
        """Get SF for lepton to tau fake."""
        # https://indico.cern.ch/event/715039/timetable/#2-lepton-tau-fake-rates-update
        # https://indico.cern.ch/event/719250/contributions/2971854/attachments/1635435/2609013/tauid_recommendations2017.pdf
        # https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Muon%20to%20tau%20fake%20rate
        eta = abs(eta)
        
        # electron -> tau (Tight for mutau)
        if genmatch==1:
          if   eta<1.460: return 1.80
          elif eta>1.558: return 1.53
        
        # muon -> tau (Loose for mutau)
        elif genmatch==2:
          if   eta<0.4: return 1.061
          elif eta<0.8: return 1.022
          elif eta<1.2: return 1.097
          elif eta<1.7: return 1.030
          else:         return 1.941
        
        # real tau (Tight)
        #elif genmatch_2==5
        #  return 0.88; // Tight
        
        return 1.0
        

