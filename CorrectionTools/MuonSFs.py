from ScaleFactorTool import ScaleFactor, ScaleFactorHTT

# /shome/ytakahas/work/Leptoquark/CMSSW_9_4_4/src/PhysicsTools/NanoAODTools/NanoTreeProducer/leptonSF
# HTT: https://github.com/CMS-HTT/LeptonEfficiencies
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
# https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#Efficiency_Scale_Factors
path    = 'CorrectionTools/leptonEfficiencies/'
pathHTT = 'CorrectionTools/leptonEfficiencies/HTT/Muon/Run2017/'

class MuonSFs:
    
    def __init__(self, year=2017):
        # Load the TH1s containing the bin by bin values
        
        # TRIGGER (Muon POG)
        self.sftool_trig = ScaleFactor(path+"EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","IsoMu27_PtEtaBins/abseta_pt_ratio",'mu_trig',ptvseta=True)
        
        ## TRIGGER (HTT)
        #self.sftool_trig = ScaleFactorHTT(pathHTT+"Muon_IsoMu24orIsoMu27.root","ZMass",'mu_idiso')
        
        # ID ISO (HTT)
        self.sftool_idiso = ScaleFactorHTT(pathHTT+"Muon_IdIso_IsoLt0p15_eff_RerecoFall17.root","ZMass",'mu_idiso')
        
        ## ID (Muon POG)
        #self.sftool_id  = ScaleFactor(path+"RunBCDEF_SF_ID.root","NUM_MediumID_DEN_genTracks",'mu_id')
        #self.sftool_iso = ScaleFactor(path+"RunBCDEF_SF_ISO.root","NUM_TightRelIso_DEN_MediumID",'mu_iso')
        
    def getTriggerSF(self, pt, eta):
        """Get SF for single muon trigger."""
        return self.sftool_trig.getSF(pt,abs(eta))
        
    def getIdIsoSF(self, pt, eta):
        """Get SF for muon identification + isolation."""
        return self.sftool_idiso.getSF(pt,eta)
        
    #def getLeptonTauFakeSF(self, genmatch, eta):
    #    """Get SF for lepton to tau fake."""
    #    # https://indico.cern.ch/event/715039/timetable/#2-lepton-tau-fake-rates-update
    #    # https://indico.cern.ch/event/719250/contributions/2971854/attachments/1635435/2609013/tauid_recommendations2017.pdf
    #    # https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Muon%20to%20tau%20fake%20rate
    #    eta = abs(eta)
    #    
    #    # electron -> tau (VLoose for etau)
    #    if genmatch==1:
    #      if   eta<1.460: return 1.09
    #      elif eta>1.558: return 1.19
    #    
    #    # muon -> tau (Tight for etau)
    #    elif genmatch==2:
    #      if   eta<0.4: return 1.165
    #      elif eta<0.8: return 1.290
    #      elif eta<1.2: return 1.137
    #      elif eta<1.7: return 0.927
    #      else:         return 1.607
    #    
    #    # real tau (Tight)
    #    #elif genmatch_2==5
    #    #  return 0.88; // Tight
    #    
    #    return 1.0


