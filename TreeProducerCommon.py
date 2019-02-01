import ROOT
import math 
import numpy as num 


var_dict = {
  'Electron_mvaFall17Iso':      'Electron_mvaFall17Iso',
  'Electron_mvaFall17Iso_WPL':  'Electron_mvaFall17Iso_WPL',
  'Electron_mvaFall17Iso_WP80': 'Electron_mvaFall17Iso_WP80',
  'Electron_mvaFall17Iso_WP90': 'Electron_mvaFall17Iso_WP90',
}

def setYear(year):
  """Help function to change the name of some variables that depend on the year."""
  if year==2018 or year==2016:
    print "setYear: setting var_dict to year %s"%(year)
    var_dict['Electron_mvaFall17Iso']      = 'Electron_mvaFall17V2Iso'
    var_dict['Electron_mvaFall17Iso_WPL']  = 'Electron_mvaFall17V2Iso_WPL'
    var_dict['Electron_mvaFall17Iso_WP80'] = 'Electron_mvaFall17V2Iso_WP80'
    var_dict['Electron_mvaFall17Iso_WP90'] = 'Electron_mvaFall17V2Iso_WP90'

def getvar(obj,var):
  """Help function to get some variable's real name from the dictionary."""
  return getattr(obj,var_dict[var])

def getVLooseTauIso(year):
  #if year==2016:
  #  return lambda e,i: ord(e.Tau_idMVAoldDM[i])>0
  #else:
  return lambda e,i: ord(e.Tau_idMVAoldDM[i])>0 or ord(e.Tau_idMVAnewDM2017v2[i])>0 or ord(e.Tau_idMVAoldDM2017v1[i])>0 or ord(e.Tau_idMVAoldDM2017v2[i])>0

def Tau_idIso(event,i):
  raw = event.Tau_rawIso[i]
  if event.Tau_photonsOutsideSignalCone[i]/event.Tau_pt[i]<0.10:
    return 0 if raw>4.5 else 1 if raw>3.5 else 3 if raw>2.5 else 7 if raw>1.5 else 15 if raw>0.8 else 31 # VVLoose, VLoose, Loose, Medium, Tight
  return 0 if raw>4.5 else 1 if raw>3.5 else 3 # VVLoose, VLoose
  
class TreeProducerCommon(object):
    
    def __init__(self, name):
        
        print 'TreeProducerCommon is called', name
        
        # create file
        self.outputfile = ROOT.TFile(name, 'RECREATE')
        self.tree = ROOT.TTree('tree','tree')
        
        # histogram for cutflow
        self.cutflow = ROOT.TH1F('cutflow', 'cutflow',  25, 0,  25)
        self.pileup  = ROOT.TH1F('pileup',  'pileup',  100, 0, 100)
        
        
        ###################
        # event variables #
        ###################
        
        self.add_branch("run",num.dtype(int))
        
        self.add_branch("run",num.dtype(int))
        self.add_branch("luminosityBlock",num.dtype(int))        
        self.add_branch("event",num.dtype(int))
        self.add_branch("MET_pt")
        self.add_branch("MET_phi")
        self.add_branch("GenMET_pt")
        self.add_branch("GenMET_phi")
        self.add_branch("PuppiMET_pt")
        self.add_branch("PuppiMET_phi")
        ###self.add_branch("MET_significance")
        ###self.add_branch("MET_covXX")
        ###self.add_branch("MET_covXY")
        ###self.add_branch("MET_covYY")
        ###self.add_branch("fixedGridRhoFastjetAll")
        self.add_branch("nPU",num.dtype(int))
        self.add_branch("nTrueInt",num.dtype(int))
        self.add_branch("npvs",num.dtype(int))
        self.add_branch("npvsGood",num.dtype(int))
        self.add_branch("LHE_Njets",num.dtype(int))
        self.add_branch("isData",num.dtype(bool))
        self.add_branch("genWeight")
        self.add_branch("weight")
        self.add_branch("trigweight")
        self.add_branch("puweight")
        self.add_branch("zptweight")
        self.add_branch("idisoweight_1")
        self.add_branch("idisoweight_2")
        self.add_branch("btagweight")
        self.add_branch("btagweight_deep")
        
        
        self.weight[0]          = 1.
        self.genweight[0]       = 1.
        self.trigweight[0]      = 1.
        self.puweight[0]        = 1.
        self.idisoweight_1[0]   = 1.
        self.idisoweight_2[0]   = 1.
        self.btagweight[0]      = 1.
        self.btagweight_deep[0] = 1.
        self.zptweight[0]       = 1.
        self.genmet[0]          = -1
        self.genmetphi[0]       = -9
        
        self.add_branch("njets",num.dtype(int))
        self.add_branch("njets50",num.dtype(int))
        self.add_branch("nfjets",num.dtype(int))
        self.add_branch("ncjets",num.dtype(int))
        self.add_branch("nbtag",num.dtype(int))
        self.add_branch("pfmt_1")
        self.add_branch("pfmt_2")
        
        self.add_branch("jpt_1")
        self.add_branch("jeta_1")
        self.add_branch("jphi_1")
        self.add_branch("jcsvv2_1")
        self.add_branch("jdeepb_1")
        
        self.add_branch("jpt_2")
        self.add_branch("jeta_2")
        self.add_branch("jphi_2")
        self.add_branch("jcsvv2_2")
        self.add_branch("jdeepb_2")
        
        self.add_branch("bpt_1")
        self.add_branch("beta_1")
        
        self.add_branch("bpt_2")
        self.add_branch("beta_2")
        
        self.add_branch("m_vis")
        self.add_branch("pt_tt")
        self.add_branch("dR_ll")
        self.add_branch("dphi_ll")
        
        self.add_branch("pzetamiss")
        self.add_branch("pzetavis")
        self.add_branch("pzeta_disc")
        
        self.add_branch("dilepton_veto",num.dtype(bool))
        self.add_branch("extraelec_veto",num.dtype(bool))
        self.add_branch("extramuon_veto",num.dtype(bool))
        
        self.add_branch("ngentauhads",num.dtype(int))
        self.add_branch("ngentaus",num.dtype(int))
        self.add_branch("m_genboson")
        self.add_branch("pt_genboson")
        
        
        self.nPU[0]         = -1
        self.nTrueInt[0]    = -1
        self.LHE_Njets[0]   = -1
        self.m_genboson[0]  = -1
        self.pt_genboson[0] = -1
    
    def add_branch(self, name, dtype=num.dtype(float)):
        #name should be a string
        #dtype should be a numpy.dtype
        
        #used to translate numpy to root nomenclature
        dtype_translation={
            '?': "O",
            'b': "B",
            'B': "b",
            'i': "I",
            'u': "i",
            'f': "D"
            # 'c':
            # 'm':
            # 'M':
            # 'O':
            # 'S':
            # 'a':
            # 'U':
            # 'V':
        }
        
        #numpy:
        # '?' 	boolean
        # 'b' 	(signed) byte
        # 'B' 	unsigned byte
        # 'i' 	(signed) integer
        # 'u' 	unsigned integer
        # 'f' 	floating-point
        # 'c' 	complex-floating point
        # 'm' 	timedelta
        # 'M' 	datetime
        # 'O' 	(Python) objects
        # 'S', 'a' 	zero-terminated bytes (not recommended)
        # 'U' 	Unicode string
        # 'V' 	raw data (void)
        
        #root:
        # C : a character string terminated by the 0 character
        # B : an 8 bit signed integer (Char_t)
        # b : an 8 bit unsigned integer (UChar_t)
        # S : a 16 bit signed integer (Short_t)
        # s : a 16 bit unsigned integer (UShort_t)
        # I : a 32 bit signed integer (Int_t)
        # i : a 32 bit unsigned integer (UInt_t)
        # F : a 32 bit floating point (Float_t)
        # D : a 64 bit floating point (Double_t)
        # L : a 64 bit signed integer (Long64_t)
        # l : a 64 bit unsigned integer (ULong64_t)
        # O : [the letter o, not a zero] a boolean (Bool_t)

        setattr(self,name,num.zeros(1, dtype=dtype))
        self.tree.Branch(name, getattr(self,name) , "{0}/{1}".format(name,dtype_translation[dtype.kind]))
        


class DiLeptonBasicClass:
    def __init__(self, id1, pt1, iso1, id2, pt2, iso2):
        self.id1  = id1
        self.id2  = id2
        self.pt1  = pt1
        self.pt2  = pt2
        self.iso1 = iso1
        self.iso2 = iso2
        
    def __gt__(self, odilep):
        """Order dilepton pairs according to the pT of both objects first, then in isolation."""
        if   self.pt1  != odilep.pt1:  return self.pt1  > odilep.pt1  # greater = higher pT
        elif self.pt2  != odilep.pt2:  return self.pt2  > odilep.pt2  # greater = higher pT
        elif self.iso1 != odilep.iso1: return self.iso1 < odilep.iso1 # greater = smaller isolation
        elif self.iso2 != odilep.iso2: return self.iso2 < odilep.iso2 # greater = smaller isolation
        return True
    
class LeptonTauPair(DiLeptonBasicClass):
    def __gt__(self, oltau):
        """Override for tau isolation."""
        if   self.pt1  != oltau.pt1:  return self.pt1  > oltau.pt1  # greater = higher pT
        elif self.pt2  != oltau.pt2:  return self.pt2  > oltau.pt2  # greater = higher pT
        elif self.iso1 != oltau.iso1: return self.iso1 < oltau.iso1 # greater = smaller lepton isolation
        elif self.iso2 != oltau.iso2: return self.iso2 > oltau.iso2 # greater = larger tau isolation
        return True
    
class DiTauPair(DiLeptonBasicClass):
    def __gt__(self, oditau):
        """Override for tau isolation."""
        if   self.pt1  != oditau.pt1:  return self.pt1  > oditau.pt1  # greater = higher pT
        elif self.pt2  != oditau.pt2:  return self.pt2  > oditau.pt2  # greater = higher pT
        elif self.iso1 != oditau.iso1: return self.iso1 > oditau.iso1 # greater = larger tau isolation
        elif self.iso2 != oditau.iso2: return self.iso2 > oditau.iso2 # greater = larger tau isolation
        return True
    


def bestDiLepton(diLeptons):
    """Take best dilepton pair."""
    if len(diLeptons)==1:
        return diLeptons[0]
    #least_iso_highest_pt = lambda dl: (-dl.tau1_pt, -dl.tau2_pt, dl.tau2_iso, -dl.tau1_iso)
    #return sorted(diLeptons, key=lambda dl: least_iso_highest_pt(dl), reverse=False)[0]
    return sorted(diLeptons, reverse=True)[0]
    

def deltaR2( e1, p1, e2, p2):
    de = e1 - e2
    dp = deltaPhi(p1, p2)
    return de*de + dp*dp


def deltaR( *args ):
    return math.sqrt( deltaR2(*args) )


def deltaPhi( p1, p2):
    """Computes delta phi, handling periodic limit conditions."""
    res = p1 - p2
    while res > math.pi:
      res -= 2*math.pi
    while res < -math.pi:
      res += 2*math.pi
    return res
    

def extraLeptonVetos(event, muon_idxs, electron_idxs, name):
    
    b_dilepton_veto_  = False;
    b_extraelec_veto_ = False;
    b_extramuon_veto_ = False;
    
    LooseMuons = [ ]
    for imuon in range(event.nMuon):
        if event.Muon_pt[imuon] < 10: continue
        if abs(event.Muon_eta[imuon]) > 2.4: continue
        if abs(event.Muon_dz[imuon]) > 0.2: continue
        if abs(event.Muon_dxy[imuon]) > 0.045: continue
        if event.Muon_pfRelIso04_all[imuon] > 0.3: continue
        if event.Muon_mediumId[imuon] > 0.5 and (imuon not in muon_idxs):
            b_extramuon_veto_ = True
        if event.Muon_pt[imuon] > 15 and event.Muon_isPFcand[imuon] > 0.5:
            LooseMuons.append(imuon)
    
    LooseElectrons = [ ]
    for ielectron in range(event.nElectron):
        if event.Electron_pt[ielectron] < 10: continue
        if abs(event.Electron_eta[ielectron]) > 2.5: continue
        if abs(event.Electron_dz[ielectron]) > 0.2: continue
        if abs(event.Electron_dxy[ielectron]) > 0.045: continue
        if event.Electron_pfRelIso03_all[ielectron] > 0.3: continue
        #if event.Electron_convVeto[ielectron] ==1 and ord(event.Electron_lostHits[ielectron]) <= 1 and event.Electron_mvaFall17Iso_WP90[ielectron] > 0.5 and (ielectron not in electron_idxs):
        if event.Electron_convVeto[ielectron] ==1 and ord(event.Electron_lostHits[ielectron]) <= 1 and getvar(event,'Electron_mvaFall17Iso_WP90')[ielectron] > 0.5 and (ielectron not in electron_idxs):
            b_extraelec_veto_ = True
        if event.Electron_pt[ielectron] > 15 and getvar(event,'Electron_mvaFall17Iso_WPL')[ielectron] > 0.5:
            LooseElectrons.append(ielectron)
    
    if name.find('mutau')!=-1:
      for idx1 in LooseMuons:
        for idx2 in LooseMuons:
            if idx1 >= idx2: continue 
            dR = deltaR(event.Muon_eta[idx1], event.Muon_phi[idx1], 
                        event.Muon_eta[idx2], event.Muon_phi[idx2])
            if event.Muon_charge[idx1] * event.Muon_charge[idx2] < 0 and dR > 0.15:
                b_dilepton_veto_ = True
    
    if name.find('eletau')!=-1:
      for idx1 in LooseElectrons:
        for idx2 in LooseElectrons:
            if idx1 >= idx2: continue 
            dR = deltaR(event.Electron_eta[idx1], event.Electron_phi[idx1], 
                        event.Electron_eta[idx2], event.Electron_phi[idx2])
            if event.Electron_charge[idx1] * event.Electron_charge[idx2] < 0 and dR > 0.15:
                b_dilepton_veto_ = True
    
    return b_extramuon_veto_, b_extraelec_veto_, b_dilepton_veto_

#  b_dilepton_veto[ch]                   = (int) b_dilepton_veto_;
#  b_extraelec_veto[ch]                  = (int) b_extraelec_veto_;
#  b_extramuon_veto[ch]                  = (int) b_extramuon_veto_;
#  b_lepton_vetos[ch]                    = ( b_dilepton_veto_ || b_extraelec_veto_ || b_extramuon_veto_ );


