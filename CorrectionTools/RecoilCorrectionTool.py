#! /bin/usr/env python
# Author: Izaak Neutelings (November 2018)
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ScaleFactorTool import ensureTFile
from ROOT import TLorentzVector

path = 'CorrectionTools/Zpt/'

class RecoilCorrectionTool:
    
    def __init__( self, year=2017 ):
        """Load Z pT weights."""
        #if year==2016:
        #  self.file = ensureTFile( path+'Zpt_weights_2017_Izaak.root', 'READ')
        #elif year==2017:
        #  self.file = ensureTFile( path+'Zpt_weights_2017_Izaak.root', 'READ')
        #else:
        #  self.file = ensureTFile( path+'Zpt_weights_2017_Izaak.root', 'READ')
        #self.file = self.file.Get('zptmass_weights')
        #self.hist.SetDirectory(0)
        #self.file.Close()
        
    
    def getZptWeight(self,Zpt,Zmass):
        """Get Z pT weight for a given Z boson pT and mass."""
        return 1.
        #weight = self.hist.GetBinContent(self.hist.GetXaxis().FindBin(Zpt),self.hist.GetXaxis().FindBin(Zmass))
        #print ">>> Warning! RecoilCorrectionTool::getZptWeight: Could not make pileup weight for npu=%s data=%s, mc=%s"%(npu,data,mc)
        #return weight
        

def getZPTMass(event):
    """Calculate Z boson pT and mass."""
    #print '-'*80
    genparticles = Collection(event,'GenPart')
    zboson = TLorentzVector()
    for id in range(event.nGenPart):
      particle = genparticles[id]
      PID      = abs(particle.pdgId)
      #if PID==23 and particle.status==62:
      #  print "%3d: PID=%3d, mass=%3.1f, pt=%3.1f, status=%2d"%(id,particle.pdgId,particle.mass,particle.pt,particle.status)
      if ((PID==11 or PID==13) and particle.status==1 and hasBit(particle.statusFlags,9)) or\
                     (PID==15  and particle.status==2 and hasBit(particle.statusFlags,9)):
        zboson += particle.p4()
        #print "%3d: PID=%3d, mass=%3.1f, pt=%3.1f, status=%2d, statusFlags=%2d (%16s), fromHardProcess=%2d"%(id,particle.pdgId,particle.mass,particle.pt,particle.status,particle.statusFlags,bin(particle.statusFlags),hasBit(particle.statusFlags,9))
    #print "tlv: mass=%3.1f, pt=%3.1f"%(zboson.M(),zboson.Pt())
    return zboson
    
def hasBit(value,bit):
  #return bin(value)[-bit]=='1'
  #return format(value,'b').zfill(bit)[-bit]=='1'
  return (value & (1 << (bit-1)))>0
