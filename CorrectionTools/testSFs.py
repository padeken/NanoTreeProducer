#! /usr/bin/env python
# Author: Izaak Neutelings (December 2018)

import time
start0 = time.time()
print
print ">>> importing modules..."

start1 = time.time()
from ROOT import TFile
print ">>>   imported ROOT classes after %.1f seconds"%(time.time()-start1)

start1 = time.time()
from ScaleFactorTool import *
print ">>>   imported ScaleFactorTool classes after %.1f seconds"%(time.time()-start1)

start1 = time.time()
from MuonSFs import *
print ">>>   imported MuonSFs classes after %.1f seconds"%(time.time()-start1)

start1 = time.time()
from ElectronSFs import *
print ">>>   imported ElectronSFs classes after %.1f seconds"%(time.time()-start1)
print ">>>   imported everything after %.1f seconds"%(time.time()-start0)
print ">>> "

# PT & ETA
ptvals  = [ 10., 27., 34., 60., 156., 223., 410., 560. ]
etavals = [ 0.0, 0.5, 1.1, 1.9, 2.0, 2.3, 2.4, 2.8, 3.4 ]
etavals = [-eta for eta in reversed(etavals)]+etavals
points  = [ ]
for pt in ptvals:
  for eta in etavals:
    #points.append((pt,-eta))
    points.append((pt,eta))

def printMatrix(name,method):
  start2 = time.time()
  print ">>>   %s:"%name
  print ">>>   %10s"%(' ')+' '.join('%10.2f'%eta for eta in etavals)
  for pt in ptvals:
    print ">>>    %10.2f"%(pt)+' '.join('%10.3f'%method(pt,eta) for eta in etavals)
  print ">>>   got %d SFs in %.1f seconds"%(len(ptvals)*len(etavals),time.time()-start2)


## TRIGGER (Muon POG)
#start1 = time.time()
#print ">>> initializing trigger SFs from Muon POG..."
#sftool_mu_trig_POG = ScaleFactor(path+"EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","IsoMu27_PtEtaBins/abseta_pt_ratio",'mu_trig',ptvseta=True)
#print ">>>   initialized in %.1f seconds"%(time.time()-start1)
#
### ID (Muon POG)
##sftool_mu_idiso_POG = ScaleFactor(path+"EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","IsoMu27_PtEtaBins/pt_abseta_ratio",'mu_trig')
#
## TRIGGER (HTT)
#start1 = time.time()
#print ">>> initializing trigger SFs from HTT..."
#sftool_mu_trig_HTT = ScaleFactorHTT(pathHTT+"Muon_IsoMu24orIsoMu27.root","ZMass",'mu_idiso')
#print ">>>   initialized in %.1f seconds"%(time.time()-start1)
# 
## ID ISO (HTT)
#start1 = time.time()
#print ">>> initializing idiso SFs from HTT..."
#sftool_mu_idiso_HTT = ScaleFactorHTT(pathHTT+"Muon_IdIso_IsoLt0p15_eff_RerecoFall17.root","ZMass",'mu_idiso')
#print ">>>   initialized in %.1f seconds"%(time.time()-start1)

sftool_ele_idiso_HTT = ScaleFactorHTT(pathHTT+"Electron_IdIso_IsoLt0.15_IsoID_eff.root","ZMass",'ele_idiso')
printMatrix('idiso HTT',sftool_ele_idiso_HTT.getSF)

sftool_ele_reco_HTT = ScaleFactor(path+"egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","EGamma_SF2D",'ele_reco',ptvseta=True)
printMatrix('reco Muon POG',sftool_ele_reco_HTT.getSF)
print ">>> "
# 
# # MUON SFs
# start1 = time.time()
# print ">>> initializing MuonSF object..."
# muonSFs = MuonSFs()
# print ">>>   initialized in %.1f seconds"%(time.time()-start1)
# 
# # GET SFs
# printMatrix('trigger',muonSFs.getTriggerSF)
# printMatrix('idiso',muonSFs.getIdIsoSF)


# MUON SFs
print ">>> "
start1 = time.time()
print ">>> initializing ElectronSFs object..."
eleSFs = ElectronSFs()
print ">>>   initialized in %.1f seconds"%(time.time()-start1)

# GET SFs
printMatrix('trigger',eleSFs.getTriggerSF)
printMatrix('idiso',eleSFs.getIdIsoSF)

print ">>> "
print ">>> done after %.1f seconds"%(time.time()-start0)
print
