#!/usr/bin/env python
import os, sys
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--infiles', dest='infiles', action='store', type=str, default=[ ])
parser.add_argument('-c', '--channel', dest='channel', action='store', choices=['tautau','mutau','eletau','muele','mumu'], type=str, default='tautau')
parser.add_argument('-t', '--type',    dest='type', action='store', choices=['data','mc'], default='mc')
parser.add_argument('-y', '--year',    dest='year', action='store', choices=[2016,2017,2018], type=int, default=2017)
parser.add_argument('-T', '--tes',     dest='tes', action='store', type=float, default=1.0)
parser.add_argument('-l', '--tag',     dest='tag', action='store', type=str, default="")
parser.add_argument('-Z', '--doZpt',   dest='doZpt', action='store_true', default=False)
args = parser.parse_args()

channel  = args.channel
year     = args.year
dataType = args.type
infiles  = args.infiles
if args.tag and args.tag[0]!='_': args.tag = '_'+args.tag
postfix  = channel + args.tag + '.root'
kwargs = {
  'year':  args.year,
  'tes':   args.tes,
  'doZpt': args.doZpt,
}

if channel=='etau':
  channel = 'eletau'
elif channel=='elemu':
  channel = 'muele'


if isinstance(infiles,str):
  infiles = infiles.split(',')
if infiles:
  dataType = 'mc'
  if infiles[0].find("/SingleMuon/")!=-1 or infiles[0].find("/Tau/")!=-1 or infiles[0].find("/SingleElectron/")!=-1:
    dataType = 'data'
else:
  if dataType=='data':
    if year==2016:
      if channel=='mutau':
        infiles = [
          ###'root://cms-xrd-global.cern.ch//store/data/Run2016B_ver1/SingleMuon/NANOAOD/Nano14Dec2018_ver1-v1/90000/87B8F064-C966-FD4F-BD32-E1FCB470AC7B.root',
          'root://cms-xrd-global.cern.ch//store/data/Run2016E/SingleMuon/NANOAOD/Nano14Dec2018-v1/10000/9DD2B819-9BE0-844E-8E12-6C0BCCBAF0E3.root',
          ###'root://cms-xrd-global.cern.ch//store/data/Run2016G/SingleMuon/NANOAOD/Nano14Dec2018-v1/20000/4BBB9D59-04A9-6241-9A92-C00B53E6D241.root',
          ###'root://cms-xrd-global.cern.ch//store/data/Run2016H/SingleMuon/NANOAOD/Nano14Dec2018-v1/90000/C6500BE6-28D9-7E46-A520-78C627F2E67E.root',
        ]
      else:
        infiles = [
          ###'root://cms-xrd-global.cern.ch//store/data/Run2016B_ver1/Tau/NANOAOD/Nano14Dec2018_ver1-v1/20000/5CB7D478-B9E8-6A48-8393-EC6C6412B2D3.root',
          'root://cms-xrd-global.cern.ch//store/data/Run2016E/Tau/NANOAOD/Nano14Dec2018-v1/10000/0ED24CE7-B408-BD48-8551-BFB7EAD2F4EB.root',
          ###'root://cms-xrd-global.cern.ch//store/data/Run2016E/Tau/NANOAOD/Nano14Dec2018-v1/10000/79A93CBE-80B7-8F4C-ADAF-0B97E1049397.root',
          ###'root://cms-xrd-global.cern.ch//store/data/Run2016H/Tau/NANOAOD/Nano14Dec2018-v1/90000/F39817D2-E1F1-4046-B45A-65B42945C313.root',
        ]
    elif year==2017:
      if channel=='mutau':
        infiles = [
          'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleMuon/NANOAOD/31Mar2018-v1/100000/EC085F70-1A64-E811-A184-003048F2E8C2.root', #    4686
          'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleMuon/NANOAOD/31Mar2018-v1/30000/3666FB57-B255-E811-BC5F-1866DA7F7AC2.root',  #   43637
          ###'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleMuon/NANOAOD/31Mar2018-v1/30000/DCA41AD0-A752-E811-A430-44A8420CC940.root',  # 1878250
        ]
      else:
        infiles = [
          'root://cms-xrd-global.cern.ch//store/data/Run2017B/Tau/NANOAOD/31Mar2018-v1/10000/04463969-D044-E811-8DC1-0242AC130002.root'
        ]
    elif year==2018: 
        infiles = [
          'root://cms-xrd-global.cern.ch//store/group/phys_tau/ProdNanoAODv4Priv/16dec18/SingleMuon/Run2018A-from_17Sep2018_ver2-NanoAODv4Priv/181216_124906/0000/myNanoRunData2018ABC_NANO_352.root',
          ###'root://cms-xrd-global.cern.ch//store/group/phys_tau/ProdNanoAODv4Priv/16dec18/SingleMuon/Run2018B-from_17Sep2018_ver1-NanoAODv4Priv/181216_124922/0000/myNanoRunData2018ABC_NANO_212.root',
          ###'root://cms-xrd-global.cern.ch//store/group/phys_tau/ProdNanoAODv4Priv/16dec18/SingleMuon/Run2018D-from_PromptReco_ver2-NanoAODv4Priv/181216_124954/0000/myNanoRunData2018D_NANO_991.root',
        ]
  else:
    if year==2016:
        infiles = [
          'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext2-v1/40000/9655C99D-1B7F-0D4C-AF8F-14B03E3732C8.root',
          ###'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext2-v1/40000/A1336349-9959-7A41-8048-CA9D272EEF7B.root',
          ###'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext2-v1/40000/58871BE4-70A5-1A4A-BE25-9C91CBC10017.root',
        ]
    elif year==2017:
#       infiles = [ 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/0A5AB04B-4B42-E811-AD7F-A4BF0112BDE6.root']
#       infiles = [ 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M500/nanoAOD/v1/nanoAOD_LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M500_1602.root']
        infiles = [
          #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B28E4243-3245-E811-B18F-001E67E71BAA.root', #  117636
          'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/90000/54362765-8948-E811-9E4A-001F29085CDE.root', #   84396
          #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/70000/002CA305-5C44-E811-9310-E0071B73B6E0.root', # 1513969
          #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B82F7904-0547-E811-BBB6-44A842CF05CC.root', 
          #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A0397F3A-7A44-E811-843F-5065F3812271.root', 
          #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/90000/1CEBB44E-4548-E811-8798-A4BF0115951C.root',
        ]
#       infiles = [ 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000_1036.root',
#                   'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000_105.root',
#                   'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000_1059.root',
#       ]
#       infiles = [ 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M500/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M500_133.root']
    elif year==2018:
      infiles = ['root://cms-xrd-global.cern.ch//store/group/phys_tau/ProdNanoAODv4Priv/16dec18/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv4Priv-from_102X_upgrade2018_realistic_v15_ver1/181216_125011/0000/myNanoRunMc2018_NANO_101.root']      
      ###infiles = ['root://cms-xrd-global.cern.ch//store/group/phys_tau/ProdNanoAODv4Priv/16dec18/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv4Priv-from_102X_upgrade2018_realistic_v15_ver2/181216_125027/0000/myNanoRunMc2018_NANO_75.root']


print ">>> %-10s = %s"%('channel',channel)
print ">>> %-10s = %s"%('dataType',dataType)
print ">>> %-10s = %s"%('year',kwargs['year'])
print ">>> %-10s = %s"%('tes',kwargs['tes'])
print ">>> %-10s = %s"%('postfix',postfix)

if channel=='tautau':
    from TauTauModule import *
    module2run = lambda : TauTauProducer(postfix, dataType, **kwargs)
elif channel=='mutau':
    from MuTauModule import *
    module2run = lambda : MuTauProducer(postfix, dataType, **kwargs)
elif channel=='eletau':
    from EleTauModule import *
    module2run = lambda : EleTauProducer(postfix, dataType, **kwargs)
elif channel=='mumu':
    from MuMuModule import *
    module2run = lambda : MuMuProducer(postfix, dataType, **kwargs)
elif channel=='muele':
    from MuEleModule import *
    module2run = lambda : MuEleProducer(postfix, dataType, **kwargs)
else:
    print 'Invalid channel name'

#p = PostProcessor(".",["../../../crab/WZ_TuneCUETP8M1_13TeV-pythia8.root"],"Jet_pt>150","keep_and_drop.txt",[exampleModule()],provenance=True)
p = PostProcessor(".", infiles, None, "keep_and_drop.txt", noOut=True, modules=[module2run()], provenance=False, postfix=postfix)
#p = PostProcessor(".",infiles,None,"keep_and_drop.txt",noOut=True, modules=[module2run()],provenance=False, jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt', postfix=postfix)

p.run()
