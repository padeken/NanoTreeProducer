#!/usr/bin/env python
import os, sys
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import optparse

parser = optparse.OptionParser()
parser.add_option('-c', '--channel', action="store", type="string", default='tautau', dest='channel')
parser.add_option('-t', '--type', action="store", choices=['data','mc'], default='mc', dest='type')
(options, args) = parser.parse_args() 
channel = options.channel
if channel=='etau':
 channel = 'eletau'
print 'channel  =', channel
#channel = 'muele'

DataType = options.type
#DataType = 'data'

print 'channel = ', channel 
print 'DataType = ', DataType

if DataType=='data':
    filelist = ['root://cms-xrd-global.cern.ch//store/data/Run2017B/Tau/NANOAOD/31Mar2018-v1/10000/04463969-D044-E811-8DC1-0242AC130002.root']
else:
#     filelist = [ 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/10000/0A5AB04B-4B42-E811-AD7F-A4BF0112BDE6.root']
#     filelist = [ 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/82C67179-0942-E811-9BA7-001E67FA3920.root']
#     filelist = [ 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/5E621211-8B42-E811-9903-001E67F8FA2E.root']
#     filelist = [ 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M500/nanoAOD/v1/nanoAOD_LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M500_1602.root']
    filelist = [ 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B28E4243-3245-E811-B18F-001E67E71BAA.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/70000/002CA305-5C44-E811-9310-E0071B73B6E0.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B82F7904-0547-E811-BBB6-44A842CF05CC.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/96D73B25-4E44-E811-9491-E0071B740D80.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A0397F3A-7A44-E811-843F-5065F3812271.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E4D51829-5745-E811-92DD-E0071B74AC00.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/3405C7E7-8444-E811-AE86-44A842CFD5F2.root', 
                 #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/90000/1CEBB44E-4548-E811-8798-A4BF0115951C.root',
                ]
#     filelist = [ 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000_1036.root',
#                  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000_105.root',
#                  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000_1059.root',
#                 ]
#     filelist = [ 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M500/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M500_133.root']


_postfix = channel + '.root'

if channel=='tautau':
    from TauTauModule import *
    module2run = lambda : TauTauProducer(_postfix, DataType)
elif channel=='mutau':
    from MuTauModule import *
    module2run = lambda : MuTauProducer(_postfix, DataType)
elif channel=='eletau':
    from EleTauModule import *
    module2run = lambda : EleTauProducer(_postfix, DataType)
elif channel=='mumu':
    from MuMuModule import *
    module2run = lambda : MuMuProducer(_postfix, DataType)
elif channel=='muele':
    from MuEleModule import *
    module2run = lambda : MuEleProducer(_postfix, DataType)
else:
    print 'Invalid channel name'

#p = PostProcessor(".",["../../../crab/WZ_TuneCUETP8M1_13TeV-pythia8.root"],"Jet_pt>150","keep_and_drop.txt",[exampleModule()],provenance=True)
p = PostProcessor(".",filelist,None,"keep_and_drop.txt",noOut=True, modules=[module2run()],provenance=False, postfix=_postfix)
#p = PostProcessor(".",filelist,None,"keep_and_drop.txt",noOut=True, modules=[module2run()],provenance=False, jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt', postfix=_postfix)

p.run()
