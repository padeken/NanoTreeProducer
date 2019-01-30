#! /usr/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from argparse import ArgumentParser

def ensureDir(directory):
  if not os.path.exists(directory):
    os.makedirs(directory)

infiles = "root://cms-xrd-global.cern.ch//store/user/arizzi/Nano01Fall17/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X-Nano01Fall17/180205_160029/0000/test94X_NANO_70.root"

parser = ArgumentParser()
parser.add_argument('-i', '--infiles', dest='infiles', action='store', type=str, default=infiles)
parser.add_argument('-o', '--outdir',  dest='outdir', action='store', type=str, default="outdir")
parser.add_argument('-N', '--outfile', dest='outfile', action='store', type=str, default="noname")
parser.add_argument('-n', '--nchunck', dest='nchunck', action='store', type=int, default='test')
parser.add_argument('-c', '--channel', dest='channel', action='store', choices=['tautau','mutau','eletau','muele','mumu'], type=str, default='tautau')
parser.add_argument('-t', '--type',    dest='type', action='store', choices=['data','mc'], default='mc')
parser.add_argument('-y', '--year',    dest='year', action='store', choices=[2016,2017,2018], type=int, default=2017)
parser.add_argument('-T', '--tes',     dest='tes', action='store', type=float, default=1.0)
args = parser.parse_args()

channel  = args.channel
dataType = args.type
infiles  = args.infiles
outdir   = args.outdir
outfile  = args.outfile
nchunck  = args.nchunck
year     = args.year
kwargs   = {
  'year':  args.year,
  'tes':   args.tes,
}

if isinstance(infiles,str):
  infiles = infiles.split(',')

if channel=='etau':
  channel = 'eletau'
elif channel=='elemu':
  channel = 'muele'

dataType = 'mc'
if infiles[0].find("/SingleMuon/")>0 or infiles[0].find("/Tau/")>0 or infiles[0].find("/SingleElectron/")>0:
  dataType = 'data'

JSON = '/shome/ineuteli/analysis/LQ_legacy/NanoTreeProducer/json/'
if year==2016:
  json = JSON+'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
  #json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
elif year==2017:
  json = JSON+'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
  #json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
else:
  json = JSON+'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
  #json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

ensureDir(outdir)

outfile = "%s_%s_%s.root"%(outfile,nchunck,channel)
postfix = "%s/%s"%(outdir,outfile)

print '-'*80
print "%-12s = %s"%('input files',infiles)
print "%-12s = %s"%('output directory',outdir)
print "%-12s = %s"%('output file',outfile)
print "%-12s = %s"%('chunck',nchunck)
print "%-12s = %s"%('channel',channel)
print "%-12s = %s"%('dataType',dataType)
print "%-12s = %s"%('year',kwargs['year'])
print "%-12s = %s"%('tes',kwargs['tes'])
print '-'*80

module2run = None
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
    module2run = lambda : MuEleProducer(postfix, dataType)
else:
    print 'Unkonwn channel !!!'
    sys.exit(0)

if dataType=='data':
    p = PostProcessor(outdir, infiles, None, "keep_and_drop.txt", noOut=True, 
                      modules=[module2run()], provenance=False, fwkJobReport=False,
                      jsonInput=json, postfix=postfix)
else:
    p = PostProcessor(outdir, infiles, None, "keep_and_drop.txt", noOut=True,
                      modules=[module2run()], provenance=False, fwkJobReport=False, postfix=postfix)

p.run()
print "DONE"
