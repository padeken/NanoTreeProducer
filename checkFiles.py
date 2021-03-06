#! /usr/bin/env python

import os, glob, sys, shlex, re
#import time
from fnmatch import fnmatch
import subprocess
from argparse import ArgumentParser
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TFile, TTree, TH1, Double

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if __name__ == '__main__':    
    description = '''Check if the job output files are valid, compare the number of events to DAS (-d), hadd them into one file per sample (-m), and merge datasets (-a).'''
    parser = ArgumentParser(prog="checkFiles",description=description,epilog="Good luck!")
    parser.add_argument('-y', '--year',     dest='years', choices=[2016,2017,2018], type=int, nargs='+', default=[2017], action='store',
                                            help="select year" )
    parser.add_argument('-c', '--channel',  dest='channels', choices=['mutau','eletau','tautau','mumu'], nargs='+', default=['tautau'], action='store' )
    parser.add_argument('-m', '--make',     dest='make', default=False, action='store_true',
                                            help="hadd all output files" )
    parser.add_argument('-a', '--hadd',     dest='haddother', default=False, action='store_true',
                                            help="hadd some samples into one (e.g. all data sets, or the extensions)" )
    parser.add_argument('-d', '--das',      dest='compareToDas', default=False, action='store_true',
                                            help="compare number of events in output to das" )
    parser.add_argument('-D', '--das-ex',   dest='compareToDasExisting', default=False, action='store_true',
                                            help="compare number of events in existing output to das" )
    parser.add_argument('-C', '--check-ex', dest='checkExisting', default=False, action='store_true',
                                            help="check existing output (e.g. 'LHE_Njets')" )
    parser.add_argument('-f', '--force',    dest='force', default=False, action='store_true',
                                            help="overwrite existing hadd'ed files" )
    parser.add_argument('-r', '--clean',    dest='cleanup', default=False, action='store_true',
                                            help="remove all output files after hadd" )
    parser.add_argument('-o', '--outdir',   dest='outdir', type=str, default=None, action='store' )
    parser.add_argument('-s', '--sample',   dest='samples', type=str, nargs='+', default=[ ], action='store',
                                            help="samples to run over, glob patterns (wildcards * and ?) are allowed." )
    parser.add_argument('-x', '--veto',     dest='veto', action='store', type=str, default=None,
                                            help="veto this sample" )
    parser.add_argument('-t', '--type',     dest='type', choices=['data','mc'], type=str, default=None, action='store',
                                            help="filter data or MC to submit" )
    parser.add_argument('-l', '--tag',      dest='tag', type=str, default="", action='store',
                                            help="add a tag to the output file" )
    parser.add_argument('-v', '--verbose',  dest='verbose', default=False, action='store_true',
                                            help="set verbose" )
    args = parser.parse_args()

subdirs = [ 'TT', 'DY', 'W*J', 'ST', 'LQ', 'Tau', 'SingleMuon', 'SingleElectron' ]
sample_dict = [
   ('DY',             "DYJetsToLL_M-10to50",              "DYJetsToLL_M-10to50_Tune*madgraph*pythia8"                ),
   ('DY',             "DYJetsToLL_M-50_ext",              "DYJetsToLL_M-50_Tune*madgraph*pythia8/*ext1"              ), # ext before reg !
   ('DY',             "DYJetsToLL_M-50_reg",              "DYJetsToLL_M-50_Tune*madgraph*pythia8"                    ),
   ('DY',             "DY1JetsToLL_M-50",                 "DY1JetsToLL_M-50_Tune*madgraph*pythia8"                   ),
   ('DY',             "DY2JetsToLL_M-50",                 "DY2JetsToLL_M-50_Tune*madgraph*pythia8"                   ),
   ('DY',             "DY3JetsToLL_M-50",                 "DY3JetsToLL_M-50_Tune*madgraph*pythia8"                   ),
   ('DY',             "DY4JetsToLL_M-50",                 "DY4JetsToLL_M-50_Tune*madgraph*pythia8"                   ),
   ('WJ',             "WJetsToLNu_ext",                   "WJetsToLNu_Tune*madgraph*pythia8/*ext"                    ), # ext before reg !
   ('WJ',             "WJetsToLNu_reg",                   "WJetsToLNu_Tune*madgraph*pythia8/RunIISummer16"           ),
   ('WJ',             "WJetsToLNu",                       "WJetsToLNu_Tune*madgraph*pythia8/RunIIFall17"             ),
   ('WJ',             "WJetsToLNu",                       "WJetsToLNu_Tune*madgraph*pythia8"                         ),
   ('WJ',             "W1JetsToLNu",                      "W1JetsToLNu_Tune*madgraph*pythia8"                        ),
   ('WJ',             "W2JetsToLNu",                      "W2JetsToLNu_Tune*madgraph*pythia8"                        ),
   ('WJ',             "W3JetsToLNu",                      "W3JetsToLNu_Tune*madgraph*pythia8"                        ),
   ('WJ',             "W4JetsToLNu",                      "W4JetsToLNu_Tune*madgraph*pythia8"                        ),
   ('ST',             "ST_t-channel_antitop",             "ST_t-channel_antitop_4f_inclusiveDecays"                  ),
   ('ST',             "ST_t-channel_top",                 "ST_t-channel_top_4f_inclusiveDecays"                      ),
   ('ST',             "ST_tW_antitop",                    "ST_tW_antitop_5f_inclusiveDecays"                         ),
   ('ST',             "ST_tW_top",                        "ST_tW_top_5f_inclusiveDecays"                             ),
   ('ST',             "ST_s-channel",                     "ST_s-channel_4f_hadronicDecays"                           ),
   #('TT',             "TT_ext",                           "TT_Tune*powheg*pythia8/*ext"                              ),
   ('TT',             "TT",                               "TT_Tune*powheg*pythia8"                                   ),
   ('TT',             "TTTo2L2Nu",                        "TTTo2L2Nu"                                                ),
   ('TT',             "TTToHadronic",                     "TTToHadronic"                                             ),
   ('TT',             "TTToSemiLeptonic",                 "TTToSemiLeptonic"                                         ),
   ('VV',             "WWTo1L1Nu2Q",                      "WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8"           ),
   ('VV',             "WWTo2L2Nu",                        "WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8" ),
   ('VV',             "WWToLNuQQ_ext",                    "WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/*ext"      ), # ext before reg !
   ('VV',             "WWToLNuQQ_reg",                    "WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8"           ),
   ('VV',             "WW",                               "WW_Tune*pythia8"                                          ),
   ('VV',             "WZ",                               "WZ_Tune*pythia8"                                          ),
   ('VV',             "ZZ",                               "ZZ_Tune*pythia8"                                          ),
   ('Tau',            "Tau_$RUN",                         "Tau/$RUN"                                                 ),
   ('Tau',            "Tau_Run2017B",                     "Tau/ytakahas-NanoTest_20180507_B"                         ),
   ('Tau',            "Tau_Run2017C",                     "Tau/ytakahas-NanoTest_20180507_C"                         ),
   ('Tau',            "Tau_Run2017D",                     "Tau/ytakahas-NanoTest_20180507_D"                         ),
   ('Tau',            "Tau_Run2017E",                     "Tau/ytakahas-NanoTest_20180507_E"                         ),
   ('Tau',            "Tau_Run2017F",                     "Tau/ytakahas-NanoTest_20180507_F"                         ),
   ('Tau',            "Tau_$RUN",                         "Tau/manzoni-$RUN"                                         ),
   ('SingleMuon',     "SingleMuon_$RUN",                  "SingleMuon/$RUN"                                          ),
   ('SingleMuon',     "SingleMuon_$RUN",                  "SingleMuon/manzoni-$RUN"                                  ),
   ('SingleElectron', "SingleElectron_$RUN",              "SingleElectron/$RUN"                                      ),
   ('SingleElectron', "SingleElectron_Run2017B",          "SingleElectron/ytakahas-Nano_SingleElectron_20180507_B"   ),
   ('SingleElectron', "SingleElectron_Run2017C",          "SingleElectron/ytakahas-Nano_SingleElectron_20180507_C"   ),
   ('SingleElectron', "SingleElectron_Run2017D",          "SingleElectron/ytakahas-Nano_SingleElectron_20180507_D"   ),
   ('SingleElectron', "SingleElectron_Run2017E",          "SingleElectron/ytakahas-Nano_SingleElectron_20180507_E"   ),
   ('SingleElectron', "SingleElectron_Run2017F",          "SingleElectron/ytakahas-Nano_SingleElectron_20180507_F"   ),
   ('LQ',             "LQ3ToTauB_t-channel_M$MASS",       "LQ3ToTauB_Fall2017_5f_Madgraph_LO_t-channel-M$MASS"       ),
   ('LQ',             "LQ3ToTauB_s-channel_M$MASS",       "LQ3ToTauB_Fall2017_5f_Madgraph_LO_s-channel-M$MASS"       ),
   ('LQ',             "LQ3ToTauB_pair_M$MASS",            "LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M$MASS"            ),
   ('LQ',             "VectorLQ3ToTauB_s-channel_M$MASS", "VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_s_channel_M$MASS" ),
   ('LQ',             "VectorLQ3ToTauB_pair_M$MASS",      "VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M$MASS"      ),
]
sample_dict = [(d,s,p.replace('*','.*').replace('$MASS','(\d+)').replace('$RUN','(Run201\d[A-H])')) for d,s,p in sample_dict] # convert to regex pattern
#sample_dict = { k: v.lstrip('/').replace('/','__') for k, v in sample_dict.iteritems() }
haddsets = [
  ('DY',             "DYJetsToLL_M-50",      [ 'DYJetsToLL_M-50_*'    ]),
  ('WJ',             "WJetsToLNu",           [ 'WJetsToLNu_*'         ]),
  ('Tau',            "Tau_$RUN",             [ 'Tau_$RUN?'            ]),
  ('SingleMuon',     "SingleMuon_$RUN",      [ 'SingleMuon_$RUN?'     ]),
  ('SingleElectron', "SingleElectron_$RUN",  [ 'SingleElectron_$RUN?' ]),
]



def main(args):
  
  years      = args.years
  channels   = args.channels
  tag        = args.tag
  if len(tag)>0 and '_' not in tag[0]:
    tag = '_'+tag
  
  for year in years:
    indir      = "output_%s/"%(year)
    samplesdir = args.outdir if args.outdir else "/scratch/ineuteli/analysis/LQ_%s"%(year)
    os.chdir(indir)
    
    # CHECK EXISTING
    if args.checkExisting:
      for channel in channels:
        infiles  = "%s/*/*_%s.root"%(channel)
        filelist = glob.glob(infiles)
        pattern  = infiles.split('/')[-1]
        for file in filelist:
          if not isValidSample(pattern): continue
          checkFiles(file,pattern)
      continue
    
    # GET LIST
    samplelist = [ ]
    for directory in sorted(os.listdir('./')):
      if not os.path.isdir(directory): continue
      if not isValidSample(directory): continue
      samplelist.append(directory)
    if not samplelist:
      print "No samples found in %s!"%(indir)
    if args.verbose:
      print 'samplelist = %s\n'%(samplelist)
    
    # CHECK samples
    for channel in channels:
      print header(year,channel)
      
      # HADD samples
      if not args.haddother or args.make:
        for directory in samplelist:
            if args.verbose:
              print directory
            
            subdir, samplename = getSampleShortName(directory)
            outdir  = "%s/%s"%(samplesdir,subdir)
            outfile = "%s/%s_%s%s.root"%(outdir,samplename,channel,tag)
            infiles = '%s/*_%s.root'%(directory,channel)
            
            if args.verbose:
              print "directory = %s"%(directory)
              print "outdir    = %s"%(outdir)
              print "outfile   = %s"%(outfile)
              print "infiles   = %s"%(infiles)
            
            #if directory.find('W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8__ytakahas-NanoTest_20180507_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be__USER')==-1: continue
            filelist = glob.glob(infiles)
            if not filelist: continue
            
            if checkFiles(filelist,directory):
              print bcolors.BOLD + bcolors.OKGREEN + '[OK] ' + directory + ' ... can be hadded ' + bcolors.ENDC
            
            if 'LQ3' not in directory:
              if args.compareToDas:
                compareEventsToDAS(filelist,directory)
              if args.compareToDasExisting and os.path.isfile(outfile):
                print '   check existing file %s:'%(outfile)
                compareEventsToDAS(outfile,directory)
              #else:
              #  print bcolors.BOLD + bcolors.OKBLUE + '   [OK] ' + directory + bcolors.ENDC
              #print
            
            # HADD
            if args.make:
                ensureDirectory(outdir)
                if os.path.isfile(outfile):
                  if args.force:
                    print bcolors.BOLD + bcolors.WARNING + "   [WN] target %s already exists! Overwriting..."%(outfile) + bcolors.ENDC
                  else:
                    print bcolors.BOLD + bcolors.FAIL + "   [NG] target %s already exists! Use --force or -f to overwrite."%(outfile) + bcolors.ENDC
                    continue
                
                haddcmd = 'hadd -f %s %s'%(outfile,infiles)
                print haddcmd
                os.system(haddcmd)
                
                if 'LQ3' not in directory:
                     compareEventsToDAS(outfile,directory)
                #    skimcmd = 'python extractTrees.py -c %s -f %s'%(channel,outfile)
                #    rmcmd = 'rm %s'%(infiles)
                #    #os.system(skimcmd)
                #    #os.system(rmcmd)
                #    continue
                
                #skimcmd = 'python extractTrees.py -c %s -f %s'%(channel,outfile)
                #os.system(skimcmd)
                
                # CLEAN UP
                if args.cleanup:
                  rmcmd = 'rm %s'%(infiles)
                  print bcolors.BOLD + bcolors.OKBLUE + "   removing %d output files..."%(len(infiles)) + bcolors.ENDC
                  if verbose:
                    print rmcmd
                  os.system(rmcmd)
                print
      
      # HADD other
      if args.haddother:
        for subdir, samplename, sampleset in haddsets:
            if args.verbose:
              print subdir, samplename, sampleset
            if args.samples and not matchSampleToPattern(samplename,args.samples): continue
            if args.veto and matchSampleToPattern(directory,args.veto): continue
            if 'SingleMuon' in subdir and channel not in ['mutau','mumu']: continue
            if 'SingleElectron' in subdir and channel!='eletau': continue
            if 'Tau' in subdir and channel!='tautau': continue
            if 'LQ3' in subdir and channel not in ['mutau','tautau','eletau']: continue
            if '2017' in samplename and year!=2017: continue
            if '2018' in samplename and year!=2018: continue
            if '$RUN' in samplename:
              samplename = samplename.replace('$RUN','Run%d'%year)
              sampleset  = [s.replace('$RUN','Run%d'%year) for s in sampleset]
            
            outdir  = "%s/%s"%(samplesdir,subdir)
            outfile = "%s/%s_%s%s.root"%(outdir,samplename,channel,tag)
            infiles = ['%s/%s_%s%s.root'%(outdir,s,channel,tag) for s in sampleset] #.replace('ele','e')
            ensureDirectory(outdir)
            
            # OVERWRITE ?
            if os.path.isfile(outfile):
              if args.force:
                pass
                #if args.verbose:
                #  print bcolors.BOLD + bcolors.WARNING + "[WN] target %s already exists! Overwriting..."%(outfile) + bcolors.ENDC
              else:
                print bcolors.BOLD + bcolors.FAIL + "[NG] target %s already exists! Use --force or -f to overwrite."%(outfile) + bcolors.ENDC
                continue
            
            # CHECK FILES
            allinfiles = [ ]
            for infile in infiles[:]:
              if '*' in infile or '?' in infile:
                files = glob.glob(infile)
                allinfiles += files
                if not files:
                  print bcolors.BOLD + bcolors.FAIL + '[NG] no match for the glob pattern %s! Removing pattern from hadd list for "%s"...'%(infile,samplename) + bcolors.ENDC
                  infiles.remove(infile)
              elif not os.path.isfile(infile):
                print bcolors.BOLD + bcolors.FAIL + '[NG] infile %s does not exists! Removing from hadd list for "%s"...'%(infile,samplename) + bcolors.ENDC
                infiles.remove(infile)
              else:
                allinfiles.append(infile)
            
            # HADD
            if args.verbose:
              print "infiles =", infiles
              print "allfiles =", allinfiles
            if len(allinfiles)==1:
              print bcolors.BOLD + bcolors.WARNING + "[WN] found only one file (%s) to hadd to %s!"%(allinfiles[0],outfile) + bcolors.ENDC 
            elif len(allinfiles)>1:
              print bcolors.BOLD + bcolors.OKGREEN + '[OK] hadding %s' %(outfile) + bcolors.ENDC
              haddcmd = 'hadd -f %s %s'%(outfile,' '.join(infiles))
              print haddcmd
              os.system(haddcmd)
            else:
              print bcolors.BOLD + bcolors.WARNING + "[WN] no files to hadd!" + bcolors.ENDC
            print
    os.chdir('..')
     


def isValidSample(pattern):
  if args.samples and not matchSampleToPattern(pattern,args.samples): return False
  if args.veto and matchSampleToPattern(pattern,args.veto): return False
  if args.type=='mc' and any(s in pattern[:len(s)+2] for s in ['SingleMuon','SingleElectron','Tau']): return False
  if args.type=='data' and not any(s in pattern[:len(s)+2] for s in ['SingleMuon','SingleElectron','Tau']): return False
  return True


indexpattern = re.compile(r".*_(\d+)_[a-z]+\.root")
def checkFiles(filelist,directory):
    if args.verbose:
      print "checkFiles: %s, %s"%(filelist,directory)
    if isinstance(filelist,str):
      filelist = [filelist]
    badfiles = [ ]
    ifound   = [ ]
    for filename in filelist:
      file  = TFile(filename, 'READ')
      isbad = False
      if file.IsZombie():
        print bcolors.FAIL + '[NG] file %s is a zombie'%(filename) + bcolors.ENDC
        isbad = True
      else:
        tree = file.Get('tree')
        if not isinstance(tree,TTree):
          print bcolors.FAIL + '[NG] no tree found in ' + filename + bcolors.ENDC
          isbad = True
        elif not isinstance(file.Get('cutflow'),TH1):
          print bcolors.FAIL + '[NG] no cutflow found in ' + filename + bcolors.ENDC
          isbad = True
        elif any(s in filename for s in ['DYJets','WJets']) and tree.GetMaximum('LHE_Njets')>10:
          print bcolors.WARNING + '[WN] LHE_Njets = %d > 10 in %s'%(tree.GetMaximum('LHE_Njets'),filename) + bcolors.ENDC
      if isbad:
        badfiles.append(filename)
        #rmcmd = 'rm %s' %filename
        #print rmcmd
        #os.system(rmcmd)
      file.Close()
      match = indexpattern.search(filename)
      if match: ifound.append(int(match.group(1)))
    
    if len(badfiles)>0:
      print bcolors.FAIL + "[NG] %s:   %d out of %d files have no tree!"%(directory,len(badfiles),len(filelist)) + bcolors.ENDC
      return False
    
    # TODO: check all chunks (those>imax)
    imax = max(ifound)+1
    if len(filelist)<imax:
      imiss = [ i for i in range(0,max(ifound)) if i not in ifound ]
      chunktext = ('chunks ' if len(imiss)>1 else 'chunk ') + ', '.join(str(i) for i in imiss)
      print bcolors.BOLD + bcolors.WARNING + "[WN] %s missing %d/%d files (%s) ?"%(directory,len(imiss),len(filelist),chunktext) + bcolors.ENDC
    
    return True
    
def compareEventsToDAS(filenames,dasname):
    """Compare a number of processed events in an output file to the available number of events in DAS."""
    dasname = dasname.replace('__', '/')
    if args.verbose:
      print "compareEventsToDAS: %s, %s"%(filenames,dasname)
      #start = time.time()
    if isinstance(filenames,str):
      filenames = [filenames]
    total_processed = 0
    nfiles = len(filenames)
    for filename in filenames:
      file = TFile(filename, 'READ')
      if file.IsZombie():
        continue
      hist = file.Get('cutflow')
      if hist:
        events_processed = hist.GetBinContent(1)
        if args.verbose:
          print "%12d events processed in %s "%(events_processed,filename)
        total_processed += events_processed
      #else:
      #  print bcolors.FAIL + '[NG] compareEventsToDAS: no cutflow found in ' + filename + bcolors.ENDC
      file.Close()
    
    instance = 'prod/phys03' if 'USER' in dasname else 'prod/global'
    dascmd   = 'das_client --limit=0 --query=\"summary dataset=/%s instance=%s\"'%(dasname,instance)
    if args.verbose:
      print dascmd
    dasargs  = shlex.split(dascmd)
    output, error = subprocess.Popen(dasargs, stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    
    if not "nevents" in output:
        print bcolors.BOLD + bcolors.FAIL + '   [NG] Did not find nevents for "%s" in DAS. Return message:'%(dasname) + bcolors.ENDC 
        print bcolors.FAIL + '     ' + output + bcolors.ENDC
        return False
    total_das = Double(output.split('"nevents":')[1].split(',')[0])
    fraction = total_processed/total_das
    
    nfiles = ", %d files"%(nfiles) if nfiles>1 else ""
    if fraction > 1.001:
        print bcolors.BOLD + bcolors.FAIL + '   [NG] DAS entries = %d, Processed in tree = %d (frac = %.2f > 1%s)'%(total_das,total_processed,fraction,nfiles) + bcolors.ENDC
    elif fraction > 0.8:
        print bcolors.BOLD + bcolors.OKBLUE + '   [OK] DAS entries = %d, Processed in tree = %d (frac = %.2f%s)'%(total_das,total_processed,fraction,nfiles) + bcolors.ENDC
    else:
        print bcolors.BOLD + bcolors.FAIL + '   [NG] DAS entries = %d, Processed in tree = %d (frac = %.2f < 0.8%s)'%(total_das,total_processed,fraction,nfiles) + bcolors.ENDC
    return True
    
def getSampleShortName(dasname):
  """Get short subdir and sample name from sample_dict."""
  #if '__nanoaod' in dasname.lower():
  #  dasname = dasname[:dasname.lower().index('__nanoaod')]
  #if '__user' in dasname.lower():
  #  dasname = dasname[:dasname.lower().index('__user')]
  if args.verbose:
    print "getSampleShortName: %s"%(dasname)
  dasname = dasname.replace('__','/').lstrip('/')
  for subdir, samplename, pattern in sample_dict:
    matches = re.findall(pattern,dasname)
    if matches:
      samplename = samplename.replace('$MASS',matches[0]).replace('$RUN',matches[0])
      if args.verbose:
         print "getSampleShortName: MATCH! subdir=%s, samplename=%s, pattern=%s"%(subdir, samplename, pattern)
      return subdir, samplename
  print bcolors.BOLD + bcolors.WARNING + '[WN] getSampleShortName: did not find subdir and short sample name for "%s"! Will save in subdir \'unknown\''%(dasname) + bcolors.ENDC 
  return "unknown", dasname.replace('/','__')
  
def getSubdir(dir):
  for subdir in subdirs:
    if '*' in subdir or '?' in subdir:
      if fnmatch(dir,subdir):
        return subdir
    else:
      if subdir==dir[:len(subdir)]:
        return subdir
  return "unknown"
  
def matchSampleToPattern(sample,patterns):
  """Match sample name to some pattern."""
  sample = sample.lstrip('/')
  if not isinstance(patterns,list):
    patterns = [patterns]
  for pattern in patterns:
    if '*' in pattern or '?' in pattern:
      if fnmatch(sample,pattern+'*'):
        return True
    else:
      if pattern in sample[:len(pattern)+1]:
        return True
  return False
  
def ensureDirectory(dirname):
  """Make directory if it does not exist."""
  if not os.path.exists(dirname):
    os.makedirs(dirname)
    print '>>> made directory "%s"'%(dirname)
    if not os.path.exists(dirname):
      print '>>> failed to make directory "%s"'%(dirname)
  return dirname
  
headeri = 0
def header(year,channel):
  global headeri
  title  = "%s, %s"%(year,channel)
  string = ("\n\n" if headeri>0 else "") +\
           "   ###%s\n"    % ('#'*(len(title)+3)) +\
           "   #  %s  #\n" % (title) +\
           "   ###%s\n"    % ('#'*(len(title)+3))
  headeri += 1
  return string


if __name__ == '__main__':
    
    print
    main(args)
    print
    


