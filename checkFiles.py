#! /usr/bin/env python
import os, glob, sys, shlex, re
from fnmatch import fnmatch
import subprocess
from ROOT import TFile, Double
from argparse import ArgumentParser

scratchdir = "/scratch/ineuteli/analysis/LQ_2017"
parser = ArgumentParser()
parser.add_argument('-m', '--make',    dest='make', default=False, action='store_true' )
parser.add_argument('-d', '--das',     dest='compareToDas', default=False, action='store_true' )
parser.add_argument('-f', '--force',   dest='force', default=False, action='store_true',
                                       help="overwrite existing hadd'ed files" )
parser.add_argument('-y', '--year',    dest='year', choices=[2017,2018], type=int, default=2017, action='store',
                                       help="select year" )
parser.add_argument('-a', '--hadd',    dest='haddother', default=False, action='store_true',
                                       help="hadd some samples (e.g. all data sets, or the extensions)" )
parser.add_argument('-c', '--channel', dest='channels', choices=['mutau','eletau','tautau'], nargs='+', default=["tautau"], action='store' )
parser.add_argument('-o', '--outdir',  dest='outdir', type=str, default=scratchdir, action='store' )
parser.add_argument('-t', '--tag',     dest='tag', type=str, default="", action='store' )
parser.add_argument('-s', '--samples', dest='samples', type=str, nargs='+', default=[ ], action='store',
                                       help="samples to run over, glob patterns (wildcards * and ?) are allowed." )
parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                                       help="set verbose" )
args = parser.parse_args()


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

subdirs = [ 'TT', 'DY', 'W*J', 'ST', 'LQ', 'Tau', 'SingleMuon', 'SingleElectron' ]
sample_dict = { 
   ('DY',  "DYJetsToLL_M-10to50" ): "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/ytakahas-Nano_20180518_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be",
   ('DY',  "DYJetsToLL_M-50_reg" ): "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('DY',  "DYJetsToLL_M-50_ext" ): "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1",
   ('DY',  "DY1JetsToLL_M-50"    ): "DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1",
   ('DY',  "DY2JetsToLL_M-50"    ): "DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1",
   ('DY',  "DY3JetsToLL_M-50"    ): "DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1",
   ('DY',  "DY4JetsToLL_M-50"    ): "DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('WJ',  "WJetsToLNu"          ): "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/ytakahas-NanoTest_20180507_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be",
   ('WJ',  "W1JetsToLNu"         ): "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/ytakahas-NanoTest_20180507_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be",
   ('WJ',  "W2JetsToLNu"         ): "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/ytakahas-NanoTest_20180507_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be",
   ('WJ',  "W3JetsToLNu"         ): "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/ytakahas-NanoTest_20180507_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be",
   ('WJ',  "W4JetsToLNu"         ): "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/ytakahas-NanoTest_20180507_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be",
   ('ST',  "ST_t-channel_antitop"): "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('ST',  "ST_t-channel_top"    ): "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('ST',  "ST_tW_antitop"       ): "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('ST',  "ST_tW_top"           ): "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('TT',  "TTTo2L2Nu"           ): "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('TT',  "TTToHadronic"        ): "TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('TT',  "TTToSemiLeptonic"    ): "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('VV',  "WWTo1L1Nu2Q"         ): "WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('VV',  "WWTo2L2Nu"           ): "WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1",
   ('VV',  "WWToLNuQQ_reg"       ): "WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('VV',  "WWToLNuQQ_ext"       ): "WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1",
   ('VV',  "WW"                  ): "WW_TuneCP5_13TeV-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('VV',  "WZ"                  ): "WZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('VV',  "ZZ"                  ): "ZZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1",
   ('LQ',  "LQ3ToTauB_t-channel_M$MASS"         ): "LQ3ToTauB_Fall2017_5f_Madgraph_LO_t-channel-M$MASS",
   ('LQ',  "LQ3ToTauB_s-channel_M$MASS"         ): "LQ3ToTauB_Fall2017_5f_Madgraph_LO_s-channel-M$MASS",
   ('LQ',  "LQ3ToTauB_pair_M$MASS"              ): "LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M$MASS",
   ('LQ',  "VectorLQ3ToTauB_s-channel_M$MASS"   ): "VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_s_channel_M$MASS",
   ('LQ',  "VectorLQ3ToTauB_pair_M$MASS"        ): "VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M$MASS",
   ('Tau',            "Tau_Run2017B"            ): "Tau/Run2017B-31Mar2018-v1",
   ('Tau',            "Tau_Run2017C"            ): "Tau/Run2017C-31Mar2018-v1",
   ('Tau',            "Tau_Run2017D"            ): "Tau/Run2017D-31Mar2018-v1",
   ('Tau',            "Tau_Run2017E"            ): "Tau/Run2017E-31Mar2018-v1",
   ('Tau',            "Tau_Run2017F"            ): "Tau/Run2017F-31Mar2018-v1",
   ('Tau',            "Tau_Run2017B"            ): "Tau/ytakahas-NanoTest_20180507_B-55055ff3316a022bb149a249662ed4c4",
   ('Tau',            "Tau_Run2017C"            ): "Tau/ytakahas-NanoTest_20180507_C-55055ff3316a022bb149a249662ed4c4",
   ('Tau',            "Tau_Run2017D"            ): "Tau/ytakahas-NanoTest_20180507_D-55055ff3316a022bb149a249662ed4c4",
   ('Tau',            "Tau_Run2017E"            ): "Tau/ytakahas-NanoTest_20180507_E-55055ff3316a022bb149a249662ed4c4",
   ('Tau',            "Tau_Run2017F"            ): "Tau/ytakahas-NanoTest_20180507_F-55055ff3316a022bb149a249662ed4c4",
   ('SingleMuon',     "SingleMuon_Run2017B"     ): "SingleMuon/Run2017B-31Mar2018-v1",
   ('SingleMuon',     "SingleMuon_Run2017C"     ): "SingleMuon/Run2017C-31Mar2018-v1",
   ('SingleMuon',     "SingleMuon_Run2017D"     ): "SingleMuon/Run2017D-31Mar2018-v1",
   ('SingleMuon',     "SingleMuon_Run2017E"     ): "SingleMuon/Run2017E-31Mar2018-v1",
   ('SingleMuon',     "SingleMuon_Run2017F"     ): "SingleMuon/Run2017F-31Mar2018-v1",
   ('SingleElectron', "SingleElectron_Run2017B" ): "SingleElectron/Run2017B-31Mar2018-v1",
   ('SingleElectron', "SingleElectron_Run2017C" ): "SingleElectron/Run2017C-31Mar2018-v1",
   ('SingleElectron', "SingleElectron_Run2017D" ): "SingleElectron/Run2017D-31Mar2018-v1",
   ('SingleElectron', "SingleElectron_Run2017E" ): "SingleElectron/Run2017E-31Mar2018-v1",
   ('SingleElectron', "SingleElectron_Run2017F" ): "SingleElectron/Run2017F-31Mar2018-v1",
   ('SingleElectron', "SingleElectron_Run2017B" ): "SingleElectron/ytakahas-Nano_SingleElectron_20180507_B-55055ff3316a022bb149a249662ed4c4/USER",
   ('SingleElectron', "SingleElectron_Run2017C" ): "SingleElectron/ytakahas-Nano_SingleElectron_20180507_C-55055ff3316a022bb149a249662ed4c4/USER",
   ('SingleElectron', "SingleElectron_Run2017D" ): "SingleElectron/ytakahas-Nano_SingleElectron_20180507_D-55055ff3316a022bb149a249662ed4c4/USER",
   ('SingleElectron', "SingleElectron_Run2017E" ): "SingleElectron/ytakahas-Nano_SingleElectron_20180507_E-55055ff3316a022bb149a249662ed4c4/USER",
   ('SingleElectron', "SingleElectron_Run2017F" ): "SingleElectron/ytakahas-Nano_SingleElectron_20180507_F-55055ff3316a022bb149a249662ed4c4/USER",
   #('SingleElectron', "SingleElectron_Run2017B" ): "/SingleElectron/Run2017B-31Mar2018-v1/NANOAOD",
   #('SingleElectron', "SingleElectron_Run2017C" ): "/SingleElectron/Run2017C-31Mar2018-v1/NANOAOD",
   #('SingleElectron', "SingleElectron_Run2017D" ): "/SingleElectron/Run2017D-31Mar2018-v1/NANOAOD",
   #('SingleElectron', "SingleElectron_Run2017E" ): "/SingleElectron/Run2017E-31Mar2018-v1/NANOAOD",
   #('SingleElectron', "SingleElectron_Run2017F" ): "/SingleElectron/Run2017F-31Mar2018-v1/NANOAOD",
}
#sample_dict = { k: v.lstrip('/').replace('/','__') for k, v in sample_dict.iteritems() }
haddsets = {
   ('DY',             "DYJetsToLL_M-50"        ): [ 'DYJetsToLL_M-50_*'       ],
   ('Tau',            "Tau_Run2017"            ): [ 'Tau_Run2017?'            ],
   ('SingleMuon',     "SingleMuon_Run2017"     ): [ 'SingleMuon_Run2017?'     ],
   ('SingleElectron', "SingleElectron_Run2017" ): [ 'SingleElectron_Run2017?' ],
}



def main():
  for channel in args.channels:
    
    # HADD samples
    if not args.haddother or args.make:
      for directory in sorted(os.listdir("./")):
          if not os.path.isdir(directory): continue
          if args.samples and not matchSample(args.samples,directory): continue
          
          subdir, samplename = getSampleName(directory)
          outdir  = "%s/%s"%(args.outdir,subdir)
          outfile = "%s/%s_%s.root"%(outdir,samplename,channel)
          infiles = '%s/*_%s.root'%(directory,channel)
          
          #if directory.find('W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8__ytakahas-NanoTest_20180507_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be__USER')==-1: continue
          filelist = glob.glob(infiles)
          if not filelist: continue
          
          flag  = False
          files = [ ]
          for file2check in filelist:
            file = TFile(file2check, 'READ')
            if not file.GetListOfKeys().Contains("tree"):
              files.append(file2check)
              flag = True
              #rmcmd = 'rm %s' %file2check
              #print rmcmd
              #os.system(rmcmd)
              print bcolors.FAIL + '[NG] no tree found in ' + file2check + bcolors.ENDC
            file.Close()
        
          if flag:
            print bcolors.FAIL + "[NG] " + directory + bcolors.ENDC
            print '   ', len(files), ' out of ', len(filelist), ' files have no tree!'
            #continue
          else:
            print bcolors.BOLD + bcolors.OKGREEN + '[OK] ' + directory + ' ... can be hadded ' + bcolors.ENDC
          
          if os.path.isfile(outfile):
            if args.compareToDas and directory.find('LQ3')<0:
              compareEventsToDAS(outfile,directory)
            else:
              print bcolors.BOLD + bcolors.OKBLUE + '    [OK] ' + directory + bcolors.ENDC
            print 
        
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
              #print haddcmd
              os.system(haddcmd)
              
              if directory.find('LQ3')!=-1:
                  skimcmd = 'python extractTrees.py -c %s -f %s'%(channel,outfile)
                  rmcmd = 'rm %s'%(infiles)
                  #os.system(skimcmd)
                  #os.system(rmcmd)
                  continue
              compareEventsToDAS(outfile,directory)
              
              skimcmd = 'python extractTrees.py -c %s -f %s'%(channel,outfile)
              #os.system(skimcmd)
              
              # cleaning up ...
              rmcmd = 'rm %s'%(infiles)
              #os.system(rmcmd)
              print
    
    # HADD other
    if args.haddother:
      for (subdir,samplename), sampleset in haddsets.iteritems():
          if args.samples and not matchSample(args.samples,samplename): continue
          
          outdir  = "%s/%s"%(args.outdir,subdir)
          outfile = "%s/%s_%s.root"%(outdir,samplename,channel)
          infiles = ['%s/%s_%s.root'%(outdir,s,channel) for s in sampleset] #.replace('ele','e')
          ensureDirectory(outdir)
          
          # OVERWRITE ?
          if os.path.isfile(outfile):
            if args.force:
              print bcolors.BOLD + bcolors.WARNING + "[WN] target %s already exists! Overwriting..."%(outfile) + bcolors.ENDC
            else:
              print bcolors.BOLD + bcolors.FAIL + "[NG] target %s already exists! Use --force or -f to overwrite."%(outfile) + bcolors.ENDC
              continue
          
          # CHECK FILES  
          for infile in infiles[:]:
            if '*' in infile or '?' in infile:
              if not glob.glob(infile):
                print bcolors.BOLD + bcolors.FAIL + '[NG] no match for the glob pattern %s! Removing pattern from hadd list for "%s"...'%(infile,samplename) + bcolors.ENDC
                infiles.remove(infile)
            elif not os.path.isfile(infile):
              print bcolors.BOLD + bcolors.FAIL + '[NG] infile %s does not exists! Removing from hadd list for "%s"...'%(infile,samplename) + bcolors.ENDC
              infiles.remove(infile)
          
          # HADD
          if len(infiles)>0:
            print bcolors.BOLD + bcolors.OKGREEN + '[OK] hadding %s' %(outfile) + bcolors.ENDC
            haddcmd = 'hadd -f %s %s'%(outfile,' '.join(infiles))
            print haddcmd
            os.system(haddcmd)
          else:
            print bcolors.BOLD + bcolors.WARNING + "[WN] no files to hadd!" + bcolors.ENDC
          
          print
        

def compareEventsToDAS(filename,dasname):
    """Compare a number of processed events in an output file to the available number of events in DAS."""
    dasname = dasname.replace('__', '/')
    
    print "compareEventsToDAS: %s"%(filename)
    f_hadd = TFile(filename)
    total_processed = Double(f_hadd.Get('h_cutflow').GetBinContent(1))
    #print 'Check number of events = ', total_processed
    
    instance = 'prod/global'
    if dasname.find('USER')!=-1:
        instance = 'prod/phys03'
    
    dascmd = 'das_client --limit=0 --query=\"summary dataset=/' + dasname + ' instance=' + instance + '\"'
    dasargs = shlex.split(dascmd)
    output, error = subprocess.Popen(dasargs, stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    
    if not "nevents" in output:
        print bcolors.BOLD + bcolors.FAIL + '   [NG] Did not find nevents for "%s" in DAS. Return message:'%(dasname) + bcolors.ENDC 
        print bcolors.FAIL + '     ' + output + bcolors.ENDC
        return False
    total_das = Double(output.split('"nevents":')[1].split(',')[0])
    
    fraction = total_processed/total_das
    
    if fraction > 1.001:
        print bcolors.BOLD + bcolors.FAIL + '   [NG] DAS entries = ' + str(int(total_das)) + ' Tree produced = ' + str(int(total_processed)) + ' (frac = {0:.2f}'.format(fraction) + ' > 1)' + bcolors.ENDC 
    elif fraction > 0.8:
        print bcolors.BOLD + bcolors.OKBLUE + '    [OK] DAS entries = ' + str(int(total_das)) + ' Tree produced = ' + str(int(total_processed)) + ' (frac = {0:.2f}'.format(fraction) + ')' + bcolors.ENDC
    else:
        print bcolors.BOLD + bcolors.FAIL + '   [NG] DAS entries = ' + str(int(total_das)) + ' Tree produced = ' + str(int(total_processed)) + ' (frac = {0:.2f}'.format(fraction) + ' < 0.8)' + bcolors.ENDC
    return True

def getSampleName(dasname):
  #if '__nanoaod' in dasname.lower():
  #  dasname = dasname[:dasname.lower().index('__nanoaod')]
  #if '__user' in dasname.lower():
  #  dasname = dasname[:dasname.lower().index('__user')]
  dasname = dasname.replace('__','/').lstrip('/')
  for (subdir,samplename), pattern in sample_dict.iteritems():
    if '$MASS' in pattern:
      matches = re.findall(pattern.replace('$MASS','(\d+)'),dasname)
      if matches:
        samplename = samplename.replace('$MASS',matches[0])
        return subdir, samplename
    elif pattern in dasname:
      return subdir, samplename
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
  
def matchSample(patterns,sample):
  if not isinstance(patterns,list):
    patterns = [patterns]
  for pattern in patterns:
    if '*' in pattern or '?' in pattern:
      if fnmatch(sample,pattern+'*'):
        return True
    else:
      if pattern in sample[:len(pattern)]:
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
    


if __name__ == '__main__':
    print
    main()



