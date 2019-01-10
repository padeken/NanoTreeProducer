#! /usr/bin/env python

import os, glob, sys, shlex
from commands import getoutput
from argparse import ArgumentParser
from checkFiles import matchSample
from submit_qsub import bcolors, createJobs, getFileListPNFS, getFileListDAS, submitJobs, split_seq
import itertools
import subprocess
from ROOT import TFile, Double

parser = ArgumentParser()
parser.add_argument('-c', '--channel', dest='channel', type=str, default="tautau", action="store")
parser.add_argument('-n', '--njob',    dest='njob', type=int, default=10, action="store")
#parser.add_argument('-m', '--make',    dest='make', default=False, action="store_true")
parser.add_argument('-m', '--mock',    dest='mock', action='store_true', default=False,
                                       help="mock submit jobs for debugging purposes" )
parser.add_argument('-y', '--year',    dest='year', choices=[2017,2018], type=int, default=2017, action='store',
                                       help="select year" )
parser.add_argument('-s', '--samples', dest='samples', type=str, nargs='+', default=[ ], action="store",
                                       help="samples to run over, glob patterns (wildcards * and ?) are allowed.")
parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                                       help="set verbose" )
args = parser.parse_args()



def main():
    
    batchSystem = 'psibatch_runner.sh'    
    year        = args.year
    channel     = args.channel
    outdir      = "output_%s/"%(year)
    os.chdir(outdir)
    chunkpattern = re.compile(r".*_(\d+)_[a-z]+\.root")
    
    # GET LIST
    samplelist = [ ]
    for directory in sorted(os.listdir('./')):
        if not os.path.isdir(directory): continue
        if args.samples and not matchSampleToPattern(args.samples,directory): continue
        samplelist.append(directory)
    if not samplelist:
      print "No samples found in %s!"%(outdir)
    if args.verbose:
      print samplelist
    
    # RESUBMIT samples
    for directory in samplelist:
        #if directory.find('W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8__ytakahas-NanoTest_20180507_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be__USER')==-1: continue
        outfilelist = glob.glob(directory + '/*_' + args.channel + '.root')
        if not outfilelist: continue
        
        # GET INPUT FILES
        if 'LQ' in directory:
          infiles = getFileListPNFS(directory)
        else:
          infiles = getFileListDAS('/' + directory.replace('__', '/'))
        
        outdir = "output_%s/%s"%(year,directory)
        infilelists = list(split_seq(infiles, args.njob))
        
        foundchunks = [ ]
        badchunks   = [ ]
        misschunks  = [ ]
        jobList = 'joblist/joblist%s_%s_retry.txt'%(directory, args.channel)
        with open(jobList, 'w') as jobslog:
          for filename in outfilelist:
              match = chunkpattern.search(filename)
              if match:
                chunk = match.group(0)
              else:
                print bcolors.BOLD + bcolors.FAIL + '[NG] did not recognize output file %s !'%(directory) + bcolors.ENDC
                exit(1)
              #chunk = filename.split('_')[-2]
              foundchunks.append(chunk)
              file = TFile(filename,'READ')
              if not file.IsZombie() and file.GetListOfKeys().Contains('tree') and file.GetListOfKeys().Contains('cutflow'):
                continue
              infiles = infilelists[chunk]
              createJobs(jobslog,infiles,outdir,directory,chunk,channel,year=year)
              badchunks.append(chunk)
          
          # BAD CHUNKS
          if len(badchunks)>0:
            chunktext = ('chunks' if len(badchunks)>1 else 'chunk') + ', '.join(str(ch) for ch in badchunks)
            print bcolors.BOLD + bcolors.WARNING + '[NG] %s, %d/%d failed ! Resubmitting %s...'%(directory,len(ids),len(outfilelist),chunktext) + bcolors.ENDC
          
          # MISSING CHUNKS
          maxchunk = max(foundchunks)+1
          if len(outfilelist)<maxchunk:
            misschunks = [ i for i in range(0,max(foundchunks)) if i not in foundchunks ]
            chunktext = ('chunks' if len(misschunks)>1 else 'chunk') + ', '.join(str(i) for i in misschunks)
            print bcolors.BOLD + bcolors.WARNING + "[WN] %s missing %d/%d files ! Resubmitting %s...%s"%(directory,imax-len(outfilelist),len(outfilelist),chunktext) + bcolors.ENDC
            for chunk in misschunks:
              infiles = infilelists[chunk]
              createJobs(jobslog,infiles,outdir,directory,chunk,channel,year=year)
        
        # RESUBMIT
        nchunks = len(badchunks)+len(misschunks)
        if nchunks>0:
            submitJobs(jobList,nchunks,directory,batchSystem)
        else:
            print bcolors.BOLD + bcolors.OKBLUE + '[OK] ' + directory + bcolors.ENDC
        
#    if flag:
#        print bcolors.FAIL + "[NG]" + directory + bcolors.ENDC
#        print '\t', len(files), ' out of ', str(total) + ' files are corrupted ... skip this sample (consider to resubmit the job)'
#
#    else:
#        print bcolors.BOLD + bcolors.OKBLUE + '[OK] ' + directory + bcolors.ENDC


if __name__ == '__main__':
    print
    main()
