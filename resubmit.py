#! /usr/bin/env python

import os, glob, sys, shlex
from commands import getoutput
from argparse import ArgumentParser
from checkFiles import matchSample
from submit_qsub import createJobs, getFileListPNFS, getFileListDAS, submitJobs
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


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def main():
    
    batchSystem = 'psibatch_runner.sh'    
    year        = args.year
    channel     = args.channel
    outdir      = "output_%s/"%(year)
    os.chdir(outdir)
    
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
        filelist = glob.glob(directory + '/*_' + args.channel + '.root')
        if not filelist: continue
        
        # GET INPUT FILES
        if 'LQ' in directory:
          files = getFileListPNFS(directory)
        else:
          files = getFileListDAS('/' + directory.replace('__', '/'))
        
        total = 0
        outdir = "output_%s/%s"%(year,directory)
        filelists = list(split_seq(files, args.njob))
        
        jobList = 'joblist/joblist%s_%s_retry.txt'%(directory, args.channel)
        with open(jobList, 'w') as jobslog:
          ids = []
          #print filelists
          
          # TODO: check if files are missing
          for file2check in filelist:
              file = TFile(file2check,'READ')
              if not file.IsZombie() and file.GetListOfKeys().Contains('tree') and file.GetListOfKeys().Contains('cutflow'):
                continue
              total += 1
              id = file2check.split('_')[-2]
              #print 'id = ', id
              ids.append(id)
              nChunks = 0
              
              for filelist in filelists:
                if int(id) == nChunks:
                  createJobs(jobslog,filelist,outdir,directory,nChunks,channel,year=year,write=True)
                else:
                  createJobs(jobslog,filelist,outdir,directory,nChunks,channel,year=year,write=False)
                nChunks = nChunks+1
                
        if len(ids)==0:
            print bcolors.BOLD + bcolors.OKBLUE + '[OK] : ' + directory + bcolors.ENDC, ids
        else:
            print bcolors.BOLD + bcolors.FAIL + '[NG] : ' + directory + ', ' + str(len(ids)) + '/' + str(len(filelist)) + ' ... broken' + bcolors.ENDC, ids

        if len(ids)!=0:
            submitJobs(jobList, total, directory, batchSystem)
    

def split_seq(iterable, size):
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size)) 
    
#    if flag:
#        print bcolors.FAIL + "[NG]" + directory + bcolors.ENDC
#        print '\t', len(files), ' out of ', str(total) + ' files are corrupted ... skip this sample (consider to resubmit the job)'
#
#    else:
#        print bcolors.BOLD + bcolors.OKBLUE + '[OK] ' + directory + bcolors.ENDC


if __name__ == '__main__':
    print
    main()
