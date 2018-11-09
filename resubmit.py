#! /usr/bin/env python

import os, glob, sys, shlex
from commands import getoutput
from argparse import ArgumentParser
from checkFiles import matchSample
import itertools
import subprocess
from ROOT import TFile, Double

parser = ArgumentParser()
parser.add_argument('-c', '--channel', dest='channel', type=str, default="tautau", action="store")
parser.add_argument('-n', '--njob',    dest='njob', type=int, default=10, action="store")
parser.add_argument('-m', '--make',    dest='make', default=False, action="store_true")
parser.add_argument('-s', '--samples', dest='samples', type=str, nargs='+', default=[ ], action="store",
                                       help="samples to run over, glob patterns (wildcards * and ?) are allowed.")
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
    
    for directory in sorted(os.listdir("./")):
        if not os.path.isdir(directory): continue
        if args.samples and not matchSample(args.samples,directory): continue
        #if directory.find('W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8__ytakahas-NanoTest_20180507_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8-a7a5b67d3e3590e4899e147be08660be__USER')==-1: continue
        filelist = glob.glob(directory + '/*_' + args.channel + '.root')
        if not filelist: continue
        total = 0
        
        ispnfs = False
        if directory.find('LQ')!=-1:
            ispnfs = True

        if ispnfs:
          #print 'Here!'
          files = getFileListPNFS(directory)
        else:
          #print 'No!'
          files = getFileListDAS('/' + directory.replace('__', '/'))
        
        filelists = list(split_seq(files, args.njob))
        
        jobList = 'joblist/joblist%s_%s_retry.txt'%(directory, args.channel)
        with open(jobList, 'w') as jobslog:
          ids = []
          #print 'ispnfs = ', ispnfs
          #print filelists
          
          for file2check in filelist:
              file = TFile(file2check, "read")
              if file.GetListOfKeys().Contains("tree"): continue
              total += 1
              id = file2check.split('_')[-2]
              #print 'id = ', id
              ids.append(id)
              nChunks = 0
              
              for filelist in filelists:
                if int(id) == nChunks:
                  createJobs(filelist,directory,directory,nChunks,args.channel,True,jobslog)
                else:
                  createJobs(filelist,directory,directory,nChunks,args.channel,False,jobslog)
                nChunks = nChunks+1
        
        if len(ids)==0:
            print bcolors.BOLD + bcolors.OKBLUE + '[OK] : ' + directory + bcolors.ENDC, ids
        else:
            print bcolors.BOLD + bcolors.FAIL + '[NG] : ' + directory + ', ' + str(len(ids)) + '/' + str(len(filelist)) + ' ... broken' + bcolors.ENDC, ids

        if args.make and len(ids)!=0:
            submitJobs(jobList, total, directory, batchSystem)


def split_seq(iterable, size):
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))

def getFileListDAS(dataset):

    instance = 'prod/global'
    if dataset.find('USER')!=-1:
        instance = 'prod/phys03'
    
    cmd='das_client --limit=0 --query="file dataset=%s instance=%s status=*"'%(dataset,instance)
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = []
    for l in tmpList:
        if l.find(".root") != -1:
            files.append(l)
	         
    return files


def getFileListPNFS(dataset):
#    instance = 'prod/global'
#    if dataset.find('USER')!=-1:
#        instance = 'prod/phys03'
    
#    cmd='das_client --limit=0 --query="file dataset=%s instance=%s"'%(dataset,instance)
    name = '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/' + dataset.replace('__','/')
    cmd='ls %s'%(name)
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = []
    for l in tmpList:
        if l.find(".root") != -1:
            files.append(name + '/' + l.rstrip())
    
#nanoAOD_LQ3ToTauB_Fall2017_5f_Madgraph_LO_pair-M2000_646.root	         
#VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_pair_M1000__nanoAOD__v1/
#/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_s_channel_M500/nanoAOD/v1/nanoAOD_VectorLQ3ToTauB_Fall2017_5f_Madgraph_LO_s_channel_M500_1.root
    
    return files 


def createJobs(filelist,outfolder,name,nchunks,channel,toWrite,jobslog):
  infiles = []
  for files in filelist:
    #infiles.append("root://cms-xrd-global.cern.ch/"+files)
    infiles.append("dcap://t3se01.psi.ch:22125/"+files)
  cmd = 'python job.py %s %s %s %i %s \n'%(','.join(infiles), outfolder,name,nchunks, channel)
  #print cmd
  
  if toWrite:
    jobslog.write(cmd)
  
  return 1


def submitJobs(jobList, nchunks, outfolder, batchSystem):
    print 'Reading joblist'
    jobListName = jobList
    print jobList
#    subCmd = 'qsub -t 1-%s -o logs nafbatch_runner_GEN.sh %s' %(nchunks,jobListName)
    subCmd = 'qsub -t 1-%s -o %s/logs/ %s %s' %(nchunks,outfolder,batchSystem,jobListName)
    print 'Going to submit', nchunks, 'jobs with', subCmd
    os.system(subCmd)

    return 1

            
#    if flag:
#        print bcolors.FAIL + "[NG]" + directory + bcolors.ENDC
#        print '\t', len(files), ' out of ', str(total) + ' files are corrupted ... skip this sample (consider to resubmit the job)'
#
#    else:
#        print bcolors.BOLD + bcolors.OKBLUE + '[OK] ' + directory + bcolors.ENDC


if __name__ == '__main__':
    print
    main()
