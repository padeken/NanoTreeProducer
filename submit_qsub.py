#!/usr/bin/env python

import os, re, glob
from commands import getoutput
from fnmatch import fnmatch
import itertools
from argparse import ArgumentParser
import checkFiles
from checkFiles import getSampleShortName, matchSampleToPattern

parser = ArgumentParser()
parser.add_argument('-f', '--force',   dest='force', action='store_true', default=False,
                                       help="do not ask for confirmation before submission of jobs" )
parser.add_argument('-c', '--channel', dest='channel', action='store', type=str, default="mutau",
                                       help="channels to submit" )
parser.add_argument('-s', '--sample',  dest='samples', type=str, nargs='+', default=[ ], action='store',
                                       help="filter these samples, glob patterns (wildcards * and ?) are allowed." )
parser.add_argument('-x', '--veto',    dest='veto', action='store', type=str, default=None,
                                       help="veto this sample" )
parser.add_argument('-y', '--year',    dest='year', choices=[2017,2018], type=int, default=2017, action='store',
                                       help="select year" )
parser.add_argument('-T', '--tes',     dest='tes', type=float, default=1.0, action='store',
                                       help="tau energy scale" )
parser.add_argument('-n', '--njob',    dest='nFilesPerJob', action='store', type=int, default=4,
                                       help="number of files per job" )
parser.add_argument('-m', '--mock',    dest='mock', action='store_true', default=False,
                                       help="mock submit jobs for debugging purposes" )
parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                                       help="set verbose" )
args = parser.parse_args()
checkFiles.args = args

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    

def checkExistingFiles(outdir,channel,njob):
    filelist = glob.glob("%s/*%s.root"%(outdir,channel))
    nfiles = len(filelist)
    if nfiles>njob:
      print bcolors.BOLD + bcolors.WARNING + "Warning! There already exist %d files, while the requested number of files per job is %d"%(nfiles,njob) + bcolors.ENDC
      remove = raw_input("Do you want to remove the extra files? [y/n] ")
      if remove.lower()=='y':
        for filename in sorted(filelist):
          matches = re.findall(r"_(\d+)_%s.root"%(channel),filename)
          if matches and int(matches[0])>njob:
            print "Removing %s..."%(filename)
            os.remove(filename)
      else:
        print "Not removing extra files. Please make sure to delete the %d last files before hadd'ing."%(nfiles-njob)
    

def split_seq(iterable, size):
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))
    

def getFileListDAS(dataset):
    instance = 'prod/global'
    if 'USER' in dataset:
        instance = 'prod/phys03'
    #cmd='das_client --limit=0 --query="file dataset=%s instance=%s"'%(dataset,instance)
    cmd = 'das_client --limit=0 --query="file dataset=%s instance=%s status=*"'%(dataset,instance)
    if args.verbose:
      print "Executing ",cmd
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = [ ]
    for line in tmpList:
        if '.root' in line:
            files.append(line)    
    return files 
    

def getFileListPNFS(dataset):
    #instance = 'prod/global'
    #if dataset.find('USER')!=-1:
    #    instance = 'prod/phys03'
    #cmd='das_client --limit=0 --query="file dataset=%s instance=%s"'%(dataset,instance)
    user = 'ytakahas'
    name = '/pnfs/psi.ch/cms/trivcat/store/user/'+user+'/' + dataset.replace('__','/')
    cmd='ls %s'%(name)
    if args.verbose:
      print "Executing ",cmd
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = []
    for line in tmpList:
        if '.root' in line:
            files.append(name+'/'+line.rstrip())
    return files


def createJobs(jobsfile, filelist, outdir, name, nchunks, channel, year=2017):
  infiles = [ ]
  for file in filelist:
      #if pattern.find('pnfs')!=-1:
      #    infiles.append("dcap://t3se01.psi.ch:22125/"+ pattern + '/' + file)
      #    infiles.append("root://cms-xrd-global.cern.ch/"+ pattern.replace('/pnfs/psi.ch/cms/trivcat','') + '/' + file)
      #else:
      if 'LQ' in file:
          infiles.append("dcap://t3se01.psi.ch:22125/"+file)
      else:
          infiles.append("root://cms-xrd-global.cern.ch/"+file)
  cmd = 'python job.py -i %s -o %s -N %s -n %i -c %s -y %s \n'%(','.join(infiles),outdir,name,nchunks,channel,args.year)
  if args.verbose:
    print cmd
  jobsfile.write(cmd)
  return 1
  

def submitJobs(jobName, jobList, nchunks, outdir, batchSystem):
    if args.verbose:
      print 'Reading joblist...'
      print jobList
    #subCmd = 'qsub -t 1-%s -o logs nafbatch_runner_GEN.sh %s' %(nchunks,jobList)
    subCmd = 'qsub -t 1-%s -N %s -o %s/logs/ %s %s'%(nchunks,jobName,outdir,batchSystem,jobList)
    print bcolors.BOLD + bcolors.OKBLUE + "Submitting %d jobs with \n    %s"%(nchunks,subCmd) + bcolors.ENDC
    if not args.mock:
      os.system(subCmd)
    return 1
    

def main():
    
    batchSystem = 'psibatch_runner.sh'
    channel     = args.channel
    year        = args.year
    tes         = args.tes
    samplelist  = "samples_%s.cfg"%(year)
    
    # READ SAMPLES
    directories = [ ]
    for line in open(samplelist, 'r'):
        line = line.rstrip().lstrip().split(' ')[0]
        if line[:2].count('#')>0: continue
        if line=='': continue
        if args.samples and not matchSampleToPattern(line,args.samples): continue
        if args.veto and matchSampleToPattern(line,args.veto): continue
        #if line.count('/')!=3:
        #    continue
        directories.append(line)
    #print directories
	
	# SUBMIT SAMPLES
    for directory in directories:
        
        if args.verbose:
          print "\ndirectory =",directory
        
        # FILTER
        if 'SingleMuon' in directory and channel not in ['mutau','mumu']: continue
        if 'SingleElectron' in directory and channel!='etau': continue
        if 'Tau' in directory and channel!='tautau': continue
        
        print bcolors.BOLD + bcolors.OKGREEN + directory + bcolors.ENDC
        files = None
        name = None
        
        if 'pnfs' in directory:
            name = directory.split('/')[8].replace('/','') + '__' + directory.split('/')[9].replace('/','') + '__' + directory.split('/')[10].replace('/','')
            #files = getFileListPNFS(directory)
            files = getFileListPNFS(name)
        else:
            files = getFileListDAS(directory)
            name = directory.split('/')[1].replace('/','') + '__' + directory.split('/')[2].replace('/','') + '__' + directory.split('/')[3].replace('/','')
        
        if not files:
          print bcolors.BOLD + bcolors.WARNING + "Warning!!! FILELIST empty" + bcolors.ENDC
          continue
        elif args.verbose:
          print "FILELIST = "+files[0]
          for file in files[1:]:
            print "           "+file
        
        # JOBLIST
        jobList = 'joblist/joblist%s_%s.txt'%(name,channel)
        print "Creating job file %s..."%(jobList)
        try: os.stat('joblist/')
        except: os.mkdir('joblist/')
        jobName = getSampleShortName(directory)[1]
        jobs    = open(jobList, 'w')
        nFilesPerJob = args.nFilesPerJob
        outdir  = "output_%s/%s"%(year,name)
        
        # NFILESPERJOBS CHECKS
        # Diboson (WW, WZ, ZZ) have very large files and acceptance,
        # and the jet-binned DY and WJ files need to be run separately because of a bug affecting LHE_Njets
        if nFilesPerJob>1 and any(vv in jobName[:4] for vv in [ 'WW', 'WZ', 'ZZ', 'DY', 'WJ', 'W1J', 'W2J', 'W3J', 'W4J' ]):
          print bcolors.BOLD + bcolors.WARNING + "[WN] setting number of files per job from %s to 1 for %s"%(nFilesPerJob,jobName) + bcolors.ENDC
          nFilesPerJob = 1
        
        try: os.stat(outdir)
        except: os.mkdir(outdir)
        try: os.stat(outdir+'/logs/')
        except: os.mkdir(outdir+'/logs/')
        
        # CREATE JOBS
        nChunks = 0
        filelists = list(split_seq(files,nFilesPerJob))
        checkExistingFiles(outdir,channel,len(filelists))
        #filelists = list(split_seq(files,1))
        for file in filelists:
        #print "FILES = ",f
            createJobs(jobs,file,outdir,name,nChunks,channel,year=year)
            nChunks = nChunks+1
        jobs.close()
        
        # SUBMIT
        jobName += "_%s_%s"%(channel,year)
        if args.force:
          submitJobs(jobName,jobList,nChunks,outdir,batchSystem)
        else:
          submit = raw_input("Do you also want to submit %d jobs to the batch system? [y/n] "%(nChunks))
          if submit.lower()=='force':
            submit = 'y'
            args.force = True
          if submit.lower()=='quit':
            exit(0)
          if submit.lower()=='y':
            submitJobs(jobName,jobList,nChunks,outdir,batchSystem)
          else:
            print "Not submitting jobs"
        print



if __name__ == "__main__":
    
    print
    main()
    print "Done\n"
		
		
		
