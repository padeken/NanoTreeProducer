#!/usr/bin/env python

import os, re, glob
from commands import getoutput
from fnmatch import fnmatch
import itertools
from argparse import ArgumentParser
from checkFiles import getSampleName

parser = ArgumentParser()
parser.add_argument('-f', '--force',   dest='force',   action="store_true", default=False,
                                       help="do not ask for confirmation before submission of jobs" )
parser.add_argument('-c', '--channel', dest='channel', action="store", type=str, default="mutau",
                                       help="channels to submit" )
parser.add_argument('-s', '--sample',  dest='sample',  action="store", type=str, default=None,
                                       help="filter this sample" )
parser.add_argument('-y', '--year',    dest='year', choices=[2017,2018], type=int, default=2017, action='store',
                                       help="select year")
parser.add_argument('-n', '--njob',    dest='njob',    action="store", type=int, default=4,
                                       help="number of files per job" )
parser.add_argument('-m', '--mock',    dest='mock',    action="store_true", default=False,
                                       help="mock submit jobs for debugging purposes" )
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
    if dataset.find('USER')!=-1:
        instance = 'prod/phys03'
    
    #cmd='das_client --limit=0 --query="file dataset=%s instance=%s"'%(dataset,instance)
    cmd='das_client --limit=0 --query="file dataset=%s instance=%s status=*"'%(dataset,instance)
    print "Executing ",cmd
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = []
    for l in tmpList:
        if l.find(".root") != -1:
            files.append(l)
    
    return files 


def getFileListPNFS(dataset):
    
    #instance = 'prod/global'
    #if dataset.find('USER')!=-1:
    #    instance = 'prod/phys03'
    #cmd='das_client --limit=0 --query="file dataset=%s instance=%s"'%(dataset,instance)
    
    user = 'ytakahas'
    name = '/pnfs/psi.ch/cms/trivcat/store/user/'+user+'/' + dataset.replace('__','/')
    cmd='ls %s'%(name)
    print "Executing ",cmd
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = []
    for l in tmpList:
        if l.find(".root") != -1:
            files.append(name + '/' + l.rstrip())
    
    return files 

   
def createJobs(jobs, filelist, outdir, name, nchunks, channel, pattern):
  infiles = [ ]
  for files in filelist:
      #if pattern.find('pnfs')!=-1:
      #    infiles.append("dcap://t3se01.psi.ch:22125/"+ pattern + '/' + files)
      #    infiles.append("root://cms-xrd-global.cern.ch/"+ pattern.replace('/pnfs/psi.ch/cms/trivcat','') + '/' + files)
      #else:
      if files.find('LQ')!=-1:
          infiles.append("dcap://t3se01.psi.ch:22125/"+files)
      else:
          infiles.append("root://cms-xrd-global.cern.ch/"+files)
  cmd = 'python job.py %s %s %s %i %s \n'%(','.join(infiles), outdir,name,nchunks, channel)
  #print cmd
  jobs.write(cmd)
  return 1


def submitJobs(jobName, jobList, nchunks, outdir, batchSystem):
    print 'Reading joblist...'
    jobListName = jobList
    print jobList
    #subCmd = 'qsub -t 1-%s -o logs nafbatch_runner_GEN.sh %s' %(nchunks,jobListName)
    subCmd = 'qsub -t 1-%s -N %s -o %s/logs/ %s %s'%(nchunks,jobName,outdir,batchSystem,jobListName)
    print bcolors.BOLD + bcolors.OKBLUE + "Submitting %d jobs with \n    %s"%(nchunks,subCmd) + bcolors.ENDC
    if not args.mock:
      os.system(subCmd)
    return 1
    

def main():
    
    batchSystem = 'psibatch_runner.sh'
    channel = args.channel
    samplelist = "samples_%s.cfg"%(args.year)
    
    # READ SAMPLES
    patterns = [ ]
    for line in open(samplelist, 'r'):
        if line.find('#')!=-1: continue
        line = line.rstrip()
        if line=='': continue
        #if line.count('/')!=3:
        #    continue
        patterns.append(line)
    #print patterns
	
	# SUBMIT SAMPLES
    for pattern in patterns:
        
        ispnfs = False
        if pattern.find('pnfs')!=-1:
            ispnfs = True
        
        if args.sample!=None:
            if '*' in args.sample or '?' in args.sample:
              if not fnmatch(pattern,'*'+args.sample+'*'): continue
            elif pattern.find(args.sample)==-1: continue
        
        if channel=='tautau':
            if pattern.find('/SingleMuon')!=-1 or pattern.find('/SingleElectron')!=-1: continue

        if channel in ['mutau', 'mumu', 'muele']:
            if pattern.find('/SingleElectron')!=-1 or pattern.find('/Tau')!=-1: continue
        
        if channel=='eletau':
            if pattern.find('/SingleMuon')!=-1 or pattern.find('/Tau')!=-1: continue
        
        print bcolors.BOLD + bcolors.OKGREEN + pattern + bcolors.ENDC
        files = None
        name = None
        
        if ispnfs:
            name = pattern.split('/')[8].replace('/','') + '__' + pattern.split('/')[9].replace('/','') + '__' + pattern.split('/')[10].replace('/','')
            #files = getFileListPNFS(pattern)
            files = getFileListPNFS(name)
        else:
            files = getFileListDAS(pattern)
            name = pattern.split('/')[1].replace('/','') + '__' + pattern.split('/')[2].replace('/','') + '__' + pattern.split('/')[3].replace('/','')
        
        #if files:
        #  print "FILELIST = "+files[0]
        #  for file in files[1:]:
        #    print "           "+file
        #else:
        if not files:
          print bcolors.BOLD + bcolors.WARNING + "Warning!!! FILELIST empty" + bcolors.ENDC
          continue
        
        # JOBLIST
        jobList = 'joblist/joblist%s_%s.txt'%(name,channel)
        print "Creating job file %s..."%(jobList)
        try: os.stat('joblist/')
        except: os.mkdir('joblist/')
        jobName = getSampleName(pattern)[1]
        jobs = open(jobList, 'w')
        njob = args.njob
        nChunks = 0
        outdir = name
        
        # NJOBS CHECKS
        if njob>1 and any(v==jobName for v in [ 'WW', 'WZ', 'ZZ' ]):
          print bcolors.BOLD + bcolors.WARNING + "Warning: Setting number of files per job from %s to 1 for %s"%(njob,jobName) + bcolors.ENDC
          njob = 1
        
        try: os.stat(outdir)
        except: os.mkdir(outdir)
        try: os.stat(outdir+'/logs/')
        except: os.mkdir(outdir+'/logs/')
        
        # CREATE JOBS
        filelists = list(split_seq(files,njob))
        checkExistingFiles(outdir,channel,len(filelists))
        #filelists = list(split_seq(files,1))
        for file in filelists:
        #print "FILES = ",f
            createJobs(jobs, file,outdir,name,nChunks,channel,pattern)
            nChunks = nChunks+1
        jobs.close()
        
        # SUBMIT
        if args.force:
          submitJobs(jobName,jobList,nChunks,outdir,batchSystem)
        else:
          submit = raw_input("Do you also want to submit %d jobs to the batch system? [y/n] "%(nChunks))
          if submit.lower()=='force':
            submit = 'y'
            args.force = True
          if submit.lower()=='y':
            submitJobs(jobName,jobList,nChunks,outdir,batchSystem)
          else:
            print "Not submitting jobs"
        print



if __name__ == "__main__":
    print
    main()
    print "Done\n"
		
		
		
