#! /usr/bin/env python

import os, glob, sys, shlex, re
#import time
from datetime import datetime
from fnmatch import fnmatch
import subprocess
from argparse import ArgumentParser
from checkFiles import getSampleShortName, matchSampleToPattern, header

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
    parser.add_argument('-y', '--year',    dest='years', choices=[2016,2017,2018], type=int, nargs='+', default=[2017], action='store',
                                           help="select year" )
    parser.add_argument('-c', '--channel', dest='channels', choices=['mutau','eletau','tautau'], nargs='+', default=['tautau'], action='store' )
    parser.add_argument('-n', '--njobs',   dest='njobs', type=int, default=1, action='store',
                                           help="number of last jobs to check" )
    parser.add_argument('-o', '--outdir',  dest='outdir', type=str, default=None, action='store' )
    parser.add_argument('-s', '--sample',  dest='samples', type=str, nargs='+', default=[ ], action='store',
                                           help="samples to run over, glob patterns (wildcards * and ?) are allowed." )
    parser.add_argument('-x', '--veto',    dest='veto', action='store', type=str, default=None,
                                           help="veto this sample" )
    parser.add_argument('-t', '--type',    dest='type', choices=['data','mc'], type=str, default=None, action='store',
                                           help="filter data or MC to submit" )
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                                           help="set verbose" )
    args = parser.parse_args()

filepattern     = re.compile(r".*\.o(\d+)\.(\d+)")
qstatpattern    = re.compile(r"(\d+) *\d\.\d+ *\w+ *\w+ *(\w+) *\d\d/\d\d/20\d\d *\d\d:\d\d:\d\d *(?:[^ ]+@[^ ]+)? *\d+ *(\d+)")
nodepattern     = re.compile(r"Running job on machine .* (t3wn\d+\.psi\.ch)")
queuepattern    = re.compile(r"(\w+\.q)@(t3wn\d+\.psi\.ch)")
rootpattern     = re.compile(r"root://.*?/.*?\.root")
RTWpattern      = "RuntimeWarning: creating executor for unknown type"
jobstartpattern = re.compile(r"job start at (\w+ \w+ \d+ \d\d:\d\d:\d\d \w+ 20\d\d)")
ppstartpattern  = re.compile(r"Pre-select \d+ entries out of \d+")
donepattern     = re.compile(r"Complete at *(\w+ \w+ \d+ \d\d:\d\d:\d\d \w+ 20\d\d)")
chunkpattern    = re.compile(r"Going to execute python job.py .* -n (\d+) .*")



class Job:
  
  def __init__(self,logfile):
    nevents   = -1
    jobid     = -1
    taskid    = -1
    runtime   = -1
    chunk     = -1
    jobstart  = None
    done      = None
    queue     = None
    node      = None
    ppstart   = False # post-processor started running
    running   = False
    stuck     = False
    failed    = False
    files     = [ ]
    
    # JOB INFO
    jobid, taskid = getJobID(logfile)
    
    # QSTAT
    process   = subprocess.Popen("qstat", stdout=subprocess.PIPE, shell=True)
    (out,err) = process.communicate()
    status    = process.wait()
    for line in out.split('\n'):
      match = qstatpattern.match(line)
      if match and jobid==int(match.group(1)) and taskid==int(match.group(3)):
        running = match.group(2)=='r'
        failed  = 'E' in match.group(2)
        match   = queuepattern.search(line)
        if match:
          queue = match.group(1)
          node  = match.group(2)
        break
    
    # READ FILE
    with open(logfile) as file:
      for line in file:
        
        match = jobstartpattern.search(line)
        if match:
          jobstart = datetime.strptime(match.group(1),'%a %b %d %X CET %Y')
          continue
        
        match = chunkpattern.search(line)
        if match:
          chunk = int(match.group(1))
          continue
        
        match = ppstartpattern.search(line)
        if match:
          ppstart = True
          continue
        
        match = donepattern.search(line)
        if match:
          done = datetime.strptime(match.group(1),'%a %b %d %X CET %Y')
          break
        
        match = nodepattern.search(line)
        if match:
          node = match.group(1)
          continue
        
        matches = rootpattern.findall(line)
        if matches and not files:
          files = matches
          continue
        
        if "segmentation violation" in line.lower() or "file probably overwritten: stopping reporting error messages" in line:
          #or "terminate called after throwing an instance of 'std::bad_alloc'" in line:
          failed = True
          done = False
          break
    
    if jobstart:
      if done:
        runtime = (done - jobstart).seconds
      else:
        if running:
          runtime = (datetime.now() - jobstart).seconds
          if not ppstart:
            stuck = runtime > 60*20 # 20 min.
        else:
          stuck = not failed
    
    self.logfile  = logfile
    self.files    = files
    self.chunk    = chunk
    self.runtime  = runtime
    self.start    = jobstart
    self.nevents  = nevents
    self.jobid    = jobid
    self.taskid   = taskid
    self.node     = node
    self.queue    = queue
    self.running  = running
    self.stuck    = stuck
    self.done     = done
    self.failed   = failed
    
  def __gt__(self, ojob):
    if self.jobid==ojob.jobid:
      #return self.taskid>ojob.taskid
      return self.chunk>ojob.chunk
    return self.jobid>ojob.jobid
  
  def __str__(self):
    return '%d.%d (%d)'%(self.jobid,self.taskid,self.chunk)
    


def getJobID(filename):
  match = filepattern.match(filename)
  if match:
    return int(match.group(1)), int(match.group(2))
  else:
    print ">>> Warning! Job::__init__: %s does not match to a job log file with a job and task ID!"
  return -1, -1
  


def printTime(seconds):
  hours, remainder = divmod(seconds, 3600)
  minutes, seconds = divmod(remainder, 60)
  string = ""
  if hours>1:     string  = "%d hour "%(hours)
  elif hours>0:   string  = "%d hours "%(hours)
  if minutes>1:   string += "%d minute "%(minutes)
  elif minutes>0: string += "%d minutes "%(minutes)
  if seconds>1:   string += "%d second"%(seconds)
  elif minutes>0: string += "%d seconds"%(seconds)
  return string
  


def main(args):
  
  years      = args.years
  channels   = args.channels
  njobs      = args.njobs
  
  for year in years:
    indir    = "output_%s/"%(year)
    os.chdir(indir)
    
    # GET LIST
    samplelist = [ ]
    for directory in sorted(os.listdir('./')):
        if not os.path.isdir(directory): continue
        if args.samples and not matchSampleToPattern(directory,args.samples): continue
        if args.veto and matchSampleToPattern(directory,args.veto): continue
        if args.type=='mc' and any(s in directory[:len(s)+2] for s in ['SingleMuon','SingleElectron','Tau']): continue
        if args.type=='data' and not any(s in directory[:len(s)+2] for s in ['SingleMuon','SingleElectron','Tau']): continue
        samplelist.append(directory)
    if not samplelist:
      print "No samples found in %s!"%(indir)
    if args.verbose:
      print 'samplelist = %s\n'%(samplelist)
    
    # CHECK samples
    for channel in channels:
      print header(year,channel)
      
      for directory in samplelist:
        print ">>> %s"%(directory)
        
        infiles  = "%s/logs/*%s_%d*.o*.*"%(directory,channel,year)
        filelist = glob.glob(infiles)
        if not filelist:
          continue
        
        jobids   = [ ]
        for filename in filelist:
          jobid, taskid = getJobID(filename)
          if jobid not in jobids:
            jobids.append(jobid)
        jobids.sort(reverse=True)
        jobids_max = jobids[:njobs]
        
        jobs    = { id: [ ] for id in jobids_max}
        stuck   = { id: [ ] for id in jobids_max}
        failed  = { id: [ ] for id in jobids_max}
        running = { id: [ ] for id in jobids_max}
        done    = { id: [ ] for id in jobids_max}
        for filename in filelist:
          if not any(".o%d."%(id) in filename for id in jobids_max):
            continue
          job = Job(filename)
          jobs[job.jobid].append(job)
          if job.stuck:
            stuck[job.jobid].append(job)
          if job.running:
            running[job.jobid].append(job)
          if job.failed:
            failed[job.jobid].append(job)
          if job.done:
            done[job.jobid].append(job)
        
        for jobid, joblist in sorted(jobs.iteritems()):
          ntot = len(joblist)
          jobs[jobid].sort()
          stuck[jobid].sort()
          failed[jobid].sort()
          running[jobid].sort()
          done[jobid].sort()
          print ">>>   %d"%(jobid)
          if running[jobid]:
            print ">>>     running: %4d /%5d "%(len(running[jobid]),ntot)
          if failed[jobid]:
            print ">>>     failed:  %4d /%5d "%(len(failed[jobid]),ntot) #+ ', '.join([str(j) for j in failed[jobid]])
          if stuck[jobid]:
            print ">>>     stuck:   %4d /%5d "%(len(stuck[jobid]),ntot)
          print ">>>     done:    %4d /%5d "%(len(done[jobid]),ntot) #+ ', '.join([str(j) for j in done[jobid]])
        
        print ">>>"



if __name__ == '__main__':    
  print
  main(args)
  print
  

