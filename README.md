# NanoTreeProducer
Produce skimmed analysis trees directly from centrally-produced NanoAOD.
This repository is meant for an analysis with a pair of tau leptons in several final states, and for 2016, 2017 and 2018 data.


## Installation

First, install NanoAODTools:
```
cmsrel CMSSW_9_4_6
cd CMSSW_9_4_6/src
cmsenv
git cms-init   #not really needed unless you later want to add some other cmssw stuff
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b
```

Then, install this package:
```
git clone https://github.com/gitytakahas/NanoTreeProducer
```


## Analysis

You want to change the analyse code in the modules and tree branches for whatever analysis you want to perform.
The examples below are for an analysis selecting for a muon and tau:
```
MuTauModule.py
TreeProducerMuTau.py
TreeProducerCommon.py
```


## Run

For a local run, do something like
```
./local.py -c mutau -y 2017
```

For job submission, you need to modify the list of samples you want to process in the config file, e.g.
```
samples_2017.cfg
```
and then do, something like
```
./submit_qsub.py -c mutau -y 2017
```

To check that the jobs were succesful, you need to ensure that all the output file contains the expected tree with the expected number of events (`-d`):
```
./submit_qsub.py -c mutau -y 2017 -d
```
If the output is fine, one can hadd (`-m`) all the output:
```
./checkFiles.py -c mutau -y 2017 -m
```
To resubmit failed jobs, do:
```
./resubmit.py -c mutau -y 2017
```
