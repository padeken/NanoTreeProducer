#! /usr/bin/env python
# Author: Izaak Neutelings (January 2019)
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2017Analysis
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2018Analysis
# https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#PileupInformation

import os, sys
from argparse import ArgumentParser
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TFile, TTree

argv = sys.argv
description = '''This script make some checks.'''
parser = ArgumentParser(prog="pileup",description=description,epilog="Succes!")
parser.add_argument('-y', '--year',     dest='years', choices=[2016,2017,2018], type=int, nargs='+', default=[2018], action='store',
                                        help="select year" )
parser.add_argument('-c', '--channel',  dest='channel', choices=['mutau','etau'], type=str, default='mutau', action='store',
                                        help="select channel" )
parser.add_argument('-v', '--verbose',  dest='verbose', default=False, action='store_true', 
                                        help="print verbose" )
args = parser.parse_args()



def getMCProfile(outfilename,indir,samples,channel,year):
    """Get pileup profile in MC by adding Pileup_nTrueInt histograms from a given list of samples."""
    print ">>> getMCProfile(%s)"%(outfilename)
    nprofiles = 0
    histname  = 'pileup'
    tothist   = None    
    for subdir, samplename in samples:
      filename = "%s/%s/%s_%s.root"%(indir,subdir,samplename,channel)
      print ">>>   %s"%(filename)
      file = TFile(filename,'READ')
      if not file or file.IsZombie():
        print ">>>   Warning! getMCProfile: Could not open %s"%(filename)
        continue
      hist      = file.Get(histname)
      if not hist:
        print ">>>   Warning! getMCProfile: Could not open histogram in %s"%(filename)      
        continue
      if tothist==None:
        tothist = hist.Clone('pileup')
        tothist.SetTitle('pileup')
        tothist.SetDirectory(0)
        nprofiles += 1
      else:
        tothist.Add(hist)
        nprofiles += 1
      file.Close()
    print ">>>   added %d MC profiles, %d entries, %.1f mean"%(nprofiles,tothist.GetEntries(),tothist.GetMean())
    
    file = TFile(outfilename,'RECREATE')
    tothist.Write('pileup')
    file.Close()
    


def getDataProfile(outfilename,JSON,pileup,bins,minbias):
    """Get pileup profile in data with pileupCalc.py tool."""
    print ">>> getDataProfile(%s,%d,%s)"%(outfilename,bins,minbias)
    minbias *= 1000
    command = "pileupCalc.py -i %s --inputLumiJSON %s --minBiasXsec %d --calcMode true --maxPileupBin %d --numPileupBins %d %s"%(JSON,pileup,bins,bins,minbias,outfilename)
    print ">>>   " + command
    os.system(command)
    
    if not os.path.isfile(outfilename):
      print ">>>   Warning! getDataProfile: Could find output file %s!"%(outfilename)
      return    
    file = TFile(outfilename,'READ')
    if not file or file.IsZombie():
      print ">>>   Warning! getDataProfile: Could not open output file %s!"%(outfilename)
      return
    hist = file.Get('pileup')
    print ">>>   pileup profile in data with min. bias has a mean of %.1f"%(minbias,hist.GetMean())
    file.Close()
    


def getGenProfile(outfilename,year):
    print ">>> getGenProfile(%s):"%(year)
    if year==2017:
      # https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
      bins = [
        3.39597497605e-05, 6.63688402133e-06, 1.39533611284e-05, 3.64963078209e-05, 6.00872171664e-05, 9.33932578027e-05,
        0.000120591524486, 0.000128694546198, 0.000361697233219, 0.000361796847553, 0.000702474896113, 0.00133766053707,
        0.00237817050805,  0.00389825605651,  0.00594546732588,  0.00856825906255,  0.0116627396044,   0.0148793350787, 
        0.0179897368379,   0.0208723871946,   0.0232564170641,   0.0249826433945,   0.0262245860346,   0.0272704617569,
        0.0283301107549,   0.0294006137386,   0.0303026836965,   0.0309692426278,   0.0308818046328,   0.0310566806228,
        0.0309692426278,   0.0310566806228,   0.0310566806228,   0.0310566806228,   0.0307696426944,   0.0300103336052,
        0.0288355370103,   0.0273233309106,   0.0264343533951,   0.0255453758796,   0.0235877272306,   0.0215627588047,
        0.0195825559393,   0.0177296309658,   0.0160560731931,   0.0146022004183,   0.0134080690078,   0.0129586991411,
        0.0125093292745,   0.0124360740539,   0.0123547104433,   0.0123953922486,   0.0124360740539,   0.0124360740539,
        0.0123547104433,   0.0124360740539,   0.0123387597772,   0.0122414455005,   0.011705203844,    0.0108187105305,
        0.00963985508986,  0.00827210065136,  0.00683770076341,  0.00545237697118,  0.00420456901556,  0.00367513566191,
        0.00314570230825,  0.0022917978982,   0.00163221454973,  0.00114065309494,  0.000784838366118, 0.000533204105387,
        0.000358474034915, 0.000238881117601, 0.0001984254989,   0.000157969880198, 0.00010375646169,  6.77366175538e-05,
        4.39850477645e-05, 2.84298066026e-05, 1.83041729561e-05, 1.17473542058e-05, 7.51982735129e-06, 6.16160108867e-06,
        4.80337482605e-06, 3.06235473369e-06, 1.94863396999e-06, 1.23726800704e-06, 7.83538083774e-07, 4.94602064224e-07,
        3.10989480331e-07, 1.94628487765e-07, 1.57888581037e-07, 1.2114867431e-07,  7.49518929908e-08, 4.6060444984e-08,
        2.81008884326e-08, 1.70121486128e-08, 1.02159894812e-08, 0.0,               #0.0,
      ]
    elif year==2018:
      # https://github.com/cms-sw/cmssw/blob/CMSSW_10_4_X/SimGeneral/MixingModule/python/mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi.py
      bins = [
          4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
          3.508135e-05, 7.12608e-05,  0.0001400641, 0.0002663403, 0.0004867473,
          0.0008469,    0.001394142,  0.002169081,  0.003198514,  0.004491138,
          0.006036423,  0.007806509,  0.00976048,   0.0118498,    0.01402411,
          0.01623639,   0.01844593,   0.02061956,   0.02273221,   0.02476554,
          0.02670494,   0.02853662,   0.03024538,   0.03181323,   0.03321895,
          0.03443884,   0.035448,     0.03622242,   0.03674106,   0.0369877,
          0.03695224,   0.03663157,   0.03602986,   0.03515857,   0.03403612,
          0.0326868,    0.03113936,   0.02942582,   0.02757999,   0.02563551,
          0.02362497,   0.02158003,   0.01953143,   0.01750863,   0.01553934,
          0.01364905,   0.01186035,   0.01019246,   0.008660705,  0.007275915,
          0.006043917,  0.004965276,  0.004035611,  0.003246373,  0.002585932,
          0.002040746,  0.001596402,  0.001238498,  0.0009533139, 0.0007282885,
          0.000552306,  0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
          0.0001238161, 8.96531e-05,  6.438087e-05, 4.585302e-05, 3.23949e-05,
          2.271048e-05, 1.580622e-05, 1.09286e-05,  7.512748e-06, 5.140304e-06,
          3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
          5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
          1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
          6.119019e-08, 5.443693e-08, 4.85036e-08,  4.31486e-08,  3.822112e-08
      ]
    else:
      print ">>>   Warning! No generator pileup profile for year %s"%(year)
    
    nbins = len(bins)
    hist  = TH1F('pileup','pileup',nbins,0,nbins)
    hist.Sumw2()
    for i, binc in enumerate(bins,1):
      hist.SetBinContent(i,binc)
    
    file = TFile(outfilename,'RECREATE')
    hist.Write('pileup')
    file.Close()
    


def main():
    
    channel = args.channel
    for year in args.years:
      filename  = "MC_PileUp_%d.root"%(year)
      indir     = "/scratch/ineuteli/analysis/LQ_%d"%(year)
      if year==2016:
        samples = [
          ( 'TT', "TT",                   ),
          ( 'DY', "DYJetsToLL_M-10to50",  ),
          ( 'DY', "DYJetsToLL_M-50_reg",  ),
          ( 'DY', "DY1JetsToLL_M-50",     ),
          ( 'DY', "DY2JetsToLL_M-50",     ),
          #( 'DY', "DY3JetsToLL_M-50",     ),
          ( 'WJ', "WJetsToLNu",           ),
          ( 'WJ', "W1JetsToLNu",          ),
          ( 'WJ', "W2JetsToLNu",          ),
          ( 'WJ', "W3JetsToLNu",          ),
          #( 'WJ', "W4JetsToLNu",          ),
          ( 'ST', "ST_tW_top",            ),
          ( 'ST', "ST_tW_antitop",        ),
          #( 'ST', "ST_t-channel_top",     ),
          ( 'ST', "ST_t-channel_antitop", ),
          #( 'ST', "ST_s-channel",         ),
          ( 'VV', "WW",                   ),
          ( 'VV', "WZ",                   ),
          ( 'VV', "ZZ",                   ),
        ]
        JSON    = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
        pileup  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt"
      elif year==2017:
        JSON    = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
        pileup  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt"    
        samples = [ 
          ( 'TT', "TTTo2L2Nu",            ),
          ( 'TT', "TTToHadronic",         ),
          ( 'TT', "TTToSemiLeptonic",     ),
          ( 'DY', "DYJetsToLL_M-10to50",  ),
          ( 'DY', "DYJetsToLL_M-50",      ),
          ( 'DY', "DY1JetsToLL_M-50",     ),
          ( 'DY', "DY2JetsToLL_M-50",     ),
          ( 'DY', "DY3JetsToLL_M-50",     ),
          ( 'DY', "DY4JetsToLL_M-50",     ),
          ( 'WJ', "WJetsToLNu",           ),
          ( 'WJ', "W1JetsToLNu",          ),
          ( 'WJ', "W2JetsToLNu",          ),
          ( 'WJ', "W3JetsToLNu",          ),
          ( 'WJ', "W4JetsToLNu",          ),
          ( 'ST', "ST_tW_top",            ),
          ( 'ST', "ST_tW_antitop",        ),
          ( 'ST', "ST_t-channel_top",     ),
          ( 'ST', "ST_t-channel_antitop", ),
          #( 'ST', "ST_s-channel",         ),
          ( 'VV', "WW",                   ),
          ( 'VV', "WZ",                   ),
          ( 'VV', "ZZ",                   ),
        ]
      else:
        JSON    = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
        pileup  = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt"
        samples = [
          ( 'TT', "TTTo2L2Nu",            ),
          ( 'TT', "TTToHadronic",         ),
          ( 'TT', "TTToSemiLeptonic",     ),
          ( 'DY', "DYJetsToLL_M-10to50",  ),
          ( 'DY', "DYJetsToLL_M-50",      ),
          ( 'DY', "DY1JetsToLL_M-50",     ),
          ( 'DY', "DY2JetsToLL_M-50",     ),
          ( 'DY', "DY3JetsToLL_M-50",     ),
          #( 'DY', "DY4JetsToLL_M-50",     ),
          ( 'WJ', "WJetsToLNu",           ),
          ( 'WJ', "W1JetsToLNu",          ),
          ( 'WJ', "W2JetsToLNu",          ),
          ( 'WJ', "W3JetsToLNu",          ),
          ( 'WJ', "W4JetsToLNu",          ),
          ( 'ST', "ST_tW_top",            ),
          ( 'ST', "ST_tW_antitop",        ),
          #( 'ST', "ST_t-channel_top",     ),
          #( 'ST', "ST_t-channel_antitop", ),
          ( 'ST', "ST_s-channel",         ),
          ( 'VV', "WW",                   ),
          ( 'VV', "WZ",                   ),
          ( 'VV', "ZZ",                   ),
        ]
      
      # MC
      getMCProfile(filename,indir,samples,channel,year)
      
      # DATA
      minbiases = [ 69.2, ] # 80.0 ] #66156, 72.383
      for minbias in minbiases:
        filename = "Data_PileUp_%d_%s_new.root"%(year,str(minbias).replace('.','p'))
        getDataProfile(filename,JSON,pileup,100,minbias)
      


if __name__ == '__main__':
    print
    main()
    print ">>> done\n"
    

