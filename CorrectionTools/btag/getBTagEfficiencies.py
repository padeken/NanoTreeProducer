#! /usr/bin/env python
# Author: Izaak Neutelings (January 2019)

import os, sys
from argparse import ArgumentParser
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TFile, TTree, TH2F

argv = sys.argv
description = '''This script extracts histograms to create b tag efficiencies.'''
parser = ArgumentParser(prog="pileup",description=description,epilog="Succes!")
parser.add_argument('-y', '--year',    dest='years', choices=[2017,2018], type=int, nargs='+', default=[2017], action='store',
                                       help="select year" )
parser.add_argument('-c', '--channel', dest='channels', choices=['eletau','mutau','tautau'], type=str, nargs='+', default=['mutau'], action='store',
                                       help="channels to submit" )
parser.add_argument('-t', '--tagger',  dest='taggers', choices=['CSVv2','DeepCSV'], type=str, nargs='+', default=['CSVvs'], action='store',
                                       help="channels to submit" )
parser.add_argument('-w', '--wp',      dest='wps', choices=['loose','medium','tight'], type=str, nargs='+', default=['medium'], action='store',
                                       help="channels to submit" )
parser.add_argument('-v', '--verbose', dest="verbose", default=False, action='store_true', 
                                       help="print verbose" )
args = parser.parse_args()



def getBTagEfficiencies(tagger,wp,outfilename,indir,samples,channel):
    """Get pileup profile in MC by adding Pileup_nTrueInt histograms from a given list of samples."""
    print ">>> getBTagEfficiencies(%s)"%(outfilename)
    nhists = 0
    hists  = { }
    for flavor in ['b','c','udsg']:
      histname = '%s_%s_%s'%(tagger,flavor,wp)
      hists[histname] = None
      hists[histname+'_all'] = None   
    for subdir, samplename, sampletitle in samples:
      filename = "%s/%s/%s_%s.root"%(indir,subdir,samplename,channel)
      print ">>>   %s"%(filename)
      file = TFile(filename,'READ')
      if not file or file.IsZombie():
        print ">>>   Warning! getBTagEfficiencies: Could not open %s"%(filename)
        continue
      for histname in hists:
        hist = file.Get(histname)
        if not hist:
          print ">>>   Warning! getBTagEfficiencies: Could not open histogram %s in %s"%(histname,filename)      
          continue
        if hists[histname]==None:
          hists[histname] = hist.Clone(histname)
          hists[histname].SetDirectory(0)
          nhists += 1
        else:
          hists[histname].Add(hist)
          nhists += 1
        file.Close()
    print ">>>   added %d MC hists"%(nhists)
    
    file = TFile(outfilename,'RECREATE')
    for histname, hist in hists.iteritems():
      if 'all' in histname:
        continue
      histname_all = histname+'_all'
      histname_eff = 'eff_'+histname
      print ">>>   writing %s..."%(histname)
      print ">>>   writing %s..."%(histname_all)
      print ">>>   writing %s..."%(histname_eff)
      hist_all = hists[histname_all]
      hist_eff = hist.Clone(histname_eff)
      hist_eff.SetTitle(histname_eff)
      hist_eff.Divide(hist_all)
      hist.Write(histname)
      hist_all.Write(histname_all)
      hist_eff.Write(histname_eff)
    file.Close()
    


def main():
    
    missing  = [ "DY4Jets", "TTToSemi", "WJets", "W1Jets", "WZ" ]
    samples0 = [
      ( 'DY', "DYJetsToLL_M-50",   "Drell-Yan 50",       ),
      ( 'DY', "DY1JetsToLL_M-50",  "Drell-Yan 1J 50",    ),
      ( 'DY', "DY2JetsToLL_M-50",  "Drell-Yan 2J 50",    ),
      ( 'DY', "DY3JetsToLL_M-50",  "Drell-Yan 3J 50",    ),
      ( 'DY', "DY4JetsToLL_M-50",  "Drell-Yan 4J 50",    ),
      ( 'TT', "TTTo2L2Nu",         "ttbar 2l2#nu",       ),
      ( 'TT', "TTToHadronic",      "ttbar hadronic",     ),
      ( 'TT', "TTToSemiLeptonic",  "ttbar semileptonic", ),
      ( 'WJ', "WJetsToLNu",        "W + jets",           ),
      ( 'WJ', "W1JetsToLNu",       "W + 1J",             ),
      ( 'WJ', "W2JetsToLNu",       "W + 2J",             ),
      ( 'WJ', "W3JetsToLNu",       "W + 3J",             ),
      ( 'WJ', "W4JetsToLNu",       "W + 4J",             ),
      ( 'ST', "ST_s-channel",      "ST s-channel",       ),
      ( 'ST', "ST_tW_top",         "ST tW",              ),
      ( 'ST', "ST_tW_antitop",     "ST atW",             ),
      ( 'VV', "WW",                "WW",                 ),
      ( 'VV', "WZ",                "WZ",                 ),
      ( 'VV', "ZZ",                "ZZ",                 ),
    ]
    samples_2018 = [(d,n,t) for (d,n,t) in samples0 if not any(s in n for s in missing)]
    
    for year in args.years:
      samples = samples_2018 if year==2018 else samples0
      for channel in args.channels:
        for taggger in args.taggers:
          for wp in args.wps:
            filename = "%s_%d_eff.root"%(tagger,year)
            indir    = "/scratch/ineuteli/analysis/LQ_%d"%(year)
            getBTagEfficiencies(tagger,wp,filename,indir,samples,channel)    
    


if __name__ == '__main__':
    print
    main()
    print ">>> done\n"
    

