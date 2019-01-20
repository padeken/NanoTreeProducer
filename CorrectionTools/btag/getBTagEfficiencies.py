#! /usr/bin/env python
# Author: Izaak Neutelings (January 2019)

import os, sys
from argparse import ArgumentParser
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gStyle, gROOT, TFile, TTree, TH2F, TCanvas
gStyle.SetOptStat(False)
gROOT.SetBatch(True)

argv = sys.argv
description = '''This script extracts histograms to create b tag efficiencies.'''
parser = ArgumentParser(prog="pileup",description=description,epilog="Succes!")
parser.add_argument('-y', '--year',    dest='years', choices=[2017,2018], type=int, nargs='+', default=[2017], action='store',
                                       help="select year" )
parser.add_argument('-c', '--channel', dest='channels', choices=['eletau','mutau','tautau'], type=str, nargs='+', default=['mutau'], action='store',
                                       help="channels to submit" )
parser.add_argument('-t', '--tagger',  dest='taggers', choices=['CSVv2','DeepCSV'], type=str, nargs='+', default=['CSVv2'], action='store',
                                       help="channels to submit" )
parser.add_argument('-w', '--wp',      dest='wps', choices=['loose','medium','tight'], type=str, nargs='+', default=['medium'], action='store',
                                       help="channels to submit" )
parser.add_argument('-p', '--plot',    dest="plot", default=False, action='store_true', 
                                       help="plot efficiencies" )
parser.add_argument('-v', '--verbose', dest="verbose", default=False, action='store_true', 
                                       help="print verbose" )
args = parser.parse_args()



def getBTagEfficiencies(tagger,wp,outfilename,indir,samples,channel,plot=False):
    """Get pileup profile in MC by adding Pileup_nTrueInt histograms from a given list of samples."""
    print ">>> getBTagEfficiencies(%s)"%(outfilename)
    
    # GET HISTOGRAMS
    nhists  = { }
    hists   = { }
    histdir = 'btag'
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
        histpath = "%s/%s"%(histdir,histname)
        hist = file.Get(histpath)
        if not hist:
          print ">>>   Warning! getBTagEfficiencies: Could not open histogram '%s' in %s"%(histpath,filename)        
          dir = file.Get(histdir)
          if dir: dir.ls()
          continue
        if hists[histname]==None:
          hists[histname] = hist.Clone(histname)
          hists[histname].SetDirectory(0)
          nhists[histname] = 1
        else:
          hists[histname].Add(hist)
          nhists[histname] += 1
      file.Close()
    if len(nhists)>0:
      print ">>>   added %d MC hists:"%(sum(nhists[n] for n in nhists))
      for histname, nhist in nhists.iteritems():
        print ">>>     %-26s%2d"%(histname+':',nhist)
    else:
      print ">>>   no histograms added !"%(sum(nhists[n] for n in nhists))
      return
    
    # SAVE HISTOGRAMS
    print ">>>   writing to %s..."%(outfilename)
    file = TFile(outfilename,'RECREATE')
    for histname, hist in hists.iteritems():
      if 'all' in histname:
        continue
      histname_all = histname+'_all'
      histname_eff = 'eff_'+histname
      print ">>>      writing %s..."%(histname)
      print ">>>      writing %s..."%(histname_all)
      print ">>>      writing %s..."%(histname_eff)
      hist_all = hists[histname_all]
      hist_eff = hist.Clone(histname_eff)
      hist_eff.SetTitle(histname_eff)
      hist_eff.Divide(hist_all)
      hist.Write(histname)
      hist_all.Write(histname_all)
      hist_eff.Write(histname_eff)
      if plot:
        plot2D(hist_eff)
    file.Close()


def plot2D(hist):
    """Plot efficiency."""
    dir    = ensureDirectory('plots')
    name   = "%s/%s"%(dir,hist.GetName())
    xtitle = 'jet p_{T} [GeV]'
    ytitle = 'jet #eta'
    ztitle = 'b tag efficiencies'
    
    canvas = TCanvas('canvas','canvas',100,100,800,700)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetTopMargin(  0.07 ); canvas.SetBottomMargin( 0.13 )
    canvas.SetLeftMargin( 0.12 ); canvas.SetRightMargin(  0.17 )
    canvas.SetTickx(0); canvas.SetTicky(0)
    canvas.SetGrid()
    canvas.cd()
    
    hist.Draw('COLZ')
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle(ytitle)
    hist.GetZaxis().SetTitle(ztitle)
    hist.GetXaxis().SetLabelSize(0.048)
    hist.GetYaxis().SetLabelSize(0.048)
    hist.GetZaxis().SetLabelSize(0.048)
    hist.GetXaxis().SetTitleSize(0.058)
    hist.GetYaxis().SetTitleSize(0.058)
    hist.GetZaxis().SetTitleSize(0.056)
    hist.GetXaxis().SetTitleOffset(1.03)
    hist.GetYaxis().SetTitleOffset(1.04)
    hist.GetZaxis().SetTitleOffset(1.02)
    hist.SetMinimum(0.0)
    hist.SetMaximum(1.0)
    
    canvas.SaveAs(name+'.pdf')
    canvas.SaveAs(name+'.png')
    

def ensureDirectory(dirname):
  """Make directory if it does not exist."""
  if not os.path.exists(dirname):
    os.makedirs(dirname)
    print '>>> made directory "%s"'%(dirname)
    if not os.path.exists(dirname):
      print '>>> failed to make directory "%s"'%(dirname)
  return dirname
  

def main():
    
    missing  = [ "DY4Jets", "TTToSemi", "WJets", "W1Jets", "WZ", "ST_t-channel" ]
    samples0 = [
      ( 'DY', "DYJetsToLL_M-50",      "Drell-Yan 50",         ),
      ( 'DY', "DY1JetsToLL_M-50",     "Drell-Yan 1J 50",      ),
      ( 'DY', "DY2JetsToLL_M-50",     "Drell-Yan 2J 50",      ),
      ( 'DY', "DY3JetsToLL_M-50",     "Drell-Yan 3J 50",      ),
      ( 'DY', "DY4JetsToLL_M-50",     "Drell-Yan 4J 50",      ),
      ( 'TT', "TTTo2L2Nu",            "ttbar 2l2#nu",         ),
      ( 'TT', "TTToHadronic",         "ttbar hadronic",       ),
      ( 'TT', "TTToSemiLeptonic",     "ttbar semileptonic",   ),
      ( 'WJ', "WJetsToLNu",           "W + jets",             ),
      ( 'WJ', "W1JetsToLNu",          "W + 1J",               ),
      ( 'WJ', "W2JetsToLNu",          "W + 2J",               ),
      ( 'WJ', "W3JetsToLNu",          "W + 3J",               ),
      ( 'WJ', "W4JetsToLNu",          "W + 4J",               ),
      ( 'ST', "ST_t-channel_top",     "ST t-channel top",     ),
      ( 'ST', "ST_t-channel_antitop", "ST t-channel antitop", ),
      #( 'ST', "ST_s-channel",         "ST s-channel",         ),
      ( 'ST', "ST_tW_top",            "ST tW",                ),
      ( 'ST', "ST_tW_antitop",        "ST atW",               ),
      ( 'VV', "WW",                   "WW",                   ),
      ( 'VV', "WZ",                   "WZ",                   ),
      ( 'VV', "ZZ",                   "ZZ",                   ),
    ]
    samples_2018 = [(d,n,t) for (d,n,t) in samples0 if not any(s in n for s in missing)]
    
    for year in args.years:
      samples = samples_2018 if year==2018 else samples0
      for channel in args.channels:
        for tagger in args.taggers:
          for wp in args.wps:
            filename = "%s_%d_eff.root"%(tagger,year)
            indir    = "/scratch/ineuteli/analysis/LQ_%d"%(year)
            getBTagEfficiencies(tagger,wp,filename,indir,samples,channel,plot=args.plot)    
    


if __name__ == '__main__':
    print
    main()
    print ">>> done\n"
    

