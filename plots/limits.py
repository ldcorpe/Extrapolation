#!/usr/bin/env python
import ROOT as r
import argparse
import os
from array import array
#r.gROOT.SetBatch(1)
r.gStyle.SetOptLogy(1)
r.gStyle.SetOptLogx(1)

#Parse command line
parser = argparse.ArgumentParser(description='Make limit plot with result of extrapolation code')
parser.add_argument('-e', '--extrapFile', help='Path to root file containing extrapolation results')
args = parser.parse_args()
parser.print_help()

#Method to get TGraphAsymmErrors from TH1's
def getTGraphAsymErrs(nominalHist,upHist,dnHist):
  result = r.TGraphAsymmErrors()
  for ibin in range(0,nominalHist.GetNbinsX()+1):
     nominalY=nominalHist.GetBinContent(ibin)
     nominalX=nominalHist.GetBinCenter(ibin)
     #XErr=nomHist.GetBinWidth(ibin)/2.
     YErrUp=upHist.GetBinContent(ibin)
     YErrDn=dnHist.GetBinContent(ibin)
     result.SetPoint(ibin,nominalX,nominalY)
     result.SetPointEYhigh(ibin,YErrUp)
     result.SetPointEYlow(ibin,YErrDn)
     print "bin ", ibin, " center ",  nominalHist.GetBinCenter(ibin) , " eUp ", YErrUp , "eDn ", YErrDn, " eAll "
  return result

def getTGraph(hist):
  result = r.TGraph()
  for ibin in range(0,hist.GetNbinsX()+1):
    histY=hist.GetBinContent(ibin)
    histX=hist.GetBinCenter(ibin)
    result.SetPoint(ibin,histX,histY)
  return result

#Open the file
inFile = r.TFile(args.extrapFile)

#Set up the canvas
canvas = r.TCanvas('canvas', 'canvas', 600, 600)

# Get histograms
centralExp = inFile.Get('xsec_BR_95CL')
plus1 = inFile.Get('xsec_BR_events__p1')
plus2 =inFile.Get('xsec_BR_events__p2')
minus1 =inFile.Get('xsec_BR_events__n1')
minus2 =inFile.Get('xsec_BR_events__n2')
centralObs = inFile.Get('xsec_BR_events__limit')

centralObs.Print()
plus1.Print()

tg_1s = getTGraphAsymErrs(centralExp,plus1,minus1)
tg_2s = getTGraphAsymErrs(centralExp,plus2,minus2)

tg_Obs = getTGraph(centralObs)

tg_1s.Draw("ape")
#tg_2s.Draw("4c same")
#tg_Obs.Draw("c same")

#canvas.Modified()
#canvas.Update()

canvas.SaveAs("test.png")
