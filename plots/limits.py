#!/usr/bin/env python
import ROOT as r
import argparse
import os
from array import array
r.gROOT.SetBatch(1)
r.gStyle.SetOptLogy(1)
r.gStyle.SetOptLogx(1)

#Parse command line
parser = argparse.ArgumentParser(description='Make limit plot with result of extrapolation code')
parser.add_argument('-e', '--extrapFile', help='Path to root file containing extrapolation results')
args = parser.parse_args()
parser.print_help()

#Open the file
inFile = r.TFile(args.extrapFile)

#Set up the canvas
#canvas = r.TCanvas('canvas', 'canvas', 600, 600)

# Get histograms
centralExp = inFile.Get('xsec_BR_95CL')
plus1 = inFile.Get('xsec_BR_events__p1')
plus2 =inFile.Get('xsec_BR_events__p2')
minus1 =inFile.Get('xsec_BR_events__n1')
minus2 =inFile.Get('xsec_BR_events__n2')
centralObs = inFile.Get('xsec_BR__limit')

#GET CENTRAL VALUES OF EACH AND FILL ARRAYS TO THEN SUBTRACT ETC

centralExp.Print()
plus1.Print()

centralExp.Draw()
r.gPad.SaveAs("testCentre.png")
plus1.Draw()
r.gPad.SaveAs("testPlus.png")

diff_plus1 = r.TH1D('diff_plus1', 'diff_plus1', centralExp.GetNbinsX(), centralExp.GetXaxis().GetXmin(), centralExp.GetXaxis().GetXmax())
diff_plus2 = r.TH1D('diff_plus2', 'diff_plus2', centralExp.GetNbinsX(), centralExp.GetXaxis().GetXmin(), centralExp.GetXaxis().GetXmax())
diff_minus1 = r.TH1D('diff_minus1', 'diff_minus2', centralExp.GetNbinsX(), centralExp.GetXaxis().GetXmin(), centralExp.GetXaxis().GetXmax())
diff_minus2 = r.TH1D('diff_minus2', 'diff_minus2', centralExp.GetNbinsX(), centralExp.GetXaxis().GetXmin(), centralExp.GetXaxis().GetXmax())

#diff_plus1.Print()

diff_plus1.Add(centralExp, plus1, -1, 1)
diff_plus2 = plus2.Add(centralExp, -1)
diff_minus1 = minus1.Add(centralExp, -1)
diff_minus2 = minus2.Add(centralExp, -1)

diff_plus1.Draw()

#g_centralExp_1 = r.TGraphAsymmErrors(centralExp)
#g_centralExp_2 = r.TGraphAsymmErrors(centralExp)
#g_centralObs = r.TGraph(centralObs)


r.gPad.SaveAs("testDiff.png")
#Take central limit as 
