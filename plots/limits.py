#!/usr/bin/env python
import ROOT as r
import argparse
import os
from array import array
r.gROOT.SetBatch(1)

#Parse command line
parser = argparse.ArgumentParser(description='Make limit plot with result of extrapolation code')
parser.add_argument('-e', '--extrapFile', help='Path to root file containing extrapolation results')
args = parser.parse_args()
parser.print_help()

#Open the file
inFile = r.TFile(args.extrapFile)


#Set up the canvas
canvas = r.TCanvas('canvas', 'canvas', 600, 600)

#Take central limit as 
