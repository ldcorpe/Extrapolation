#!/usr/bin/env python
import ROOT as r
import argparse
import os
from array import array
r.gROOT.SetBatch(1)

#Import all functions from the selection macro 
r.gSystem.Load('CalRSelection2017_C.so')

parser = argparse.ArgumentParser(description='Slim MC file to use in extrapolation code.')
#parser.add_argument('-s', metavar='sample', dest='samplename', nargs=1, help='Sample to process, e.g. mH600_mS150_lt5m')

parser.add_argument('-s', '--sample', help='Sample to process, e.g. mH600_mS150_lt5m.root')
parser.add_argument('-o', '--outfile', help='Name of output root file')
parser.add_argument('-n', '--nevents', help='Number of events to process (default -1)', default='-1', type=int)
parser.add_argument('-m', '--mH', help='Boson mass, defines which selection to use', type=int)
args = parser.parse_args()
parser.print_help()

#Open the right file
treeName = 'recoTree'
inFile = r.TFile(args.sample)
tree = inFile.Get(treeName)

#Define the output file and tree
f = r.TFile( args.outfile, 'recreate' ) 
t = r.TTree( 'extrapTree', 'Used as input for the extrapolation' ) 

#Set up the branches in the output tree
eventNumber = array( 'i', [ 0 ] )
PassedCalRatio = array( 'i', [ 0 ] )

llp1_pt = array( 'd', [ 0. ] )
llp2_pt = array( 'd', [ 0. ] )
llp1_eta = array( 'd', [ 0. ] )
llp2_eta = array( 'd', [ 0. ] )
llp1_phi = array( 'd', [ 0. ] )
llp2_phi = array( 'd', [ 0. ] )
llp1_E = array( 'd', [ 0. ] )
llp2_E = array( 'd', [ 0. ] )
llp1_Lxy = array( 'd', [ 0. ] )
llp2_Lxy = array( 'd', [ 0. ] )

event_weight = array( 'd', [ 0. ] )

RegionA = array( 'i', [ 0 ] )
RegionB = array( 'i', [ 0 ] )
RegionC = array( 'i', [ 0 ] )
RegionD = array( 'i', [ 0 ] )

t.Branch( 'eventNumber', eventNumber, 'eventNumber/I')
t.Branch( 'PassedCalRatio', PassedCalRatio, 'PassedCalRatio/I')

t.Branch( 'llp1_pt', llp1_pt, 'llp1_pt/D')
t.Branch( 'llp2_pt', llp2_pt, 'llp2_pt/D')
t.Branch( 'llp1_eta', llp1_eta, 'llp1_eta/D')
t.Branch( 'llp2_eta', llp2_eta, 'llp2_eta/D')
t.Branch( 'llp1_phi', llp1_phi, 'llp1_phi/D')
t.Branch( 'llp2_phi', llp2_phi, 'llp2_phi/D')
t.Branch( 'llp1_E', llp1_E, 'llp1_E/D')
t.Branch( 'llp2_E', llp2_E, 'llp2_E/D')
t.Branch( 'llp1_Lxy', llp1_Lxy, 'llp1_Lxy/D')
t.Branch( 'llp2_Lxy', llp2_Lxy, 'llp2_Lxy/D')

t.Branch( 'event_weight', event_weight, 'event_weight/D')

t.Branch( 'RegionA', RegionA, 'RegionA/I')
t.Branch( 'RegionB', RegionB, 'RegionB/I')
t.Branch( 'RegionC', RegionC, 'RegionC/I')
t.Branch( 'RegionD', RegionD, 'RegionD/I')

#Keep track of the number of events
counter = 0
maxcount = args.nevents

#For each event:
# - check if it passes the trigger and selection
# - where it lies in the ABCD plane
# - write out to new file
for ev in tree:
  
  counter += 1
  if maxcount > 0:
    if counter > maxcount: break
  
  passedTrigger = False
  passedTrigger = ev.event_passCalRatio_cleanLLP_TAU60
  
  isSelected = True

  PassBibJetVeto=True
  for i in range (ev.CalibJet_BDT3weights_bib.size()):
    if (ev.CalibJet_BDT3weights_bib[i] > 0.6 and ev.CalibJet_pT[i] > 40 and abs(ev.CalibJet_eta[i]) < 2.5 and ev.CalibJet_isGoodLLP[i] ): PassBibJetVeto = False

  if ev.event_NCleanJets>1 and PassBibJetVeto and args.mH>399: 
    isSelected = r.event_selection(2, ev.CalibJet_pT[ev.BDT3weights_signal_cleanJet_index[0]],r.readVecBool(ev.CalibJet_isGoodLLP,ev.BDT3weights_signal_cleanJet_index[0]),r.readVecBool(ev.CalibJet_isGoodLLP,ev.BDT3weights_signal_cleanJet_index[1]),ev.CalibJet_time[ev.BDT3weights_signal_cleanJet_index[0]],ev.CalibJet_time[ev.BDT3weights_signal_cleanJet_index[1]],ev.CalibJet_time[ev.BDT3weights_bib_cleanJet_index[0]],ev.CalibJet_time[ev.BDT3weights_bib_cleanJet_index[1]],ev.event_MHToHT,ev.event_sumMinDR,ev.CalibJet_minDRTrkpt2[ev.BDT3weights_signal_cleanJet_index[0]],ev.CalibJet_minDRTrkpt2[ev.BDT3weights_signal_cleanJet_index[1]],r.readVecBool(ev.CalibJet_isCRHLTJet,ev.BDT3weights_signal_cleanJet_index[0]),r.readVecBool(ev.CalibJet_isCRHLTJet,ev.BDT3weights_signal_cleanJet_index[1]),ev.CalibJet_logRatio[ev.BDT3weights_signal_cleanJet_index[0]],ev.CalibJet_logRatio[ev.BDT3weights_signal_cleanJet_index[1]],ev.eventBDT_value)

  if ev.event_NCleanJets>1 and PassBibJetVeto and args.mH<399 : isSelected = r.event_selection(1, ev.CalibJet_pT[ev.BDT3weights_signal_cleanJet_index[0]],r.readVecBool(ev.CalibJet_isGoodLLP,ev.BDT3weights_signal_cleanJet_index[0]),r.readVecBool(ev.CalibJet_isGoodLLP,ev.BDT3weights_signal_cleanJet_index[1]),ev.CalibJet_time[ev.BDT3weights_signal_cleanJet_index[0]],ev.CalibJet_time[ev.BDT3weights_signal_cleanJet_index[1]],ev.CalibJet_time[ev.BDT3weights_bib_cleanJet_index[0]],ev.CalibJet_time[ev.BDT3weights_bib_cleanJet_index[1]],ev.event_MHToHT,ev.event_sumMinDR,ev.CalibJet_minDRTrkpt2[ev.BDT3weights_signal_cleanJet_index[0]],ev.CalibJet_minDRTrkpt2[ev.BDT3weights_signal_cleanJet_index[1]],r.readVecBool(ev.CalibJet_isCRHLTJet,ev.BDT3weights_signal_cleanJet_index[0]),r.readVecBool(ev.CalibJet_isCRHLTJet,ev.BDT3weights_signal_cleanJet_index[1]),ev.CalibJet_logRatio[ev.BDT3weights_signal_cleanJet_index[0]],ev.CalibJet_logRatio[ev.BDT3weights_signal_cleanJet_index[1]],ev.eventBDT_value)

  region = 0

  if (ev.event_NCleanJets>1): region = r.event_ABCD_plane(ev.eventBDT_value,ev.event_sumMinDR)

  llp1_E[0] = ev.LLP_E[0]
  llp2_E[0] = ev.LLP_E[1]
  llp1_eta[0] = ev.LLP_eta[0]
  llp2_eta[0] = ev.LLP_eta[1]
  llp1_pt[0] = ev.LLP_pT[0]
  llp2_pt[0] = ev.LLP_pT[1]
  llp1_phi[0] = ev.LLP_phi[0]
  llp2_phi[0] = ev.LLP_phi[1]
  llp1_Lxy[0] = ev.LLP_Lxy[0]
  llp2_Lxy[0] = ev.LLP_Lxy[1]

  PassedCalRatio[0] = passedTrigger
  eventNumber[0] = ev.eventNumber
    
  event_weight[0] = ev.eventWeight * abs(ev.pileupEventWeight)

  RegionA[0],RegionB[0],RegionC[0],RegionD[0] = 0,0,0,0
  if (region == 1 and passedTrigger and isSelected): RegionA[0] = 1
  if (region == 2 and passedTrigger and isSelected): RegionB[0] = 1
  if (region == 3 and passedTrigger and isSelected): RegionC[0] = 1
  if (region == 4 and passedTrigger and isSelected): RegionD[0] = 1
#  print RegionA,RegionB,RegionC,RegionD

  t.Fill()


f.Write()
f.Close()
