#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
bool readVecBool(std::vector<bool>vec ,int index){
  for (size_t i =0; i < vec.size(); i++){
  }
  return vec[index];
}

float deltaPhi (float phi1, float phi2){
  float dPhi = fabs(phi1 - phi2);
  if(dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
  return dPhi;
}

vector<float> selectEvents(TString fn, TString treeName) {
  TFile *f = TFile::Open(fn);
  TTree *tree  = (TTree*) f->Get(treeName);
  int allEntries = tree->GetEntries();
  float sel_sum =0;
  int sel_count =0;
  bool selected = true;

  int    event_NJets=0;
  double event_HTMiss=0;
  double event_HT=0;
  double event_sumMinDR=0;
  double eventBDT_value=0;
  vector<double> *  CalibJet_BDT3weights_bib=0;
  vector<double> *  CalibJet_eta=0;
  vector<double> *  CalibJet_phi=0;
  vector<double> *  CalibJet_pT=0;
  vector<double> *  CalibJet_time=0;
  vector<double> *  CalibJet_EMF=0;
  vector<bool> * CalibJet_isGoodLLP=0;
  vector<bool> * CalibJet_isCRHLTJet=0;
  int    n_event_NJets=0;
  bool    n_event_passL1TAU60=0;
  bool    n_event_passCalRatio_cleanLLP_TAU60=0;
  double n_event_HTMiss=0;
  double n_event_HT=0;
  double n_event_sumMinDR=0;
  double n_eventBDT_value=0;
  double n_pileupEventWeight=0;
  double n_eventWeight=0;
  vector<double> *  n_CalibJet_BDT3weights_bib=0;
  vector<double> *  n_CalibJet_eta=0;
  vector<double> *  n_CalibJet_phi=0;
  vector<double> *  n_CalibJet_pT=0;
  vector<double> *  n_CalibJet_time=0;
  vector<double> *  n_CalibJet_EMF=0;
  vector<bool> * n_CalibJet_isGoodLLP=0;
  vector<bool> * n_CalibJet_isCRHLTJet=0;
  vector<int> * BDT3weights_signal_cleanJet_index=0;
  
  tree->SetBranchAddress("CalibJet_isGoodLLP",&CalibJet_isGoodLLP);
  tree->SetBranchAddress("CalibJet_isCRHLTJet",&CalibJet_isCRHLTJet);
  tree->SetBranchAddress("CalibJet_time",&CalibJet_time);
  tree->SetBranchAddress("CalibJet_pT",&CalibJet_pT);
  tree->SetBranchAddress("CalibJet_phi",&CalibJet_phi);
  tree->SetBranchAddress("CalibJet_eta",&CalibJet_eta);
  tree->SetBranchAddress("CalibJet_EMF",&CalibJet_EMF);
  tree->SetBranchAddress("CalibJet_BDT3weights_bib",&CalibJet_BDT3weights_bib);
  tree->SetBranchAddress("event_NJets",&event_NJets);
  tree->SetBranchAddress("event_sumMinDR",&event_sumMinDR);
  tree->SetBranchAddress("event_HT",&event_HT);
  tree->SetBranchAddress("event_HTMiss",&event_HTMiss);
  tree->SetBranchAddress("eventBDT_value",& eventBDT_value);
  tree->SetBranchAddress("BDT3weights_signal_cleanJet_index",&BDT3weights_signal_cleanJet_index);
  
  for (int i =0 ; i<allEntries; i++){
	tree->GetEntry(i);
	float w=n_eventWeight*n_pileupEventWeight;
	bool trigger = n_event_passL1TAU60 * n_event_passCalRatio_cleanLLP_TAU60 ;
	if (trigger){
		if (!(event_NJets > 1)) continue;
	  if (!(eventBDT_value >0.05)) continue;
		if (!(((CalibJet_pT)->at(BDT3weights_signal_cleanJet_index->at(0))>40) && (fabs((CalibJet_eta)->at(BDT3weights_signal_cleanJet_index->at(0)))<2.5) && ((CalibJet_isGoodLLP)->at(BDT3weights_signal_cleanJet_index->at(0))))) continue; 
		if (!(((CalibJet_pT)->at(BDT3weights_signal_cleanJet_index->at(1))>40) && (fabs((CalibJet_eta)->at(BDT3weights_signal_cleanJet_index->at(1)))<2.5) && ((CalibJet_isGoodLLP)->at(BDT3weights_signal_cleanJet_index->at(1))))) continue; 
		if (!((CalibJet_isCRHLTJet)->at(BDT3weights_signal_cleanJet_index->at(0)) || (CalibJet_isCRHLTJet)->at(BDT3weights_signal_cleanJet_index->at(1))))  continue;
		if (!(fabs(deltaPhi((CalibJet_phi)->at(BDT3weights_signal_cleanJet_index->at(0)),(CalibJet_phi)->at(BDT3weights_signal_cleanJet_index->at(1)))> 0.75))) continue; 
		if (!(((CalibJet_time)->at(BDT3weights_signal_cleanJet_index->at(0))>-3 && (CalibJet_time)->at(BDT3weights_signal_cleanJet_index->at(0))<15) )) continue; 
		if (!(((CalibJet_time)->at(BDT3weights_signal_cleanJet_index->at(1))>-3 && (CalibJet_time)->at(BDT3weights_signal_cleanJet_index->at(1))<15) )) continue; 
		if (!( (event_HTMiss/event_HT) < 0.95 )) continue; 
		int nextEvent=0;
		for (size_t i=0; i < CalibJet_BDT3weights_bib->size(); i++){if ((CalibJet_BDT3weights_bib)->at(i) >0.6) nextEvent=1; }
		if (nextEvent) continue;
		if (!(((CalibJet_pT)->at(BDT3weights_signal_cleanJet_index->at(0))>120))) continue;
		if (!(event_sumMinDR > 0.5)) continue; 
		sel_sum+=w;
		sel_count+=1;
	}
}
  vector<float> result;
  result.push_back(sel_sum);
  result.push_back(sel_count);
  delete tree;
  return result;
}
