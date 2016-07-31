// All code for reading through the muon tree TTree object.
#include "muon_tree_processor.h"

using namespace std;

// Open up and fetch the root file
muon_tree_processor::muon_tree_processor(const string &filename)
{
	_file = unique_ptr<TFile>(TFile::Open(filename.c_str(), "READ"));
	if (!_file->IsOpen()) {
		throw runtime_error("Unable to open ROOT input file " + filename + "!");
	}
	_tree = static_cast<TTree*>(_file->Get("extrapTree"));
	if (!_tree) {
		throw runtime_error("Unable to find extrapTree TTree in file " + filename + "!");
	}

	// Next, link everything up
	_tree->SetBranchAddress("PassedCalRatio", &(_tree_data.PassedCalRatio));
	_tree->SetBranchAddress("vpi1_pt", &(_tree_data.vpi1_pt));
	_tree->SetBranchAddress("vpi1_eta", &(_tree_data.vpi1_eta));
	_tree->SetBranchAddress("vpi1_phi", &(_tree_data.vpi1_phi));
	_tree->SetBranchAddress("vpi1_E", &(_tree_data.vpi1_E));
	_tree->SetBranchAddress("vpi1_Lxy", &(_tree_data.vpi1_Lxy));
	_tree->SetBranchAddress("vpi2_pt", &(_tree_data.vpi2_pt));
	_tree->SetBranchAddress("vpi2_eta", &(_tree_data.vpi2_eta));
	_tree->SetBranchAddress("vpi2_phi", &(_tree_data.vpi2_phi));
	_tree->SetBranchAddress("vpi2_E", &(_tree_data.vpi2_E));
	_tree->SetBranchAddress("vpi2_Lxy", &(_tree_data.vpi2_Lxy));
	_tree->SetBranchAddress("event_weight", &(_tree_data.weight));
	_tree->SetBranchAddress("RegionA", &(_tree_data.RegionA));
	_tree->SetBranchAddress("RegionB", &(_tree_data.RegionB));
	_tree->SetBranchAddress("RegionC", &(_tree_data.RegionC));
	_tree->SetBranchAddress("RegionD", &(_tree_data.RegionD));
}


muon_tree_processor::~muon_tree_processor()
{
}
