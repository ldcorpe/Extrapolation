#ifndef __muon_tree_processor__
#define __muon_tree_processor__

#include <TTree.h>
#include <TFile.h>

#include <string>
#include <memory>

class muon_tree_processor
{
public:
	muon_tree_processor(const std::string &filename);
	~muon_tree_processor();

	struct eventInfo {
		bool PassedCalRatio;
		double vpi1_pt;
		double vpi1_eta;
		double vpi1_phi;
		double vpi1_E;
		double vpi1_Lxy;
		double vpi2_pt;
		double vpi2_eta;
		double vpi2_phi;
		double vpi2_E;
		double vpi2_Lxy;
		double weight;
		bool IsInSignalRegion;
	};

	//Call f for each entry in the ntuple
	template<class UnaryFunction>
	void process_all_entries(UnaryFunction f) const
	{
		auto n_entries = _tree->GetEntries();
		for (decltype(n_entries) i = 0; i < n_entries; i++) {
			_tree->GetEntry(i);
			f(_tree_data);
		}
	}

private:
	std::unique_ptr<TTree> _tree;
	std::unique_ptr<TFile> _file;

	eventInfo _tree_data;
};

#endif
