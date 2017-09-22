#ifndef __muon_tree_processor__
#define __muon_tree_processor__

#include <TTree.h>
#include <TFile.h>

#include <string>
#include <memory>
#include <vector>
#include <functional>

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
		bool RegionA;
		bool RegionB;
		bool RegionC;
		bool RegionD;
	};

	// This function will be called before the entries are processed. Only if it returns true will
	// your process function be called. If you want to avoid calling them, pass a special argument
	// to process_all_entries.
	void add_preselection(const std::function<bool(const muon_tree_processor::eventInfo&)> &func)
	{
		_preselection_list.push_back(func);
	}

	//Call f for each entry in the ntuple
	template<class UnaryFunction>
	void process_all_entries(UnaryFunction f, bool apply_preselection = true) const
	{
		auto n_entries = _tree->GetEntries();
		for (decltype(n_entries) i = 0; i < n_entries; i++) {
			_tree->GetEntry(i);
			bool good_event = true;
			if (apply_preselection) {
				for (auto &f : _preselection_list) {
					good_event = good_event && f(_tree_data);
				}
			}
			if (good_event) {
				f(_tree_data);
			}
		}
	}

private:
	// The tree we are scanning over
	// The tree is owned by the file, so we don't need to delete it; it will be deleted
	// for us when the file is removed.
	TTree *_tree;
	std::unique_ptr<TFile> _file;

	eventInfo _tree_data;
	std::vector <std::function<bool(const eventInfo&)> > _preselection_list;
};

#endif
