// Cache values of beta that are calculated.
// This is a way of making the calculations more efficient, but
// also keeping the flow of the code clean.
#pragma once

#include "muon_tree_processor.h"
#pragma warning (push)
#pragma warning (disable: 4244)
#include "TLorentzVector.h"
#pragma warning (pop)

class beta_cache {
public:
	inline beta_cache(const muon_tree_processor::eventInfo &entry)
		: _gotit(false), _entry(entry)
	{}

	inline double beta1() { calc();  return _b1; }
	inline double beta2() { calc();  return _b2; }

private:
	bool _gotit;
	double _b1, _b2;
	const muon_tree_processor::eventInfo &_entry;

	void calc()
	{
		if (_gotit)
			return;
		_gotit = true;

		TLorentzVector vpi1, vpi2;
		vpi1.SetPtEtaPhiE(_entry.vpi1_pt, _entry.vpi1_eta, _entry.vpi1_phi, _entry.vpi1_E);
		vpi2.SetPtEtaPhiE(_entry.vpi2_pt, _entry.vpi2_eta, _entry.vpi2_phi, _entry.vpi2_E);

		_b1 = vpi1.Beta();
		_b2 = vpi2.Beta();
	}
};
