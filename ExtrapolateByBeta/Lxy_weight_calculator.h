//
// Calculate and look up the probability that a particular lxy will pass
// the analysis.
//
#pragma once

#include "muon_tree_processor.h"

#include "TH2D.h"

#include <memory>

class Lxy_weight_calculator
{
public:
	// Look up and calculate the weight table
	Lxy_weight_calculator(const muon_tree_processor &reader);
	~Lxy_weight_calculator();

	// Region can be 0 == A, 1 == B, 2 == C, 3 == D
	double operator() (int region, double lxy1, double lxy2) const;

	// Return a copy, and the caller will own it
	std::unique_ptr<TH2D> clone_weight(int region) const {
		if (region < 0 || region > 3) {
			throw std::runtime_error("Illegial region number");
		}
		auto a = std::unique_ptr<TH2D>(static_cast<TH2D*>(_pass_weight[region]->Clone()));
		a->SetDirectory(nullptr);
		return a;
	}

private:
	std::unique_ptr<TH2D> _pass_weight[4];
};

