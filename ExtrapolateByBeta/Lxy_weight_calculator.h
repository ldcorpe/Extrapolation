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

	double operator() (double lxy1, double lxy2) const;

	// Return a copy, and the caller will own it
	std::unique_ptr<TH2D> clone_weight() const {
		return std::unique_ptr<TH2D>(static_cast<TH2D*>(_pass_weight->Clone()));
	}

private:
	std::unique_ptr<TH2D> _pass_weight;
};

