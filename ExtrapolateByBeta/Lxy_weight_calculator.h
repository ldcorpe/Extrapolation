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
	inline virtual ~Lxy_weight_calculator() {};

	// Region can be 0 == A, 1 == B, 2 == C, 3 == D
	virtual double operator() (int region, double lxy1, double lxy2) const = 0;

	// Return a copy, and the caller will own it
	virtual std::unique_ptr<TH1> clone_weight(int region) const = 0;

private:
};


class Lxy_weight_calculator1D : public Lxy_weight_calculator
{
public:
	// Look up and calculate the weight table
	Lxy_weight_calculator1D(const muon_tree_processor &reader);
	~Lxy_weight_calculator1D();

	// Region can be 0 == A, 1 == B, 2 == C, 3 == D
	double operator() (int region, double lxy1, double lxy2) const;

	// Return a copy, and the caller will own it
	std::unique_ptr<TH1> clone_weight(int region) const {
		if (region < 0 || region > 3) {
			throw std::runtime_error("Illegial region number");
		}
		auto a = std::unique_ptr<TH1>(static_cast<TH1*>(_pass_weight[region]->Clone()));
		a->SetDirectory(nullptr);
		return a;
	}

private:
	std::unique_ptr<TH1D> _pass_weight[4];
};


class Lxy_weight_calculator2D : public Lxy_weight_calculator
{
public:
	// Look up and calculate the weight table
	Lxy_weight_calculator2D(const muon_tree_processor &reader);
	~Lxy_weight_calculator2D();

	// Region can be 0 == A, 1 == B, 2 == C, 3 == D
	double operator() (int region, double lxy1, double lxy2) const;

	// Return a copy, and the caller will own it
	std::unique_ptr<TH1> clone_weight(int region) const {
		if (region < 0 || region > 3) {
			throw std::runtime_error("Illegial region number");
		}
		auto a = std::unique_ptr<TH1>(static_cast<TH1*>(_pass_weight[region]->Clone()));
		a->SetDirectory(nullptr);
		return a;
	}

private:
	std::unique_ptr<TH2D> _pass_weight[4];
};

