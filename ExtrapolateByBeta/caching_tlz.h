#ifndef __caching_tlz__
#define __caching_tlz__
// Caching TLorentzVector. Assumes you never change or update it (those methods aren't implemented!).
// To protect against this, this isn't a TLZ.

#include "TLorentzVector.h"
#include "cache_object.h"

class caching_tlz {
public:
	caching_tlz(const TLorentzVector &v)
		: _beta(make_cache_func<double>([&v]() {return v.Beta(); })),
		_gamma(make_cache_func<double>([&v]() {return v.Gamma(); })),
		_theta(make_cache_func<double>([&v]() {return v.Theta(); })),
		_phi(make_cache_func<double>([&v]() {return v.Phi(); })),
		_beta_transverse(make_cache_func<double>([&v]() {auto g = v.Pt() / v.M(); return sqrt(1 - 1/(g*g)); }))
	{
	}

	double Beta(void) const { return _beta.value(); }
	double Gamma(void) const { return _gamma.value(); }
	double Phi(void) const { return _phi.value(); }
	double Theta(void) const { return _theta.value(); }
	double BetaTransverse(void) const { return _beta_transverse.value(); }

private:
	cache_object<double> _beta;
	cache_object<double> _gamma;
	cache_object<double> _theta;
	cache_object<double> _phi;
	cache_object<double> _beta_transverse;
};

#endif
