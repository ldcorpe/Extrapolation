#ifndef __SimulABCD__
#define __SimulABCD__

#include "limit_datastructures.h"
#include "HypoTestInvTool.h"

#include "TROOT.h"

#include <iostream>

/*
Two variables x,y create 4 regions: A,B,C,D with A signal dominated region and
B,C,D control regions (sidebands).

^ y (sumpt)
|
|--------------------------+
|    C         |     B     |
|              |           |
|              |           |
|              |           |
|              |           |
|--------------+-----------|
|    D         |     A     |
|              |           |
|              |  Signal   |
|              |           |
+-----------------------------> x (deltaPhi)
^

ABCD Ansatz A:B = D:C -->  A = B * D/C

*/

HypoTestInvTool::LimitResults simultaneousABCD(const Double_t n[4], const Double_t s[4], const Double_t b[4], const Double_t c[4],
	TString out_filename = "ABCD_ws.root",
	Bool_t useB = kFALSE, // Use background as estimated in MC
	Bool_t useC = kFALSE, // Use other background events (do subtraction of c above
	Bool_t blindA = kTRUE, // Assume no signal, so we get expected limits
	Int_t calcType = 0 // 0 for toys, 2 for asym fit
);

inline HypoTestInvTool::LimitResults simultaneousABCD(const std::vector<double> &n,
	const std::vector<double> &s,
	const std::vector<double> &b,
	const std::vector<double> &c,
	TString out_filename = "ABCD_ws.root",
	Bool_t useB = kFALSE, // Use background as estimated in MC
	Bool_t useC = kFALSE, // Use other background events (do subtraction of c above
	Bool_t blindA = kTRUE, // Assume no signal, so we get expected limits
	Int_t calcType = 0 // 0 for toys, 2 for asym fit
)
{
	if (n.size() != 4
		|| s.size() != 4
		|| b.size() != 4
		|| c.size() != 4) {
		throw std::runtime_error("n, s, b, and c must have 4 elements exactly");
	}

	return simultaneousABCD(&(n[0]), &(s[0]), &(b[0]), &(c[0]),
		out_filename, useB, useC, blindA, calcType);
}

inline std::vector<double> ABCD_as_vector(const ABCD &events)
{
	std::vector<double> n;
	n.push_back(events.A);
	n.push_back(events.B);
	n.push_back(events.C);
	n.push_back(events.D);
	return n;
}

inline limit_result do_abcd_limit(const ABCD &data, const signal_lifetime &expected_signal, const abcd_limit_config &config)
{
	std::vector<double> dummy(4);

	auto limit = simultaneousABCD(ABCD_as_vector(data), ABCD_as_vector(expected_signal.signalEvents),
		dummy, dummy,
		"limit_calc.root",
		false, false,
		data.A == 0,
		config.useToys ? 0 : 2);

	std::cout << "Limit. data: " << data << "  expected signal: " << expected_signal << std::endl;
	std::cout << "  -> " << limit << std::endl;

	limit_result r;
	r.cl_p1sigma = limit.sigma_plus_1;
	r.cl_p2sigma = limit.sigma_plus_2;
	r.cl_n1sigma = limit.sigma_minus_1;
	r.cl_n2sigma = limit.sigma_minus_2;
	r.cl_95 = limit.median;
	r.signal = expected_signal;
	r.observed_data = data;
	return r;
}
#endif

