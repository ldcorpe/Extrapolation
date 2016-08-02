#ifndef __SimulABCD__
#define __SimulABCD__

#include "limit_datastructures.h"

#include "TROOT.h"

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

Double_t simultaneousABCD(const Double_t n[4], const Double_t s[4], const Double_t b[4], const Double_t c[4],
	TString out_filename = "ABCD_ws.root",
	Bool_t useB = kFALSE, // Use background as estimated in MC
	Bool_t useC = kFALSE, // Use other background events (do subtraction of c above
	Bool_t blindA = kTRUE, // Assume no signal, so we get expected limits
	Int_t calcType = 0 // 0 for toys, 2 for asym fit
);

inline Double_t simultaneousABCD(const std::vector<double> &n,
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

inline limit_result do_abcd_limit(const ABCD &data, const ABCD &expected_signal, const abcd_limit_config &config);
#endif

