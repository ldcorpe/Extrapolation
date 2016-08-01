#ifndef __SimulABCD__
#define __SimulABCD__

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

Double_t simultaneousABCD(Double_t n[4], Double_t s[4], Double_t b[4], Double_t c[4],
	TString out_filename = "ABCD_ws.root",
	Bool_t useB = kFALSE, // Use background as estimated in MC
	Bool_t useC = kFALSE, // Use other background events (do subtraction of c above
	Bool_t blindA = kTRUE, // Assume no signal, so we get expected limits
	Int_t calcType = 0 // 0 for toys, 2 for asym fit
);

#endif

