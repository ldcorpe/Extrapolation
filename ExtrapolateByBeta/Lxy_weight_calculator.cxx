#include "Lxy_weight_calculator.h"

using namespace std;

// Defien the constants for the Lxy1 and Lxy2 binning
const int n_bins = 100; // 100 bins
const double lxy_min = 0.0; // lxy minimum in meters
const double lxy_max = 5.0; // lxy maximum in meters


//
// Calcualte the weight histogram
// Note: In Run1 this was done during the actual running on MC that generates this ntuple. But by doing it here
//       we better keep cuts consistent by keeping them close to each other.
Lxy_weight_calculator::Lxy_weight_calculator(const muon_tree_processor &reader)
{
	// Events where they were generated and where they passed the analysis selection
	unique_ptr<TH2D> generated(new TH2D("generated", "generated", n_bins, lxy_min, lxy_max, n_bins, lxy_min, lxy_max));
	unique_ptr<TH2D> passed(new TH2D("passed", "passed", n_bins, lxy_min, lxy_max, n_bins, lxy_min, lxy_max));
	generated->Sumw2();
	passed->Sumw2();

	reader.process_all_entries([&generated, &passed](const muon_tree_processor::eventInfo &entry) {
		generated->Fill(entry.vpi1_Lxy/1000.0, entry.vpi2_Lxy/1000.0);
		if (entry.RegionA) {
			passed->Fill(entry.vpi1_Lxy/1000.0, entry.vpi2_Lxy/1000.0);
		}
	});

	// The key is the ratio.
	// We can't use ROOT sumw2 errors because they assume independent histograms. So we use
	// binomial errors (the "B" option).
	_pass_weight = unique_ptr<TH2D>(new TH2D("_pass_weight", "_pass_weight", n_bins, lxy_min, lxy_max, n_bins, lxy_min, lxy_max));
	_pass_weight->Divide(passed.get(), generated.get(), 1.0, 1.0, "B");
}


Lxy_weight_calculator::~Lxy_weight_calculator()
{
}

// Do the lookup to make the calculation
double Lxy_weight_calculator::operator()(double lxy1, double lxy2) const
{
	int xbin = _pass_weight->FindBin(lxy1, lxy2);
	return _pass_weight->GetBinContent(xbin);
}
