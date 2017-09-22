#pragma once

#include <vector>
#include <TH1D.h>

// Helper class to make making variable binning easy
class variable_binning_builder {
public:
	// Extract the binning from the histogram's x-axis.
	variable_binning_builder(const TH1D *h)
	{
		auto axis = h->GetXaxis();
		for (int i = 1; i <= axis->GetNbins(); i++) {
			_v.push_back(axis->GetBinLowEdge(i));
		}
		_v.push_back(axis->GetBinUpEdge(axis->GetNbins()));
	}

	// Start from a point.
	explicit variable_binning_builder(double low)
	{
		_v.push_back(low);
	}

	// Add number bins, teach bin_width apart.
	void add_bins(int number, double bin_width)
	{
		for (int i = 0; i < number; i++) {
			_v.push_back(*(_v.end() - 1) + bin_width);
		}
	}
	void add_binning_interval(double interval, double bin_width)
	{
		for (double count = 0.0; count <= bin_width; count++)
		{
			_v.push_back(*(_v.end() - 1) + bin_width);
		}
	}
	void bin_up_to(double max_value, double bin_width)
	{
		while (*(_v.end() - 1) < max_value) {
			_v.push_back(*(_v.end() - 1) + bin_width);
		}
	}
	// Number of bins to delcare to ROOT - note it is one less, as we have to
	// give it the upper adn lower boundaries of the last and first bin!
	int nbin() const { return _v.size() - 1; }
	double *bin_list() const { return (double*)&(_v[0]); }

private:
	std::vector<double> _v;
};
