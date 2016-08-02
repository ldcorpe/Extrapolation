// Make access to items in an extrapolation file easy so we don't have to pepper
// our code with weird ROOT calls.

#ifndef __extrap_file_wrapper__
#define __extrap_file_wrapper__

#include "limit_datastructures.h"

#include <TFile.h>
#include <TH1D.h>

#include <string>
#include <memory>
#include <vector>

class extrap_file_wrapper {
public:
	inline extrap_file_wrapper(const std::string &filename)
		: _loaded_generated_info(false)
	{
		_file = std::unique_ptr<TFile>(TFile::Open(filename.c_str(), "READ"));
		if (!_file->IsOpen()) {
			throw std::runtime_error((std::string("Unable to locate file ") + filename + " - can't continue").c_str());
		}
	}

	// Return the generated lifetime and A, B, C, and D for this MC input.
	inline signal_lifetime generated_lifetime() const {
		if (!_loaded_generated_info) {
			_generated_lifetime.efficiency = load_histo_as_vector("eff_as_generated");
			_generated_lifetime.lifetime = load_histo_as_vector("generated_ctau")[0];
			auto abcd = load_histo_as_vector("n_passed_as_generated");
			_generated_lifetime.signalEvents.A = abcd[A];
			_generated_lifetime.signalEvents.B = abcd[B];
			_generated_lifetime.signalEvents.C = abcd[C];
			_generated_lifetime.signalEvents.D = abcd[D];

			_loaded_generated_info = true;
		}
		return _generated_lifetime;
	}

	// Return a list of the lifetimes that we have efficiencies and #'s of events for.
	// We will look at the center of each bin.
	inline std::vector<double> list_of_lifetimes() const
	{
		load_efficiencies();
		std::vector<double> result;

		for (int i = 1; i <= _efficiencies[0]->GetNbinsX(); i++) {
			result.push_back(_efficiencies[0]->GetBinCenter(i));
		}

		return result;
	}

	// Return what the signal sample woudl look like at a particular point. Do this
	// by scaling from the original generated lifetime.
	inline signal_lifetime lifetime(double ctau) const
	{
		const auto &original = generated_lifetime();
		auto result = original;
		result.lifetime = ctau;
		
		// Update the efficiencies & scale by # of evenets
		load_efficiencies();
		auto i_bin = _efficiencies[0]->FindBin(ctau);
		for (int i = 0; i < 4; i++) {
			result.efficiency[i] = _efficiencies[i]->GetBinContent(i_bin);
		}
		result.signalEvents.A *= result.efficiency[A] / original.efficiency[A];
		result.signalEvents.B *= result.efficiency[B] / original.efficiency[B];
		result.signalEvents.C *= result.efficiency[C] / original.efficiency[C];
		result.signalEvents.D *= result.efficiency[D] / original.efficiency[D];

		return result;
	}

private:
	std::unique_ptr<TFile> _file;

	mutable bool _loaded_generated_info;
	mutable signal_lifetime _generated_lifetime;

	// Fetch a histogram from the file. NOTE: Do not delete this!!
	// ROOT memory management will deal with it automatically.
	TH1D* load_histo(const std::string &name) const
	{
		auto h = static_cast<TH1D*>(_file->Get(name.c_str()));
		if (h == nullptr) {
			throw std::runtime_error(std::string("Unable to load histogram ") + name + " from extrapolation file");
		}
		return h;
	}

	// Load a histogram from the file and turn it into a vector.
	std::vector<double> load_histo_as_vector(const std::string &name) const
	{
		auto h = load_histo(name);
		std::vector<double> result;
		for (int i = 1; i <= h->GetNbinsX(); i++) {
			result.push_back(h->GetBinContent(i));
		}

		return result;
	}

	mutable std::vector<TH1D*> _efficiencies;
	void load_efficiencies() const
	{
		if (_efficiencies.size() == 0) {
			_efficiencies.push_back(load_histo("h_res_eff_A"));
			_efficiencies.push_back(load_histo("h_res_eff_B"));
			_efficiencies.push_back(load_histo("h_res_eff_C"));
			_efficiencies.push_back(load_histo("h_res_eff_D"));
		}
	}
};

#endif
