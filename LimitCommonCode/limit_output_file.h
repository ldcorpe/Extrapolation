// Write out a limit file
#ifndef __limit_output_file__
#define __limit_output_file__

#include "limit_datastructures.h"
#include "variable_binning_builder.h"

#include "TFile.h"
#include "TH1D.h"

#include <string>
#include <vector>
#include <memory>

inline double protected_val(double val) {
	return std::isfinite(val) ? val 
		: 1000;
}

template<typename T95, typename T1, typename T2, typename T3, typename T4, typename T5>
class group_plotter
{
public:
	group_plotter(const variable_binning_builder &binning, const std::string &cl_name, const std::string &sigma_prefix,
		T95 func95, T1 funcP1S, T2 funcN1S, T3 funcP2S, T4 funcN2S, T5 funcLimit)
		: _func95(func95), _funcP1S(funcP1S), _funcP2S(funcP2S), _funcN1S(funcN1S), _funcN2S(funcN2S), _funcLimit(funcLimit)
	{
		_h_95 = new TH1D(cl_name.c_str(), cl_name.c_str(), binning.nbin(), binning.bin_list());
		_h_p1sigma = new TH1D((sigma_prefix + "_p1").c_str(), (sigma_prefix + "_p1").c_str(), binning.nbin(), binning.bin_list());
		_h_p2sigma = new TH1D((sigma_prefix + "_p2").c_str(), (sigma_prefix + "_p2").c_str(), binning.nbin(), binning.bin_list());
		_h_n1sigma = new TH1D((sigma_prefix + "_n1").c_str(), (sigma_prefix + "_n1").c_str(), binning.nbin(), binning.bin_list());
		_h_n2sigma = new TH1D((sigma_prefix + "_n2").c_str(), (sigma_prefix + "_n2").c_str(), binning.nbin(), binning.bin_list());
		_h_limit = new TH1D((sigma_prefix + "_limit").c_str(), (sigma_prefix + "_limit").c_str(), binning.nbin(), binning.bin_list());
	}

	void operator() (const std::vector<limit_result> &results) {
		int i_bin = 1;
		for (const auto &r : results) {
			_h_95->SetBinContent(i_bin, protected_val(_func95(r)));
			_h_p1sigma->SetBinContent(i_bin, protected_val(_funcP1S(r)));
			_h_p2sigma->SetBinContent(i_bin, protected_val(_funcP2S(r)));
			_h_n1sigma->SetBinContent(i_bin, protected_val(_funcN1S(r)));
			_h_n2sigma->SetBinContent(i_bin, protected_val(_funcN2S(r)));
			_h_limit->SetBinContent(i_bin, protected_val(_funcLimit(r)));
			i_bin++;
		}
	}
private:
	TH1D *_h_95;
	TH1D *_h_p1sigma;
	TH1D *_h_p2sigma;
	TH1D *_h_n1sigma;
	TH1D *_h_n2sigma;
	TH1D *_h_limit;
	
	T95 _func95;
	T1 _funcP1S;
	T3 _funcP2S;
	T2 _funcN1S;
	T4 _funcN2S;
	T5 _funcLimit;
};

template<typename T95, typename T1, typename T2, typename T3, typename T4, typename T5>
group_plotter<T95, T1, T2, T3, T4, T5> make_group_plotter(const variable_binning_builder &binning, const std::string &cl_name, const std::string &sigma_prefix,
	T95 func,
	T1 funcP1S, T2 funcN1S, T3 funcP2S, T4 funcN2S,
	T5 funcLimit)
{
	return group_plotter<T95, T1, T2, T3, T4, T5>(binning, cl_name, sigma_prefix, func, funcP1S, funcN1S, funcP2S, funcN2S, funcLimit);
}

inline void write_limit_output_file(const abcd_limit_config &lconfig, const std::vector<limit_result> &results,
	const variable_binning_builder &binning)
{
	auto f = std::unique_ptr<TFile>(TFile::Open(lconfig.fileName.c_str(), "RECREATE"));

	double lumi = lconfig.luminosity * 1000; //pb^-1

	// Plot signal strength
	make_group_plotter(binning, "mu_95", "mu_sigma_",
		[](const limit_result &r) { return r.cl_95; },
		[](const limit_result &r) { return r.cl_p1sigma; },
		[](const limit_result &r) { return r.cl_n1sigma; },
		[](const limit_result &r) { return r.cl_p2sigma; },
		[](const limit_result &r) { return r.cl_n2sigma; },
		[](const limit_result &r) { return r.cl_limit; }
	)(results);

	// Plot number of events that are expected at this lifetime
	make_group_plotter(binning, "expected_events_95CL", "expected_allowed_events_",
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; }
	)(results);

	// Plot number of events that are allowed by the limit setting proceedure
	make_group_plotter(binning, "allowed_events_95CL", "sigma_allowed_events_",
		[](const limit_result &r) { return r.cl_95 * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_p1sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_n1sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_p2sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_n2sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_limit * r.signal.signalEvents.A; }
	)(results);
		
	// Plot x-sec X BR
	make_group_plotter(binning, "xsec_BR_95CL", "xsec_BR_events_",
		[lumi](const limit_result &r) { return r.cl_95 * r.signal.signalEvents.A / r.signal.efficiency[A] / lumi; },
		[lumi](const limit_result &r) { return r.cl_p1sigma * r.signal.signalEvents.A / r.signal.efficiency[A] / lumi; },
		[lumi](const limit_result &r) { return r.cl_n1sigma * r.signal.signalEvents.A / r.signal.efficiency[A] / lumi; },
		[lumi](const limit_result &r) { return r.cl_p2sigma * r.signal.signalEvents.A / r.signal.efficiency[A] / lumi; },
		[lumi](const limit_result &r) { return r.cl_n2sigma * r.signal.signalEvents.A / r.signal.efficiency[A] / lumi; },
		[lumi](const limit_result &r) { return r.cl_limit * r.signal.signalEvents.A / r.signal.efficiency[A] / lumi; }
	)(results);

	f->Write();
}

#endif
