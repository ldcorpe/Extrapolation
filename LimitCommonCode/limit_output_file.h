// Write out a limit file
#ifndef __limit_output_file__
#define __limit_output_file__

#include "limit_datastructures.h"

#include "TFile.h"
#include "TH1D.h"

#include <string>
#include <vector>
#include <memory>

inline double protected_val(double val) {
	return std::isfinite(val) ? val 
		: 1000;
}

template<typename T95, typename T1, typename T2, typename T3, typename T4>
class group_plotter
{
public:
	group_plotter(int nbins, const std::string &cl_name, const std::string &sigma_prefix,
		T95 func95, T1 funcP1S, T2 funcN1S, T3 funcP2S, T4 funcN2S)
		: _func95(func95), _funcP1S(funcP1S), _funcP2S(funcP2S), _funcN1S(funcN1S), _funcN2S(funcN2S)
	{
		_h_95 = new TH1D(cl_name.c_str(), cl_name.c_str(), nbins, 0.0, nbins);
		_h_p1sigma = new TH1D((sigma_prefix + "_p1").c_str(), (sigma_prefix + "_p1").c_str(), nbins, 0.0, nbins);
		_h_p2sigma = new TH1D((sigma_prefix + "_p2").c_str(), (sigma_prefix + "_p2").c_str(), nbins, 0.0, nbins);
		_h_n1sigma = new TH1D((sigma_prefix + "_p3").c_str(), (sigma_prefix + "_p3").c_str(), nbins, 0.0, nbins);
		_h_n2sigma = new TH1D((sigma_prefix + "_p4").c_str(), (sigma_prefix + "_p4").c_str(), nbins, 0.0, nbins);
	}

	void operator() (const std::vector<limit_result> &results) {
		int i_bin = 1;
		for (const auto &r : results) {
			_h_95->SetBinContent(i_bin, protected_val(_func95(r)));
			_h_p1sigma->SetBinContent(i_bin, protected_val(_funcP1S(r)));
			_h_p2sigma->SetBinContent(i_bin, protected_val(_funcP2S(r)));
			_h_n1sigma->SetBinContent(i_bin, protected_val(_funcN1S(r)));
			_h_n2sigma->SetBinContent(i_bin, protected_val(_funcN2S(r)));
			i_bin++;
		}
	}
private:
	TH1D *_h_95;
	TH1D *_h_p1sigma;
	TH1D *_h_p2sigma;
	TH1D *_h_n1sigma;
	TH1D *_h_n2sigma;
	
	T95 _func95;
	T1 _funcP1S;
	T3 _funcP2S;
	T2 _funcN1S;
	T4 _funcN2S;
};

template<typename T95, typename T1, typename T2, typename T3, typename T4>
group_plotter<T95, T1, T2, T3, T4> make_group_plotter(int nbins, const std::string &cl_name, const std::string &sigma_prefix,
	T95 func,
	T1 funcP1S, T2 funcN1S, T3 funcP2S, T4 funcN2S)
{
	return group_plotter<T95, T1, T2, T3, T4>(nbins, cl_name, sigma_prefix, func, funcP1S, funcN1S, funcP2S, funcN2S);
}

inline void write_limit_output_file(const abcd_limit_config &lconfig, const std::vector<limit_result> &results)
{
	auto f = std::unique_ptr<TFile>(TFile::Open(lconfig.fileName.c_str(), "RECREATE"));

	// Plot signal strength
	make_group_plotter(results.size(), "mu_95", "mu_sigma_",
		[](const limit_result &r) { return r.cl_95; },
		[](const limit_result &r) { return r.cl_p1sigma; },
		[](const limit_result &r) { return r.cl_n1sigma; },
		[](const limit_result &r) { return r.cl_p2sigma; },
		[](const limit_result &r) { return r.cl_n2sigma; }
	)(results);

	// Plot number of events that are expected
	make_group_plotter(results.size(), "expected_events_95CL", "expected_allowed_events_",
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.signal.signalEvents.A; }
	)(results);

	// Plot number of events that are allowed by the limit setting proceedure
	make_group_plotter(results.size(), "allowed_events_95CL", "sigma_allowed_events_",
		[](const limit_result &r) { return r.cl_95 * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_p1sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_n1sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_p2sigma * r.signal.signalEvents.A; },
		[](const limit_result &r) { return r.cl_n2sigma * r.signal.signalEvents.A; }
		)(results);
		
		//auto h_95 = new TH1D("CL_95", "CL_95", results.size(), 0.0, results.size());
	//auto h_p1sigma = new TH1D("sigma_p1", "sigma_p1", results.size(), 0.0, results.size());
	//auto h_p2sigma = new TH1D("sigma_p2", "sigma_p2", results.size(), 0.0, results.size());
	//auto h_n1sigma = new TH1D("sigma_n1", "sigma_n1", results.size(), 0.0, results.size());
	//auto h_n2sigma = new TH1D("sigma_n2", "sigma_n2", results.size(), 0.0, results.size());

	//for (int i = 0; i < results.size(); i++) {
	//	h_95->SetBinContent(i + 1, protected_val (results[i].cl_95 * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
	//	h_p1sigma->SetBinContent(i + 1, protected_val(results[i].cl_p1sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
	//	h_p2sigma->SetBinContent(i + 1, protected_val(results[i].cl_p2sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
	//	h_n1sigma->SetBinContent(i + 1, protected_val(results[i].cl_n1sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
	//	h_n2sigma->SetBinContent(i + 1, protected_val(results[i].cl_n2sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
	//}

	f->Write();
}

#endif
