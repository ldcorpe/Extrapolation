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
	return val > 10000 ? 10000 
		: val < 0.0 ? 10000
		: val;
}

inline void write_limit_output_file(const abcd_limit_config &lconfig, const std::vector<limit_result> &results)
{
	auto f = std::unique_ptr<TFile>(TFile::Open(lconfig.fileName.c_str(), "RECREATE"));

	auto h_95 = new TH1D("CL_95", "CL_95", results.size(), 0.0, results.size());
	auto h_p1sigma = new TH1D("sigma_p1", "sigma_p1", results.size(), 0.0, results.size());
	auto h_p2sigma = new TH1D("sigma_p2", "sigma_p2", results.size(), 0.0, results.size());
	auto h_n1sigma = new TH1D("sigma_n1", "sigma_n1", results.size(), 0.0, results.size());
	auto h_n2sigma = new TH1D("sigma_n2", "sigma_n2", results.size(), 0.0, results.size());

	for (int i = 0; i < results.size(); i++) {
		h_95->SetBinContent(i + 1, protected_val (results[i].cl_95 * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
		h_p1sigma->SetBinContent(i + 1, protected_val(results[i].cl_p1sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
		h_p2sigma->SetBinContent(i + 1, protected_val(results[i].cl_p2sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
		h_n1sigma->SetBinContent(i + 1, protected_val(results[i].cl_n1sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
		h_n2sigma->SetBinContent(i + 1, protected_val(results[i].cl_n2sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]));
	}

	f->Write();
}

#endif
