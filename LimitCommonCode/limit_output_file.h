// Write out a limit file
#ifndef __limit_output_file__
#define __limit_output_file__

#include "limit_datastructures.h"

#include "TFile.h"
#include "TH1D.h"

#include <string>
#include <vector>
#include <memory>

inline void write_limit_output_file(const std::string &fname, const std::vector<limit_result> &results)
{
	auto f = std::unique_ptr<TFile>(TFile::Open(fname.c_str(), "RECREATE"));

	auto h_95 = new TH1D("CL_95", "CL_95", results.size(), 0.0, results.size());
	auto h_1sigma = new TH1D("sigma_1", "sigma_1", results.size(), 0.0, results.size());
	auto h_2sigma = new TH1D("sigma_2", "sigma_2", results.size(), 0.0, results.size());

	for (int i = 0; i < results.size(); i++) {
		h_95->SetBinContent(i + 1, results[i].cl_95 * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]);
		h_1sigma->SetBinContent(i + 1, results[i].cl_1sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]);
		h_2sigma->SetBinContent(i + 1, results[i].cl_2sigma * results[i].signal.signalEvents.A / results[i].signal.efficiency[A]);
	}

	f->Write();
}

#endif
