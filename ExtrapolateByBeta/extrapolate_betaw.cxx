// Given a limit for a sample and an input file, extrapolate to a range of lifetimes.
//
// This code was lifted (and modified) from the Run 1 CalRatio analysis.
//
// Main goal: Calculate the relative efficiency of the analysis as a function of proper life-time.
//            - It calculates the absolute efficiency.
//
// Inputs are taken on the command line (to make it easy to use in a tool chain).
//   arg1: ROOT File that contains the MuonTree.
//   arg2: proper lifetime of generation of the sample we are looking at (e.g. 5 or a ctau = 5 meters).

#include "doubleError.h"
#include "muon_tree_processor.h"
#include "variable_binning_builder.h"
#include "Lxy_weight_calculator.h"
#include "beta_cache.h"

#include "TApplication.h"
#pragma warning (push)
#pragma warning (disable: 4244)
#include "TLorentzVector.h"
#pragma warning (pop)
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <memory>

using namespace std;

// Some config constants

// How many loops in tau should we do?
size_t n_tau_loops = 5;

// Helper methods
struct extrapolate_config {
	string _muon_tree_root_file;
	string _output_filename;
	double _tau_gen;
};
extrapolate_config parse_command_line(int argc, char **argv);
variable_binning_builder PopulateTauTable();
pair<vector<unique_ptr<TH2F>>, unique_ptr<TH2F>> GetFullBetaShape(double tau, int ntauloops, const muon_tree_processor &mc_entries, const Lxy_weight_calculator &lxyWeight);
template<class T> vector<unique_ptr<T>> DivideShape(
	const pair<vector<unique_ptr<T>>, unique_ptr<T>> &r,
	const string &name, const string &title);
vector<doubleError> CalcPassedEvents(const muon_tree_processor &reader, const vector<unique_ptr<TH2F>> &weightHist, bool eventCountOnly = false);
std::pair<Double_t, Double_t> getBayes(const doubleError &num, const doubleError &den);
void SetAsymError(unique_ptr<TGraphAsymmErrors> &g, int bin, double tau, double bvalue, const pair<double, double> &assErrors);

unique_ptr<TH1D> save_as_histo(const string &name, double number);
unique_ptr<TH1D> save_as_histo(const string &name, const vector<double> &number);
unique_ptr<TH1D> save_as_histo(const string &name, const vector<doubleError> &number);

////////////////////////////////
// Main entry point.
////////////////////////////////
int main(int argc, char**argv)
{
	int dummy_argc = 0;
	auto a = make_unique<TApplication>("extrapolate_betaw", &dummy_argc, argv);
	try {
		// Pull out the various command line arguments
		auto config = parse_command_line(argc, argv);

		// Create the muon tree reader object and calculate the lxy weighting histogram
		muon_tree_processor reader (config._muon_tree_root_file);
		Lxy_weight_calculator lxy_weight(reader);

		// How often, for the generated sample, a pair of beta1, beta2 vpions reaches the HCal.
		// This is done at generation lifetime, so this will be the baseline which we scale against
		// in the tau loop below.
		auto r = GetFullBetaShape(config._tau_gen, n_tau_loops, reader, lxy_weight);
		auto h_gen_ratio = DivideShape(r, "h_Ngen_ratio", "Fraction of events in beta space at raw generated ctau");
		auto passedEventsAtGen = CalcPassedEvents(reader, vector<unique_ptr<TH2F>> (), false);
		auto totalEventsAtGen = CalcPassedEvents(reader, vector<unique_ptr<TH2F>>(), true);

		// Create the histograms we will use to store the raw results.
		auto tau_binning = PopulateTauTable();
		vector<TH1F*> h_res_eff;
		vector<unique_ptr<TGraphAsymmErrors>> g_res_eff;

		for (int i_region = 0; i_region < 4; i_region++) {
			//Efficiency VS lifetime
			ostringstream name_h;
			name_h << "h_res_eff_" << (char) ('A' + i_region);
			h_res_eff.push_back(new TH1F(name_h.str().c_str(), name_h.str().c_str(), tau_binning.nbin(), tau_binning.bin_list()));

			g_res_eff.push_back(make_unique<TGraphAsymmErrors>(tau_binning.nbin()));
			ostringstream name_g;
			name_g << "g_res_eff_" << (char) ('A' + i_region);
			g_res_eff[i_region]->SetName(name_g.str().c_str());
			ostringstream title_g;
			title_g << "Absolute efficiency to the generated analysis for region " << (char) ('A' + i_region);
			g_res_eff[i_region]->SetTitle(title_g.str().c_str());
		}

		// Loop over proper lifetime
		for (unsigned int i_tau = 0; i_tau < tau_binning.nbin(); i_tau++) {
			auto tau = h_res_eff[0]->GetBinCenter(i_tau+1); // Recal ROOT indicies bins at 1

			// Get the full Beta shape
			auto rtau = GetFullBetaShape(tau, n_tau_loops, reader, lxy_weight);
			auto h_caut_ratio = DivideShape(rtau, "h_caut_ratio", "h_caut_ratio");

			// Now, create a weighting histogram. This is just the differece between the numerators at the
			// extrapolated ctau and at the generated ctau
			decltype(h_caut_ratio) h_Nratio;
			for (int i_region = 0; i_region < 4; i_region++) {
				auto h = unique_ptr<TH2F>(static_cast<TH2F*>(h_caut_ratio[i_region]->Clone()));
				h->Divide(h_gen_ratio[i_region].get());
				h_Nratio.push_back(move(h));
			}

			// The the number of events that passed for this lifetime.
			auto passedEventsAtTau = CalcPassedEvents(reader, h_Nratio, false);

			// Calculate proper asymmetric errors and save the extrapolation result for the change in efficency.
			for (int i_region = 0; i_region < 4; i_region++) {
				doubleError eff = passedEventsAtTau[i_region] / totalEventsAtGen[i_region];
				std::pair<Double_t, Double_t> bayerr_sig_perc = getBayes(passedEventsAtTau[i_region], totalEventsAtGen[i_region]);
				Double_t erro = (bayerr_sig_perc.first + bayerr_sig_perc.second)*0.5;

				h_res_eff[i_region]->SetBinContent(i_tau + 1, eff.value());
				h_res_eff[i_region]->SetBinError(i_tau + 1, erro);
				SetAsymError(g_res_eff[i_region], i_tau, tau, eff.value(), bayerr_sig_perc);
			}

			cout << " tau = " << tau << " npassed = " << passedEventsAtTau[0] << " passed tau/gen = " << passedEventsAtTau[0] / passedEventsAtGen[0] << " global eff = " << passedEventsAtTau[0] / totalEventsAtGen[0] << endl;
		}

		// Save plots in the output file
		auto output_file = unique_ptr<TFile>(TFile::Open(config._output_filename.c_str(), "RECREATE"));
		for (int i_region = 0; i_region < 4; i_region++) {
			output_file->Add(h_res_eff[i_region]);
		}

		// Save basic information for the generated sample.
		output_file->Add(save_as_histo("generated_ctau", config._tau_gen).release());
		output_file->Add(save_as_histo("n_passed_as_generated", passedEventsAtGen).release());
		output_file->Add(save_as_histo("n_as_generated", totalEventsAtGen).release());

		vector<doubleError> effAtGen;
		for (int i = 0; i < 4; i++) {
			effAtGen.push_back(passedEventsAtGen[i] / totalEventsAtGen[i]);
		}
		output_file->Add(save_as_histo("eff_as_generated", effAtGen).release());

		output_file->Write();
	}
	catch (exception &e)
	{
		cout << "Failed!!" << endl;
		cout << " --> " << e.what() << endl;
		return 1;
	}
	return 0;
}

// Parse the command line.
extrapolate_config parse_command_line(int argc, char **argv)
{
	if (argc != 4) {
		cout << "Usage: extrapolate_betaw <muonTree-file> <output-filename> <generated-ctau>" << endl;
		cout << endl;
		cout << "    <muonTree-file>           TFile containing the muonTree root file generated from a particular MC sample." << endl;
		cout << "    <output-filename>		   ROOT filename for output" << endl;
		cout << "    <generated-ctau>          The lifetime (in meters) where the sample was generated." << endl;

		throw runtime_error("Wrong number of arguments");
	}

	extrapolate_config r;
	r._muon_tree_root_file = argv[1];
	r._output_filename = argv[2];
	r._tau_gen = atof(argv[3]);
	return r;
}

// Initalize and populate the tau decay table.
// Due to the fact we run out of stats, this is, by its very nature, not equal binning.
variable_binning_builder PopulateTauTable()
{
	variable_binning_builder r(0.0);
#ifdef TEST_RUN
	r.bin_up_to(1.0, 0.05);
#else
	r.bin_up_to(0.6, 0.005);
	r.bin_up_to(4.0, 0.05);
	r.bin_up_to(50.0, 0.2);
#endif
	return r;
}

// Binning we will use for beta histograms
variable_binning_builder PopulateBetaBinning()
{
	variable_binning_builder beta_binning(0.0);
	beta_binning.bin_up_to(0.8, 0.2);
	beta_binning.bin_up_to(1.0, 0.05);
	return beta_binning;
}

// Sample from the proper lifetime tau for a specific lifetime, and then do the special relativity
// calculation to understand where it ended up.
void doSR(TLorentzVector vpi1, TLorentzVector vpi2, Double_t tau, Double_t &beta1, Double_t &beta2, Double_t &L2D1, Double_t &L2D2) {

	beta1 = vpi1.Beta();
	beta2 = vpi2.Beta();
	Double_t gamma1 = vpi1.Gamma();
	Double_t gamma2 = vpi2.Gamma();

	Double_t ct1 = gRandom->Exp(tau);  //get ct of the two vpions 
	Double_t ct2 = gRandom->Exp(tau);

	Double_t ct1prime = gamma1 * ct1;
	Double_t ct2prime = gamma2 * ct2;
	Double_t lxy1 = beta1 * ct1prime;  //construct Lxy of the two vpions 
	Double_t lxy2 = beta2 * ct2prime;

	Double_t lighttime1 = lxy1 / 2.9979E8 * 1E9;
	Double_t lighttime2 = lxy2 / 2.9979E8 * 1E9;
	Double_t timing1 = ct1prime / 2.9979E8 * 1E9;
	Double_t timing2 = ct2prime / 2.9979E8 * 1E9;
	Double_t deltat1 = timing1 - lighttime1;
	Double_t deltat2 = timing2 - lighttime2;

	Double_t theta1 = vpi1.Theta();
	Double_t theta2 = vpi2.Theta();

	TVector3 vpixyz1, vpixyz2;
	vpixyz1.SetMagThetaPhi(lxy1, theta1, vpi1.Phi());
	vpixyz2.SetMagThetaPhi(lxy2, theta2, vpi2.Phi());

	L2D1 = vpixyz1.Perp();
	L2D2 = vpixyz2.Perp();

	return;
}

// Generate a lxy1 and lxy2 set of histograms. The denominator (first item) is
// the generated ones. The numerator is modified by the acceptance histogram
// built from the input file.
pair<vector<unique_ptr<TH2F>>, unique_ptr<TH2F>> GetFullBetaShape(double tau, int ntauloops, const muon_tree_processor &mc_entries, const Lxy_weight_calculator &lxyWeight)
{
	// Create numerator and denominator histograms.
	// To avoid annoying ROOT error messages, make a unique name for each.
	ostringstream dname, nname;
	dname << "tau_" << tau << "_den";
	nname << "tau_" << tau << "_num";
	auto beta_binning(PopulateBetaBinning());
	unique_ptr<TH2F> den(make_unique<TH2F>(dname.str().c_str(), dname.str().c_str(), beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list()));
	vector<unique_ptr<TH2F>> num;
	for (int i_region = 0; i_region < 4; i_region++) {
		nname << "A";
		num.push_back(make_unique<TH2F>(nname.str().c_str(), nname.str().c_str(), beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list()));
		num[i_region]->Sumw2();
	}
	den->Sumw2();

	// Loop over each MC entry, and generate tau's at several different places
	mc_entries.process_all_entries([&den, &num, ntauloops, tau, &lxyWeight](const muon_tree_processor::eventInfo &entry) {
		TLorentzVector vpi1, vpi2;
		vpi1.SetPtEtaPhiE(entry.vpi1_pt, entry.vpi1_eta, entry.vpi1_phi, entry.vpi1_E);
		vpi2.SetPtEtaPhiE(entry.vpi2_pt, entry.vpi2_eta, entry.vpi2_phi, entry.vpi2_E);

		for (Int_t maketaus = 0; maketaus < ntauloops; maketaus++) { // tau loop to generate toy events

			Double_t beta1 = -1, beta2 = -1, L2D1 = -1, L2D2 = -1;

			doSR(vpi1, vpi2, tau, beta1, beta2, L2D1, L2D2); // Do special relativity

			den->Fill(beta1, beta2, entry.weight);
			for (int i_region = 0; i_region < 4; i_region++) {
				num[i_region]->Fill(beta1, beta2, entry.weight * lxyWeight(i_region, L2D1, L2D2));
			}
		}
	});

	return make_pair(move(num), move(den));
}

// Calculate a 2D efficiency given a denominator and the numerator selected from the denominator
template<class T>
vector<unique_ptr<T>> DivideShape(const pair<vector<unique_ptr<T>>, unique_ptr<T>> &r, const string &name, const string &title)
{
	vector<unique_ptr<T>> result;
	for (auto &info : r.first) {
		unique_ptr<T> ratio(static_cast<T*>(info->Clone()));
		ratio->Divide(info.get(), r.second.get(), 1.0, 1.0, "B");
		ratio->SetNameTitle(name.c_str(), title.c_str());
		result.push_back(move(ratio));
	}

	return result;
}

// Calculate the number of events that pass our cuts (possibly weighted).
vector<doubleError> CalcPassedEvents(const muon_tree_processor &reader, const vector<unique_ptr<TH2F>> &weightHist, bool eventCountOnly)
{
	vector<doubleError> nEvents(4); // This will be the number of events in the TTree (without gluons)

	// Calculate the event weight. A combination of the pile up reweighting from the ntuple and perhaps
	// the beta re-weighting from the input histogram.
	reader.process_all_entries([&weightHist, eventCountOnly, &nEvents](const muon_tree_processor::eventInfo &entry) {

		for (int i_region = 0; i_region < 4; i_region++) {

			// TODO: is this the right way to do an error here? Does it propagate so do we care?
			doubleError weight(entry.weight, entry.weight);

			beta_cache b(entry);
			if (weightHist.size() > 0) {
				int nbin = weightHist[i_region]->FindBin(b.beta1(), b.beta2());
				weight *= doubleError(weightHist[i_region]->GetBinContent(nbin), weightHist[i_region]->GetBinError(nbin));
			}

			// Count the event and populate the output histogram, if
			// We are doing an event count only (e.g. the denominator) or
			// it passes our analysis cuts (e.g. the numerator).
			// TODO: when eventCountOnly is true, every single event goes into each region, no matter what.
			//       make sure that is what we want.
			auto inRegion =
				i_region == 0 ? entry.RegionA
				: i_region == 1 ? entry.RegionB
				: i_region == 2 ? entry.RegionC
				: entry.RegionD;
			if (eventCountOnly
				|| inRegion
				) {

				nEvents[i_region] += weight;
			}
		}
	});

	return nEvents;
}

// Calc error via bayes
std::pair<Double_t, Double_t> getBayes(const doubleError &num, const doubleError &den) {
	// returns the Bayesian uncertainty over the num/den ratio
	std::pair<Double_t, Double_t> result(0., 0.);

	auto h_num = make_unique<TH1D>("h_num", "", 1, 0, 1);
	auto h_den = make_unique<TH1D>("h_den", "", 1, 0, 1);

	h_num->SetBinContent(1, num.value());
	h_den->SetBinContent(1, den.value());
	h_num->SetBinError(1, num.err());
	h_den->SetBinError(1, den.err());

	auto h_eff = make_unique<TGraphAsymmErrors>();
	h_eff->BayesDivide(h_num.get(), h_den.get());//, "w"); 

	result.first = h_eff->GetErrorYlow(0);
	result.second = h_eff->GetErrorYhigh(0);

	return result;
}

// Set the asymmetric error simply
void SetAsymError(unique_ptr<TGraphAsymmErrors> &g, int bin, double tau, double bvalue, const pair<double, double> &assErrors)
{
	g->SetPoint(bin, tau, bvalue);
	g->SetPointEYlow(bin, assErrors.first);
	g->SetPointEYhigh(bin, assErrors.second);
}

// Helper func to save a number into a histogram.
unique_ptr<TH1D> save_as_histo(const string &name, double number)
{
	auto h = make_unique<TH1D>(name.c_str(), name.c_str(), 1, 0.0, 1.0);
	h->SetDirectory(nullptr);
	h->SetBinContent(1, number);
	return move(h);
}

// Save a vector into a histogram.
unique_ptr<TH1D> save_as_histo(const string &name, const vector<double> &number)
{
	auto h = make_unique<TH1D>(name.c_str(), name.c_str(), number.size(), 0.0, 1.0);
	h->SetDirectory(nullptr);
	for (int i = 0; i < number.size(); i++) {
		h->SetBinContent(i + 1, number[i]);
	}
	return move(h);
}

// Save a vector into a histogram.
unique_ptr<TH1D> save_as_histo(const string &name, const vector<doubleError> &number)
{
	auto h = make_unique<TH1D>(name.c_str(), name.c_str(), number.size(), 0.0, 1.0);
	h->SetDirectory(nullptr);
	for (int i = 0; i < number.size(); i++) {
		h->SetBinContent(i + 1, number[i].value());
		h->SetBinError(i + 1, number[i].err());
	}
	return move(h);
}
