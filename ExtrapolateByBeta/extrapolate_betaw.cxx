// Given a limit for a sample and an input file, extrapolate to a range of lifetimes.
//
// This code was lifted (and modified) from the Run 1 CalRatio analysis.
//
// Main goal: Calculate the relative efficiency of the analysis as a function of proper life-time.
//            - It calculates the absolute efficiency.
//
// Inputs are taken on the command line (to make it easy to use in a tool chain).
//   arg1: ROOT File that contains the MuonTree.
//   arg2: proper lifetime of generation of the sample we are looking at.

#include "doubleError.h"
#include "muon_tree_processor.h"
#include "variable_binning_builder.h"
#include "Lxy_weight_calculator.h"
#include "beta_cache.h"
#include "caching_tlz.h"

#include "Wild/CommandLine.h"

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
using namespace Wild::CommandLine;

// Some config constants

// How many loops in tau should we do? This is a bit dynamic
#ifdef TEST_RUN
double n_tau_loops_at_gen = 200;
size_t tau_loops(double ctau) {
	return n_tau_loops_at_gen;
}
#else
double n_tau_loops_at_gen = 200;
size_t tau_loops(double ctau) {
	return n_tau_loops_at_gen;
}
#endif
// For the study for the number of loops, see the logbook. But this will affect if the extrap
// at each lifetime stablieses, so change it with care.

enum BetaShapeType
{
	FromMC,		// Uses the MC input files to derive the MC beta shape in each region
	Unity		// Does weighting by only Lxy
};

// Helper methods
struct extrapolate_config {
	string _muon_tree_root_file;
	string _output_filename;
	double _tau_gen;
	BetaShapeType _beta_type;
};
extrapolate_config parse_command_line(int argc, char **argv);
variable_binning_builder PopulateTauTable();
pair<vector<unique_ptr<TH2F>>, unique_ptr<TH2F>> GetFullBetaShape(double tau, int ntauloops, const muon_tree_processor &mc_entries, const Lxy_weight_calculator &lxyWeight);
template<class T> vector<unique_ptr<T>> DivideShape(
	const pair<vector<unique_ptr<T>>, unique_ptr<T>> &r,
	const string &name, const string &title);
vector<doubleError> CalcPassedEvents(const muon_tree_processor &reader, const vector<unique_ptr<TH2F>> &weightHist, bool eventCountOnly = false);
vector<doubleError> CalcPassedEventsLxy(const muon_tree_processor &reader, double lifetime, Lxy_weight_calculator &lxyWeight);
vector<doubleError> GenericCalcPassedEvents(const muon_tree_processor &reader, bool ignore_preselection = false);
std::pair<Double_t, Double_t> getBayes(const doubleError &num, const doubleError &den);
bool doMCPreselection(const muon_tree_processor::eventInfo &entry);
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
		cout << "TTree input: " << config._muon_tree_root_file << endl;
		cout << "Output file: " << config._output_filename << endl;
		cout << "Generated lifetime: " << config._tau_gen << endl;
		cout << "We are using beta-based extrapolation: " << (config._beta_type == BetaShapeType::FromMC ? "yes" : "no") << endl;

		// Create the muon tree reader object.
		muon_tree_processor reader (config._muon_tree_root_file);
		reader.add_preselection([](const muon_tree_processor::eventInfo &e) {return doMCPreselection(e); });

		// The Lxy efficiency histogram has to be calculated next. We do it with a "denominator", to make sure that we
		// only are looking at things that are possible.
		Lxy_weight_calculator2D lxy_weight(reader);

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

		// How often, for the generated sample, a pair of beta1, beta2 vpions reaches the HCal.
		// This is done at generation lifetime, so this will be the baseline which we scale against
		// in the tau loop below.
		vector<unique_ptr<TH2F>> h_gen_ratio;
		if (config._beta_type == BetaShapeType::FromMC) {
			auto r = GetFullBetaShape(config._tau_gen, n_tau_loops_at_gen, reader, lxy_weight);
			h_gen_ratio = DivideShape(r, "h_Ngen_ratio", "Fraction of events in beta space at raw generated ctau");
		}

		// And how many events actually are in the signal regions at generation?
		auto passedEventsAtGen = CalcPassedEvents(reader, vector<unique_ptr<TH2F>>(), false);

		// Count the total number of events, taking into account all weighting (like pileup, etc.).
		auto generatedEventsWithWeightsInRegions = CalcPassedEvents(reader, vector<unique_ptr<TH2F>>(), true);
		auto totalGeneratedEvents = generatedEventsWithWeightsInRegions[0];
		cout << " Total Generated Events: " << totalGeneratedEvents << endl;
		for (int i = 0; i < 4; i++) {
			cout << " Total Events in Region " << i << ": " << passedEventsAtGen[i] << endl;
		}
		
		// Next, do a double check to make sure our preselection isn't eliminating any of our signal.
		auto crossCheckNumberOfEventsInRegions = GenericCalcPassedEvents(reader, true);
		for (int i = 0; i < 4; i++) {
			if (crossCheckNumberOfEventsInRegions[i] != passedEventsAtGen[i]) {
				cout << " ** ERROR - in region " << i << " the number of events passed " << crossCheckNumberOfEventsInRegions[i] << " does not match number after preselection " << passedEventsAtGen[i] << endl;
			}
		}

		// Loop over proper lifetime
		vector<vector<unique_ptr<TH2F> > > ctau_cache; // Cache of ctau beta plots to be written out later.
		for (unsigned int i_tau = 0; i_tau < tau_binning.nbin(); i_tau++) {
			auto tau = h_res_eff[0]->GetBinCenter(i_tau+1); // Recal ROOT indicies bins at 1

			vector<doubleError> passedEventsAtTau;
			if (config._beta_type == BetaShapeType::FromMC) {
				// Get the full Beta shape
				auto rtau = GetFullBetaShape(tau, tau_loops(tau), reader, lxy_weight);
				ostringstream ctau_ratio_name;
				ctau_ratio_name << "h_ctau_ratio_" << tau << "_";
				auto h_caut_ratio = DivideShape(rtau, ctau_ratio_name.str(), ctau_ratio_name.str());

				// Now, create a weighting histogram. This is just the differece between the numerators at the
				// extrapolated ctau and at the generated ctau
				decltype(h_caut_ratio) h_Nratio;
				for (int i_region = 0; i_region < 4; i_region++) {
					auto h = unique_ptr<TH2F>(static_cast<TH2F*>(h_caut_ratio[i_region]->Clone()));
					h->Divide(h_gen_ratio[i_region].get());
					h_Nratio.push_back(move(h));
				}

				if (i_tau % 1 == 0) {
					ctau_cache.push_back(move(h_caut_ratio));
				}

				// The the number of events that passed for this lifetime.
				passedEventsAtTau = CalcPassedEvents(reader, h_Nratio, false);
			}
			else {
				// Just do Lxy scaling
				passedEventsAtTau = CalcPassedEventsLxy(reader, tau, lxy_weight);
			}

			// Calculate proper asymmetric errors and save the extrapolation result for the change in efficency.
			for (int i_region = 0; i_region < 4; i_region++) {
				doubleError eff = passedEventsAtTau[i_region] / totalGeneratedEvents;
				std::pair<Double_t, Double_t> bayerr_sig_perc = getBayes(passedEventsAtTau[i_region], totalGeneratedEvents);
				Double_t erro = (bayerr_sig_perc.first + bayerr_sig_perc.second)*0.5;

				h_res_eff[i_region]->SetBinContent(i_tau + 1, eff.value());
				h_res_eff[i_region]->SetBinError(i_tau + 1, erro);
				SetAsymError(g_res_eff[i_region], i_tau, tau, eff.value(), bayerr_sig_perc);
			}

			cout << " tau = " << tau << " npassed = " << passedEventsAtTau[0] << " passed tau/gen = " << passedEventsAtTau[0] / passedEventsAtGen[0] << " global eff = " << passedEventsAtTau[0] / totalGeneratedEvents << endl;
		}

		// Save plots in the output file
		auto output_file = unique_ptr<TFile>(TFile::Open(config._output_filename.c_str(), "RECREATE"));
		for (int i_region = 0; i_region < 4; i_region++) {
			output_file->Add(h_res_eff[i_region]);
		}

		// Save the Lxy efficiency plot
		for (int i = 0; i < 4; i++) {
			output_file->Add(lxy_weight.clone_weight(i).release());
		}

		// The default as-generated beta shape, along with a few check-points.
		if (config._beta_type == BetaShapeType::FromMC) {
			for (int i = 0; i < 4; i++) {
				auto h = static_cast<TH2F*>(h_gen_ratio[i]->Clone());
				h->SetDirectory(nullptr);
				output_file->Add(h);
			}
			for (const auto &ctaus : ctau_cache) {
				for (int i = 0; i < 4; i++) {
					auto h = static_cast<TH2F*>(ctaus[i]->Clone());
					h->SetDirectory(nullptr);
					output_file->Add(h);
				}
			}
		}

		// Save basic information for the generated sample.
		output_file->Add(save_as_histo("generated_ctau", config._tau_gen).release());
		output_file->Add(save_as_histo("n_passed_as_generated", passedEventsAtGen).release());
		output_file->Add(save_as_histo("n_as_generated", generatedEventsWithWeightsInRegions).release());

		vector<doubleError> effAtGen;
		for (int i = 0; i < 4; i++) {
			effAtGen.push_back(passedEventsAtGen[i] / totalGeneratedEvents);
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
	// Setup the command line arguments
	Args args({
		// Files that drive this proces
		Arg("muonTreeFile", "m", "The input MC file containing LLP and pass/fail info", Arg::Ordinality::Required),
		Arg("output", "f", "Output file to write results out", Arg::Ordinality::Required),

		// CTau that is required
		Arg("ctau", "c", "The ctau that the MC file was generated at", Arg::Ordinality::Required),

		// Options
		Flag("UseFlatBeta", "b", "Use only Lxy to do the extrapolation", Arg::Is::Optional),
	});

	// Make sure we got all the command line arguments we need
	if (argc == 1 || !args.Parse(argc, argv)) {
		cout << args.Usage("PlotSingleLimit") << endl;
		throw runtime_error("Bad command line arguments - exiting");
	}

	// Extract the info.

	extrapolate_config r;
	r._muon_tree_root_file = args.Get("muonTreeFile");
	r._output_filename = args.Get("output");
	r._tau_gen = args.GetAsFloat("ctau");
	r._beta_type = args.IsSet("UseFlatBeta") ? BetaShapeType::Unity : BetaShapeType::FromMC;

	return r;
}

// Initalize and populate the tau decay table.
// Due to the fact we run out of stats, this is, by its very nature, not equal binning.
variable_binning_builder PopulateTauTable()
{
	variable_binning_builder r(0.0);
#ifdef TEST_RUN
	r.bin_up_to(0.6, 0.6);
	r.bin_up_to(0.8, 0.2);
	r.bin_up_to(1.68, 0.88);
#else
	r.bin_up_to(0.6, 0.6);
	r.bin_up_to(0.8, 0.2);
	r.bin_up_to(1.68, 0.88);
	//r.bin_up_to(0.6, 0.005);
	//r.bin_up_to(4.0, 0.05);
	//r.bin_up_to(10.0, 0.2);
	//r.bin_up_to(50.0, 1.0);
#endif
	return r;
}

// Binning we will use for beta histograms
variable_binning_builder PopulateBetaBinning()
{
	variable_binning_builder beta_binning(0.0);
	beta_binning.bin_up_to(0.8, 0.2);
	beta_binning.bin_up_to(0.95, 0.05);
	beta_binning.bin_up_to(1.0, 0.005);
	return beta_binning;

	// for 50 GeV scalar, at 300 GeV, beta = 0.986
	// for 50 GeV scalar, at 200 GeV, beta = 0.968
	// for 50 GeV scalar, at 100 GeV, beta = 0.866
}
variable_binning_builder PopulateBetaBinningForUnity()
{
	variable_binning_builder beta_binning(0.0);
	beta_binning.bin_up_to(1.0, 1.0);
	return beta_binning;
}

// Window out events that will never contribute to the lifetime no matter what ctau they are
// re-simulated at. This main point for this is too prevent us from rolling the dice on the tau's
// and it will significantly speed this up.
bool doMCPreselection(const muon_tree_processor::eventInfo &entry)
{
	// These cuts derived by looking at eff for signalA region.
	// https://1drv.ms/u/s!AnlM9ZYrD4WgtXmIm59h9tbdK9WG?wd=target%282015%20Analysis%2FAnalysis%20Topics.one%7C2291C3ED-E8E5-49AA-9C56-882D0513A351%2FClosure%20Test%7CBCD03136-03E7-442E-8CC5-92CFCBA29830%2F%29

	return abs(entry.vpi1_eta) <= 3.0
		&& abs(entry.vpi2_eta) <= 3.0
		&& entry.vpi1_pt / 1000.0 > 0.0
		&& entry.vpi2_pt / 1000.0 > 0.0;
}

// Sample from the proper lifetime tau for a specific lifetime, and then do the special relativity
// calculation to understand where it ended up.
// tau - is in units of meters.
bool doSR(const caching_tlz &vpi1, const caching_tlz &vpi2, Double_t tau, Double_t &beta1, Double_t &beta2, Double_t &L2D1, Double_t &L2D2) {

	beta1 = vpi1.Beta();
	beta2 = vpi2.Beta();
	Double_t gamma1 = vpi1.Gamma();
	Double_t gamma2 = vpi2.Gamma();

	// Get ctau of the two we are to simulate, in meters.
	Double_t ct1 = gRandom->Exp(tau);
	Double_t ct2 = gRandom->Exp(tau);

	// What is the decay length in the lab frame (in meters)?
	Double_t ct1prime = gamma1 * ct1;
	Double_t ct2prime = gamma2 * ct2;
	Double_t lxy1 = beta1 * ct1prime;
	Double_t lxy2 = beta2 * ct2prime;

	// Calculate the timing in nano-seconds.
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

	// Timing restrictions. This first test should never fire
	// because there is no way for the particle to go faster than "c", and
	// that is the only way to have negative timing. But we keep it here for completness
	// sake.
	if (deltat1 < -3.0 || deltat2 < -3.0) {
		return false;
	}
	if (deltat1 > 15.0 || deltat2 > 15.0) {
		return false;
	}

	return true;
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

	auto beta_binning = PopulateBetaBinning();

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
		TLorentzVector vpi1_tlz, vpi2_tlz;
		vpi1_tlz.SetPtEtaPhiE(entry.vpi1_pt / 1000.0, entry.vpi1_eta, entry.vpi1_phi, entry.vpi1_E / 1000.0);
		vpi2_tlz.SetPtEtaPhiE(entry.vpi2_pt / 1000.0, entry.vpi2_eta, entry.vpi2_phi, entry.vpi2_E / 1000.0);

		auto vpi1 = caching_tlz(vpi1_tlz);
		auto vpi2 = caching_tlz(vpi2_tlz);

		for (Int_t maketaus = 0; maketaus < ntauloops; maketaus++) { // tau loop to generate toy events

			Double_t beta1 = -1, beta2 = -1, L2D1 = -1, L2D2 = -1;

			// Do SR, apply SR related cuts (like timing).
			if (doSR(vpi1, vpi2, tau, beta1, beta2, L2D1, L2D2)) {
				den->Fill(beta1, beta2, entry.weight);
				for (int i_region = 0; i_region < 4; i_region++) {
					num[i_region]->Fill(beta1, beta2, entry.weight * lxyWeight(i_region, L2D1, L2D2));
				}
			}
		}
	});
	return make_pair(move(num), move(den));
}

vector<doubleError> CalcPassedEventsLxy(const muon_tree_processor &mc_entries, double tau, Lxy_weight_calculator &lxyWeight)
{
	// The resulting sums are just all the weights added together.
	vector<doubleError> results(4);

	size_t nloops = 100;

	// Loop over each MC entry, and generate tau's at several different places
	int count = 0;
	mc_entries.process_all_entries([&count, &results, nloops, tau, &lxyWeight](const muon_tree_processor::eventInfo &entry) {
#ifdef notyet
		for (int i_region = 0; i_region < 4; i_region++) {
			results[i_region] += lxyWeight(i_region, entry.vpi1_Lxy/1000.0, entry.vpi2_Lxy/1000.0);
		}
#else
		TLorentzVector vpi1_tlz, vpi2_tlz;
		vpi1_tlz.SetPtEtaPhiE(entry.vpi1_pt / 1000.0, entry.vpi1_eta, entry.vpi1_phi, entry.vpi1_E / 1000.0);
		vpi2_tlz.SetPtEtaPhiE(entry.vpi2_pt / 1000.0, entry.vpi2_eta, entry.vpi2_phi, entry.vpi2_E / 1000.0);

		auto vpi1 = caching_tlz(vpi1_tlz);
		auto vpi2 = caching_tlz(vpi2_tlz);

		for (Int_t maketaus = 0; maketaus < nloops; maketaus++) { // tau loop to generate toy events

			Double_t beta1 = -1, beta2 = -1, L2D1 = -1, L2D2 = -1;

			// Do special relativity, apply cuts as needed.
			if (doSR(vpi1, vpi2, tau, beta1, beta2, L2D1, L2D2)) {
				for (int i_region = 0; i_region < 4; i_region++) {
					results[i_region] += entry.weight * lxyWeight(i_region, L2D1, L2D2);
				}
			}
		}
#endif
	});

#ifndef notyet
	for (int i_region = 0; i_region < 4; i_region++) {
		results[i_region] = results[i_region] / nloops;
	}
#endif
	return results;
}

// Calculate a 2D efficiency given a denominator and the numerator selected from the denominator
template<class T>
vector<unique_ptr<T>> DivideShape(const pair<vector<unique_ptr<T>>, unique_ptr<T>> &r, const string &name, const string &title)
{
	vector<unique_ptr<T>> result;
	char region = 'A';
	for (auto &info : r.first) {
		unique_ptr<T> ratio(static_cast<T*>(info->Clone()));
		ratio->Divide(info.get(), r.second.get(), 1.0, 1.0, "B");
		ratio->SetNameTitle((name + region).c_str(), (title + " " + region + "; Beta; Beta").c_str());
		result.push_back(move(ratio));
		region++;
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
	}, !eventCountOnly);

	return nEvents;
}

// Do a very generic calculation on the number of events that have passed.
vector<doubleError> GenericCalcPassedEvents(const muon_tree_processor &reader, bool ignore_preselection)
{
	vector<doubleError> nEvents(4);

	reader.process_all_entries([&nEvents](const muon_tree_processor::eventInfo &entry) {

		doubleError weight(entry.weight, entry.weight);

		int i_region = entry.RegionA ? 0
			: entry.RegionB ? 1
			: entry.RegionC ? 2
			: entry.RegionD ? 3
			: -1;
		if (i_region >= 0) {
			nEvents[i_region] += weight;
		}
	}, !ignore_preselection);

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
