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

////////////////////////////////
// Main entry point.
////////////////////////////////
int main(int argc, char**argv)
{
	int dummy_argc = 0;
	auto a = unique_ptr<TApplication>(new TApplication("extrapolate_betaw", &dummy_argc, argv));
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

			g_res_eff.push_back(unique_ptr<TGraphAsymmErrors>(new TGraphAsymmErrors(tau_binning.nbin())));
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
		auto output_file = unique_ptr<TFile>(TFile::Open("extrapolate_betaw_results.root", "RECREATE"));
		for (int i_region = 0; i_region < 4; i_region++) {
			output_file->Add(h_res_eff[i_region]);
		}
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
	if (argc != 3) {
		cout << "Usage: extrapolate_betaw <muonTree-file> <generated-ctau>" << endl;
		cout << endl;
		cout << "    <muonTree-file>           TFile containing the muonTree root file generated from a particular MC sample." << endl;
		cout << "    <generated-ctau>          The lifetime (in meters) where the sample was generated." << endl;

		throw runtime_error("Wrong number of arguments");
	}

	extrapolate_config r;
	r._muon_tree_root_file = argv[1];
	r._tau_gen = atof(argv[2]);
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

class bogus
{
public:
	bogus() {}

};

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
	unique_ptr<TH2F> den(new TH2F(dname.str().c_str(), dname.str().c_str(), beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list()));
	vector<unique_ptr<bogus>> dork;
	vector<unique_ptr<TH2F>> num;
	for (int i_region = 0; i_region < 4; i_region++) {
		nname << "A";
		num.push_back(unique_ptr<TH2F>(new TH2F(nname.str().c_str(), nname.str().c_str(), beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list())));
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

	auto h_num = unique_ptr<TH1D>(new TH1D("h_num", "", 1, 0, 1));
	auto h_den = unique_ptr<TH1D>(new TH1D("h_den", "", 1, 0, 1));

	h_num->SetBinContent(1, num.value());
	h_den->SetBinContent(1, den.value());
	h_num->SetBinError(1, num.err());
	h_den->SetBinError(1, den.err());

	auto h_eff = unique_ptr<TGraphAsymmErrors>(new TGraphAsymmErrors());
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

#ifdef notyet
TH1F *h_res_eff, *h_res_ev;

TH2F *h_Npass_gen = 0;
TH2F *h_Ngen_num = 0;
TH2F *h_Ngen_den = 0;

TGraphAsymmErrors *g_res_eff;
TGraphAsymmErrors *g_res_ev;
TGraphAsymmErrors *g_res_landau;
void setParams(std::string mp, Double_t &mass, Double_t &xsec, Double_t &effi_real, int &k_gen, Double_t &tau_gen, Double_t &scale_factor);
void doSR(TLorentzVector vpi1, TLorentzVector vpi2, Double_t tau, Double_t &beta1, Double_t &beta2, Double_t &L2D1, Double_t &L2D2);

TH2F *ef12;

ofstream logfile;


// Hold onto variable binning for beta.
variable_binning_builder beta_binning(0.0);

std::pair<Double_t, Double_t> getBayes(Double_t num, Double_t den) {
	// returns the Bayesian uncertainty over the num/den ratio
	std::pair<Double_t, Double_t> result(0., 0.);

	TH1D *h_num = new TH1D("h_num", "", 1, 0, 1);
	TH1D *h_den = new TH1D("h_den", "", 1, 0, 1);

	h_num->SetBinContent(1, num);
	h_den->SetBinContent(1, den);

	TGraphAsymmErrors *h_eff = new TGraphAsymmErrors();
	h_eff->BayesDivide(h_num, h_den);//, "w"); 

	result.first = h_eff->GetErrorYlow(0);
	result.second = h_eff->GetErrorYhigh(0);

	delete h_num;
	delete h_den;
	delete h_eff;

	return result;
}

void extrapolate_betaw(string mp, string file = "testing", Double_t timecutshift = 0.0, Double_t ETcutshift = 0.0, int ncount = -1.0, bool triggerOnly = false, bool use_lxyeff = false, string inputFileSuffix = "_Rachel_3D")
{
	lxyUseEffWeight = use_lxyeff;
	int ntauloops = 5;
	bool onlyGenTau = false;
	TString outputroot = mp + "/" + mp + "_" + file + ".root";
	std::string output = mp + "/" + mp + "_" + file + ".txt";
	std::string instr = "input/" + mp + inputFileSuffix + ".root";
	std::string sfhist;

	logfile.open(output.c_str(), ofstream::out | ofstream::trunc);

	logfile << "***Processing: "
		<< "  -> Mass = " << mp << endl
		<< "  -> Input file " << instr << endl
		<< "  -> Time cut shift = " << timecutshift << endl
		<< "  -> ET cut shift = " << ETcutshift << endl
		<< "  -> Tau Loops: " << ntauloops << endl
		<< "  -> # Events in file to use (-1 is all): " << ncount << endl 
		<< "  -> Using lxy efficiency: " << lxyUseEffWeight << endl
		<< endl;

	// We want uneven binning for beta.

	beta_binning.add_bins(4, 0.2);
	beta_binning.add_bins(4, 0.05);
	h_Npass_gen = new TH2F("h_Npass_gen", "h_Npass_gen", beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list());

	// Proper lifetimes used in this simulation (in meters)
	Int_t nbin(0);
	Double_t *xbin = PopulateTauTable(nbin);

	// Load the file and attach to testtree.
	LoadAndAttachFile(instr, ncount);

	// Get everything configured from our inputs.

	Double_t lumi = 20300.; //Integrated luminosity in pb-1 for N events VS lifetime  
	Double_t mass = -1, xsec = -1, effi_real = -1, tau_gen = -1, scale_factor = -1;
	int k_gen = -1;
	setParams(mp, mass, xsec, effi_real, k_gen, tau_gen, scale_factor);  // Set parameters for mass point mp
	Double_t Lxy_min = 0, Lxy_max = 10000 / 1000.;

	// Dump info about the sample so we have a log of what we are doing.

	logfile << "***Generated Sample: "
		<< "  -> Tau = " << tau_gen << endl
		<< "  -> Efficiency =  " << effi_real << endl
		<< "  -> Cross-section = " << xsec << endl
		<< "  -> # of events = " << nentries << endl << endl;

	// ************************************************************************************** //
	// First, we want the number of events passing all cuts as a function of beta1 and beta2.
	// ************************************************************************************** //

	// Get the number of events that have passed, and keep track with a histo...
	doubleError passedEventsAtGen = CalcPassedEvents(nullptr, false, h_Npass_gen, triggerOnly);

	// Get the total number of events
	doubleError totalEventsAtGen = CalcPassedEvents(nullptr, true, nullptr, triggerOnly);

	logfile << "Finished initial scan of MC: " << endl
		<< "  -> Total number of events passing: " << passedEventsAtGen << endl
		<< "  -> Total number of generated events: " << totalEventsAtGen << endl
		<< "  -> Global eff fraction at generated ctau: " << passedEventsAtGen / totalEventsAtGen << endl;

	// **************************************************************************************************** //
	// Second, we want how often, for the generated sample, a pair of beta1, beta2 vpions reaches the HCal.
	// **************************************************************************************************** //

	pair<TH2F*, TH2F*> r = GetFullBetaShape(tau_gen, ntauloops, "h_Ngen");
	TH2F *h_gen_ratio = DivideShape(r, "h_Ngen_ratio", "Fraction of events in beta space at raw generated ctau");
	h_Ngen_num = r.first;
	h_Ngen_den = r.second;

	// **************************************************************************************** //
	// Third, we want how often, for some ctau, a pair of beta1, beta2 vpions reaches the HCal.
	// **************************************************************************************** //

	vector<TH2F*> toSave;
	h_res_eff = new TH1F("h_res_eff", "h_res_eff", nbin, xbin); //Efficiency VS lifetime
	h_res_ev = new TH1F("h_res_ev", "h_res_ev", nbin, xbin);   //Number of events VS lifetime

	g_res_eff = new TGraphAsymmErrors(nbin);
	g_res_eff->SetName("g_res_eff");
	g_res_eff->SetTitle("Absolute efficiency to the generated analysis");

	g_res_ev = new TGraphAsymmErrors(nbin); //Number of events VS lifetime
	g_res_ev->SetName("g_res_ev");
	g_res_ev->SetTitle("Number of events as a function of lifetime");

	TH1F *h_relative_eff = new TH1F("h_relative_eff", "Relative Efficience compared to generated ctau", nbin, xbin);

	for (int k = 0; k < h_res_eff->GetNbinsX(); k++) {

		Double_t tau = h_res_eff->GetBinCenter(k + 1);

		TH2F *h_Nctau_num = 0;
		TH2F *h_Nctau_den = 0;
		TH2F *h_Nratio = 0;

		// Get the beta shape for these guys that pass at this proper lifetime.

		ostringstream name;
		name << "h_Ncatu_" << tau;
		pair<TH2F*, TH2F*> r = GetFullBetaShape(tau, ntauloops, name.str());
		h_Nctau_num = r.first;
		h_Nctau_den = r.second;

		ostringstream rname;
		rname << "h_ratio_ctau_" << tau;

		TH2F *h_caut_ratio = DivideShape(r, rname.str(), rname.str());

		// Now, create a weighting histogram. This is just the differece between the numerators at the
		// extrapolated ctau and at the generated ctau

		h_Nratio = (TH2F*)h_caut_ratio->Clone();
		h_Nratio->Divide(h_gen_ratio);
		ostringstream rrname;
		rrname << "h_Nratio_" << tau;
		h_Nratio->SetName(rrname.str().c_str());
		h_Nratio->SetTitle(rrname.str().c_str());

		// Next, re-calc the efficiency

		doubleError passedEventsAtTau = CalcPassedEvents(h_Nratio, false, nullptr, triggerOnly);
		doubleError relativeEff = passedEventsAtTau / passedEventsAtGen;
		doubleError eff = passedEventsAtTau / totalEventsAtGen;

		std::pair<Double_t, Double_t> bayerr_sig_perc = getBayes(passedEventsAtTau, totalEventsAtGen);
		Double_t erro = (bayerr_sig_perc.first + bayerr_sig_perc.second)*0.5;

		h_res_eff->SetBinContent(k + 1, eff.value());
		h_res_eff->SetBinError(k + 1, erro);
		SetAsymError(g_res_eff, k, tau, eff.value(), bayerr_sig_perc);

		h_res_ev->SetBinContent(k + 1, (eff * xsec * lumi).value());
		h_res_ev->SetBinError(k + 1, erro * xsec * lumi);
		pair<double, double> bayerr_nev = make_pair(bayerr_sig_perc.first * xsec * lumi, bayerr_sig_perc.second * xsec * lumi);
		SetAsymError(g_res_ev, k, tau, (eff * xsec * lumi).value(), bayerr_nev);

		h_relative_eff->SetBinContent(k + 1, relativeEff.value());
		h_relative_eff->SetBinError(k + 1, relativeEff.err());

		logfile << "k = " << k << " tau = " << tau << " npassed = " << passedEventsAtTau << " passed tau/gen = " << passedEventsAtTau / passedEventsAtGen << " global eff = " << eff << endl;

		if (k % 50 == 0) {
			toSave.push_back(h_Nctau_den);
			toSave.push_back(h_Nctau_num);
			toSave.push_back(h_Nratio);
			toSave.push_back(h_caut_ratio);
		}
		else {
			delete h_Nctau_num;
			delete h_Nctau_den;
			delete h_Nratio;
		}
	} // k loop

	g_res_landau = new TGraphAsymmErrors(nbin); //Number of events VS lifetime: landau fit
	g_res_landau->SetName("g_res_landau");

	TFile *outf = TFile::Open(outputroot, "RECREATE");

	g_res_eff->Write();

	h_res_eff->Write();
	h_res_ev->Write();
	g_res_ev->Write();

	h_Npass_gen->Write();
	h_Ngen_num->Write();
	h_Ngen_den->Write();
	h_gen_ratio->Write();
	h_relative_eff->Write();

	for (size_t i = 0; i < toSave.size(); i++) {
		toSave[i]->Write();
	}

	outf->Close();
	logfile.close();
}

//
// Parse the command line arguments
//
extrapolate_config parse_command_line(int argc, char ** argv)
{
	if (argc != 2) {
		throw exception("Incorrect number of arguments.");
	}
	auto r = extrapolate_config();
	r._muon_tree_root_file = argv[1];
	return r;
}

void setParams(std::string mp, Double_t &mass, Double_t &xsec, Double_t &effi_real, int &k_gen, Double_t &tau_gen, Double_t &scale_factor) {

	if (mp.compare("100_10") == 0) {
		mass = 10.0;  //Mass in GeV of the vpion (159220)  
		xsec = 29.68;  //X-section in pb for N events VS lifetime 159220 (100_10)
		effi_real = 0.001134;  //At generated lifetime 450
		k_gen = 87;
		tau_gen = 0.45;
		scale_factor = 1; //1.00013;
	}
	else if (mp.compare("100_25") == 0) {
		mass = 25.0;  //Mass in GeV of the vpion (159221)  
		xsec = 29.68;  //X-section in pb for N events VS lifetime 159221 (100_25)
		effi_real = 0.000630;  //At generated lifetime 1250
		k_gen = 132;
		tau_gen = 1.25;
		scale_factor = 1; //1.02033;
	}
	else if (mp.compare("126_10") == 0) {
		mass = 10.0;  //Mass in GeV of the vpion (159222)  
		xsec = 18.97;  //X-section in pb for N events VS lifetime 159222 (126_10)
		effi_real = 0.002655;  //At generated lifetime 350
		k_gen = 68;
		tau_gen = 0.35;
		scale_factor = 1; //1.00000;
	}
	else if (mp.compare("126_25") == 0) {
		mass = 25.0;  //Mass in GeV of the vpion (159223)  
		xsec = 18.97;  //X-section in pb for N events VS lifetime 159223 (126_25)
		effi_real = 0.002145;  //At generated lifetime 900
		k_gen = 124;
		tau_gen = 0.9;
		scale_factor = 1; //1.00298;
	}
	else if (mp.compare("126_40") == 0) {
		mass = 40.0;  //Mass in GeV of the vpion (159224)  
		xsec = 18.97;  //X-section in pb for N events VS lifetime 159224 (126_40)
		effi_real = 0.001116;  //At generated lifetime 1850
		k_gen = 148;
		tau_gen = 1.85;
		scale_factor = 1; //1.11395;  //!! LARGE!
	}
	else if (mp.compare("140_10") == 0) {
		mass = 10.0;  //Mass in GeV of the vpion (159225)  
		xsec = 15.42;  //X-section in pb for N events VS lifetime 159225 (140_10)
		effi_real = 0.003243;  //At generated lifetime 275
		k_gen = 54;
		tau_gen = 0.275;
		scale_factor = 1; //1.00000;
	}
	else if (mp.compare("140_20") == 0) {
		mass = 20.0;  //Mass in GeV of the vpion (158346)  
		xsec = 15.42;  //X-section in pb for N events VS lifetime 158346 (140_20)
		effi_real = 0.003077;  //At generated lifetime 630
		k_gen = 117;
		tau_gen = 0.63;
		scale_factor = 1; //1.00016;
	}
	else if (mp.compare("140_40") == 0) {
		mass = 40.0;  //Mass in GeV of the vpion (159226)  
		xsec = 15.42;  //X-section in pb for N events VS lifetime 159226 (140_40)
		effi_real = 0.002060;  //At generated lifetime 1500
		k_gen = 139;
		tau_gen = 1.5;
		scale_factor = 1; //1.04871;
	}
	else if (mp.compare("300_50") == 0) {
		mass = 50.0;  //Mass in GeV of the vpion (300_50)  
		xsec = 3.594;  //X-section in pb for N events VS lifetime (300_50)
		effi_real = 0.0060467;  //At generated lifetime 800
		k_gen = 121;
		tau_gen = 0.8;
		scale_factor = 1; //1.00095;
	}
	else if (mp.compare("600_50") == 0) {
		mass = 50.0;  //Mass in GeV of the vpion (600_50)  
		xsec = 0.523;  //X-section in pb for N events VS lifetime 159226 (600_50)
		effi_real = 0.005001;  //At generated lifetime 500
		k_gen = 96;
		tau_gen = 0.5;
		scale_factor = 1; //1.00000;
	}
	else if (mp.compare("600_150") == 0) {
		mass = 150.0;  //Mass in GeV of the vpion (600_150)  
		xsec = 0.523;  //X-section in pb for N events VS lifetime (600_150)
		effi_real = 0.004066;  //At generated lifetime 1700
		k_gen = 144;
		tau_gen = 1.7;
		scale_factor = 1; //1.07390;
	}
	else if (mp.compare("900_50") == 0) {
		mass = 50.0;  //Mass in GeV of the vpion (600_50)  
		xsec = 0.0571;  //X-section in pb for N events VS lifetime 159226 (600_50)
		effi_real = 0.004746;  //At generated lifetime 500
		k_gen = 96;
		tau_gen = 0.5;
		scale_factor = 1; //1.00000;
	}
	else if (mp.compare("900_150") == 0) {
		mass = 150.0;  //Mass in GeV of the vpion (900_150)  
		xsec = 0.0571;  //X-section in pb for N events VS lifetime (900_150)
		effi_real = 0.003711;   //At generated lifetime 1300
		k_gen = 134;
		tau_gen = 1.3;
		scale_factor = 1; //1.03640;
	}
	else {
		std::cout << mp << " is not a valid mass point." << std::endl;
	}
	return;
}
#endif
