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

#pragma warning (push)
#pragma warning (disable: 4244)
#include "TLorentzVector.h"
#pragma warning (pop)

#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include <fstream>
#include "TRandom2.h"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include <memory>

using namespace std;

// Some config constants

// How many loops in tau should we do?
size_t n_tau_loops = 5;

// Helper methods
struct extrapolate_config {
	string _muon_tree_root_file;
};
extrapolate_config parse_command_line(int argc, char **argv);
variable_binning_builder PopulateTauTable();

////////////////////////////////
// Main entry point.
////////////////////////////////
int main(int argc, char**argv)
{
	try {
		// Pull out the various command line arguments
		auto config = parse_command_line(argc, argv);

		// Create the muon tree reader object
		auto reader = new muon_tree_processor(config._muon_tree_root_file);

		// Create the histograms we will use to store the raw results.
		auto tau_binning = PopulateTauTable();
		auto h_res_eff = new TH1F("h_res_eff", "h_res_eff", tau_binning.nbin(), tau_binning.bin_list()); //Efficiency VS lifetime
		auto h_res_ev = new TH1F("h_res_ev", "h_res_ev", tau_binning.nbin(), tau_binning.bin_list());   //Number of events VS lifetime

		// Loop over proper lifetime
		for (unsigned int i_tau; i_tau < tau_binning.first; i_tau++) {
			auto tau = h_res_eff->GetBinCenter(i_tau+1); // Recal ROOT indicies bins at 1

			// Get the full Beta shape
			auto r = GetFullBetaShape(tau, n_tau_loops);
		}
	}
	catch (exception &e)
	{
		cout << "Failed!!" << endl;
		cout << " --> " << e.what() << endl;
		return 1;
	}
	return 0;
}

// Initalize and populate the tau decay table.
// Due to the fact we run out of stats, this is, by its very nature, not equal binning.
variable_binning_builder PopulateTauTable()
{
	variable_binning_builder r(0.0);
	r.bin_up_to(0.6, 0.005);
	r.bin_up_to(4.0, 0.05);
	r.bin_up_to(50.0, 0.2);
	return r;
}

// Generate a lxy1 and lxy2 set of histograms. The denominator (first item) is
// the generated ones. The numerator is modified by the weight histogram.
pair<unique_ptr<TH2F>, unique_ptr<TH2F>> GetFullBetaShape(double tau, int ntauloops)
{
	// Create numerator and denominator histograms.
	// To avoid annoying ROOT error messages, make a unique name for each.
	ostringstream dname, nname;
	dname << "tau_" << tau << "_den";
	nname << "tau_" << tau << "_num";
	TH2F *den = new TH2F(dname.str().c_str(), dname.str().c_str(), beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list());
	TH2F *num = new TH2F(dname.str().c_str(), dname.str().c_str(), beta_binning.nbin(), beta_binning.bin_list(), beta_binning.nbin(), beta_binning.bin_list());
	den->Sumw2();
	num->Sumw2();

	for (Int_t n = 0; n < nentries; n++) {

		testtree->GetEntry(n);
		if (hasGluons) continue;

		TLorentzVector vpi1, vpi2;
		vpi1.SetPtEtaPhiE(vpi1_pt, vpi1_eta, vpi1_phi, vpi1_E);
		vpi2.SetPtEtaPhiE(vpi2_pt, vpi2_eta, vpi2_phi, vpi2_E);

		for (Int_t maketaus = 0; maketaus < ntauloops; maketaus++) { // tau loop to generate toy events

			Double_t beta1 = -1, beta2 = -1, L2D1 = -1, L2D2 = -1;

			doSR(vpi1, vpi2, tau, beta1, beta2, L2D1, L2D2); // Do special relativity

			den->Fill(beta1, beta2, puweight);
			double wt = lxyWeight(L2D1, L2D2);
			num->Fill(beta1, beta2, puweight * wt);
		}  // maketaus to get h_Ngen_ctau
	}  // nentries to get the same

	return make_pair(num, den);
}

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


// Load the file and attach the branches. We do this as global variables (yeah, yuck), but that means
// we can also separate out code and make it easier to re-run.
Long64_t nentries = 0;
bool hasGluons;
int PassedCalRatio;
bool hasSecondJet_case1, hasSecondJet_case2;
float vpi1_pt, vpi1_eta, vpi1_phi, vpi1_E;
float vpi2_pt, vpi2_eta, vpi2_phi, vpi2_E;
float puweight;
TTree *testtree = 0;
void LoadAndAttachFile(const string &instr, const Long64_t maxEventsToRun) {
	TFile *input = TFile::Open(instr.c_str());
	if (!input->IsOpen()) {
		cout << "ERROR: Unable to open file!" << endl;
		abort();
	}

	// The 2D lxy efficiency histogram (normalize it)
	ef12 = (TH2F*)input->Get("Final_events/effi_Lxy1_Lxy2");
	if (ef12 == nullptr){
		cout << "Unable to load the effi_Lxy1_Lxy2 histogram!" << endl;
		abort();
	}
	double totprop = ef12->Integral();
	ef12->Scale(1.0 / totprop);

	testtree = (TTree*)input->Get("muonTree");
	nentries = testtree->GetEntries();
	if (maxEventsToRun > 0) {
		nentries = min(maxEventsToRun, nentries);
	}

	TBranch *b_hasGluons = testtree->GetBranch("hasGluons");
	TBranch *b_PassedCalRatio = testtree->GetBranch("PassedCalRatio");
	TBranch *b_hasSecondJet_case1 = testtree->GetBranch("hasSecondJet_case1");
	TBranch *b_hasSecondJet_case2 = testtree->GetBranch("hasSecondJet_case2");
	TBranch *b_vpi1_pt = testtree->GetBranch("vpi1_pt");
	TBranch *b_vpi1_eta = testtree->GetBranch("vpi1_eta");
	TBranch *b_vpi1_phi = testtree->GetBranch("vpi1_phi");
	TBranch *b_vpi1_E = testtree->GetBranch("vpi1_E");
	TBranch *b_vpi2_pt = testtree->GetBranch("vpi2_pt");
	TBranch *b_vpi2_eta = testtree->GetBranch("vpi2_eta");
	TBranch *b_vpi2_phi = testtree->GetBranch("vpi2_phi");
	TBranch *b_vpi2_E = testtree->GetBranch("vpi2_E");
	TBranch *b_puweight = testtree->GetBranch("weight");
	b_hasGluons->SetAddress(&hasGluons);
	b_PassedCalRatio->SetAddress(&PassedCalRatio);
	b_hasSecondJet_case1->SetAddress(&hasSecondJet_case1);
	b_hasSecondJet_case2->SetAddress(&hasSecondJet_case2);
	b_vpi1_pt->SetAddress(&vpi1_pt);
	b_vpi1_eta->SetAddress(&vpi1_eta);
	b_vpi1_phi->SetAddress(&vpi1_phi);
	b_vpi1_E->SetAddress(&vpi1_E);
	b_vpi2_pt->SetAddress(&vpi2_pt);
	b_vpi2_eta->SetAddress(&vpi2_eta);
	b_vpi2_phi->SetAddress(&vpi2_phi);
	b_vpi2_E->SetAddress(&vpi2_E);
	b_puweight->SetAddress(&puweight);
}

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

std::pair<Double_t, Double_t> getBayes(const doubleError &num, const doubleError &den) {
	// returns the Bayesian uncertainty over the num/den ratio
	std::pair<Double_t, Double_t> result(0., 0.);

	TH1D *h_num = new TH1D("h_num", "", 1, 0, 1);
	TH1D *h_den = new TH1D("h_den", "", 1, 0, 1);

	h_num->SetBinContent(1, num.value());
	h_den->SetBinContent(1, den.value());
	h_num->SetBinError(1, num.err());
	h_den->SetBinError(1, den.err());

	TGraphAsymmErrors *h_eff = new TGraphAsymmErrors();
	h_eff->BayesDivide(h_num, h_den);//, "w"); 

	result.first = h_eff->GetErrorYlow(0);
	result.second = h_eff->GetErrorYhigh(0);

	delete h_num;
	delete h_den;
	delete h_eff;

	return result;
}

// Controls how we calc in/out.
bool lxyUseEffWeight = false;

// Do a square weight calculation.
double lxy1DWeight(double L2D)
{
	return L2D < 3.88 && L2D > 2.28 ? 1.0 : 0.0; // From an email from Daniela
	//return L2D < 4.25 && L2D > 2.25 ? 1.0 : 0.0; // Original
}

// Calc the weight for an lxy guy to be in the proper region.
double lxyWeight(double L2D1, double L2D2) {
	if (lxyUseEffWeight) {
		// look up in the weight histogram
		int xbin = ef12->FindBin(L2D1, L2D2);
		return ef12->GetBinContent(xbin);
	}
	else {
		return lxy1DWeight(L2D1) * lxy1DWeight(L2D2);
	}
}


// Calculate the beta's for the event only once. This is done because it is
// very slow to do this (the SetPtEtaPhiE is a very expensive call).
class beta_cache {
public:
	inline beta_cache()
		: _gotit(false)
	{}

	inline double beta1() { calc();  return _b1; }
	inline double beta2() { calc();  return _b2; }

private:
	bool _gotit;
	double _b1, _b2;

	void calc()
	{
		if (_gotit)
			return;
		_gotit = true;

		TLorentzVector vpi1, vpi2;
		vpi1.SetPtEtaPhiE(vpi1_pt, vpi1_eta, vpi1_phi, vpi1_E);
		vpi2.SetPtEtaPhiE(vpi2_pt, vpi2_eta, vpi2_phi, vpi2_E);

		_b1 = vpi1.Beta();
		_b2 = vpi2.Beta();
	}
};

// Calculate the number of events that pass our cuts (possibly weighted).
doubleError CalcPassedEvents(TH2F *weightHist = 0, bool eventCountOnly = false, TH2F *passHist = 0, bool triggerOnly = false)
{
	doubleError nEvents; // This will be the number of events in the TTree (without gluons)

	if (passHist) {
		passHist->Sumw2();
	}

	for (Int_t n = 0; n < nentries; n++) {

		testtree->GetEntry(n);
		if (hasGluons) continue;

		// Calculate the event weight. A combination of the pile up reweighting from the ntuple and perhaps
		// the beta re-weighting from the input histogram.
		doubleError weight(puweight, puweight);

		beta_cache b;
		if (weightHist) {
			int nbin = weightHist->FindBin(b.beta1(), b.beta2());
			weight *= doubleError(weightHist->GetBinContent(nbin), weightHist->GetBinError(nbin));
		}

		// Count the event and populate the output histogram, if
		// We are doing an event count only (e.g. the denominator) or
		// it passes our analysis cuts (e.g. the numerator).
		if (eventCountOnly
			|| (hasSecondJet_case1 || hasSecondJet_case2)
			|| (triggerOnly && PassedCalRatio)
			) {

			nEvents += weight;
			if (passHist) passHist->Fill(b.beta1(), b.beta2(), weight.value());
		}
	}

	return nEvents;
}

void SetAsymError(TGraphAsymmErrors *g, int bin, double tau, double bvalue, const pair<double, double> &assErrors)
{
	g->SetPoint(bin, tau, bvalue);
	g->SetPointEYlow(bin, assErrors.first);
	g->SetPointEYhigh(bin, assErrors.second);
}

// Given the result from the beta shape above, do a ratio,a nd rename
TH2F *DivideShape(const pair<TH2F*, TH2F*> &r, const string &name, const string &title)
{
	TH2F *ratio = (TH2F*)r.first->Clone();
	ratio->Divide(r.second);
	ratio->SetNameTitle(name.c_str(), title.c_str());

	return ratio;
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
	vpixyz1.SetMagThetaPhi(lxy1, theta1, vpi1_phi);
	vpixyz2.SetMagThetaPhi(lxy2, theta2, vpi2_phi);

	L2D1 = vpixyz1.Perp();
	L2D2 = vpixyz2.Perp();

	return;
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
