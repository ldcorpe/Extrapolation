///
/// Given input efficiency curves as a function of lifetime, and the result of a limit setting
/// (the 95% CL, 1 sigma, 2 sigma), generate a final limit plot with 1 and 3 sigma bands
///

#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "Wild/CommandLine.h"

#include <iostream>
#include <memory>
#include <string>

using namespace std;
using namespace Wild::CommandLine;

// Forward decls
unique_ptr<TGraphAsymmErrors> get_efficiency_vs_lifetime(const string &filename);

int main(int argc, char *argv[])
{
	// Setup the command line arguments
	Args args({
		// Arguments that must be present on the command line, one after the other
		Arg("extrapfile", "e", "ROOT file containing the efficiency vs ctau extrapolation curve", Arg::Ordinality::Required)
	});

	// Make sure we got all the command line arguments we need
	if (argc == 1 || !args.Parse(argc, argv)) {
		cout << args.Usage("PlotSingleLimit") << endl;
		return 1;
	}

	// Load up the input data.
	auto eff_vs_lifetime = get_efficiency_vs_lifetime(args.Get("extrapfile"));

	// Ready to go!
	cout << args.Get("extrapolationfile") << endl;
	return 1;
}

// Load up from the input file the extrapolation curve
unique_ptr<TGraphAsymmErrors> get_efficiency_vs_lifetime(const string & filename)
{
	auto file = unique_ptr<TFile>(TFile::Open(filename.c_str(), "READ"));
	auto r = unique_ptr<TFile>(static_cast<TGraphAsymmErrors*>(file->Get("")));
	return unique_ptr<TGraphAsymmErrors>();
}
