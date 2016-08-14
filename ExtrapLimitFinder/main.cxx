//
// Given an extrapolation file as input, derive the limits as a function of
// lifetime for output.
//

#include "extrap_file_wrapper.h"
#include "limitSetting.h"

#include "Wild/CommandLine.h"
#include "TApplication.h"
#include "TROOT.h"

#include <iostream>
#include <memory>
#include <string>

using namespace std;
using namespace Wild::CommandLine;

// Helper functions, etc.
struct config {
	string extrapolate_filename; // Input ROOT file to guide this whole thing

	ABCD observed_data;
	abcd_limit_config limit_settings;
};
config parse_command_line(int argc, char **argv);

// Main entry point
int main(int argc, char **argv)
{
	// Get ROOT setup so we don't have funny problems with it later.
	int dummy_argc = 0;
	auto a = unique_ptr<TApplication>(new TApplication("extrapolate_betaw", &dummy_argc, argv));
	gROOT->SetBatch(true);

	// Protect against anything weird going wrong so we get a sensible error message.
	try {
		auto c = parse_command_line(argc, argv);
		extrap_file_wrapper input_file(c.extrapolate_filename);

		// If we are going to do it by setting the limit once and calculating the efficiency.
		if (c.limit_settings.scaleLimitByEfficiency) {
			extrapolate_limit_to_lifetime_by_efficency (input_file, c.observed_data, c.limit_settings);
		}
		else {
			extrapolate_limit_to_lifetime (input_file, c.observed_data, c.limit_settings);
		}
	}
	catch (exception &e) {
		cout << "Total failure - exception thrown: " << e.what() << endl;
		return 1;
	}
	return 0;
}

// Parse command line arguments
config parse_command_line(int argc, char **argv)
{
	// Setup the command line arguments
	Args args({
		// The input MC file which contains just about all the info we could need from MC.
		Arg("extrapFile", "e", "The extrapolation root file", Arg::Is::Required),

		// We need to know how many data were in each of the A, B, C, and D regions.
		Arg("nA", "A", "How many events observed in data in region A", Arg::Is::Required),
		Arg("nB", "B", "How many events observed in data in region B", Arg::Is::Required),
		Arg("nC", "C", "How many events observed in data in region C", Arg::Is::Required),
		Arg("nD", "D", "How many events observed in data in region D", Arg::Is::Required),

		// Output options
		Arg("OutputFile", "f", "Name of output root file for limit results. By default based on input filename", Arg::Is::Optional),

		// How is the limit set?
		Flag("UseAsym", "a", "Do asymtotic fit rather than using toys (toys are slow!)", Arg::Is::Optional),
		Flag("ExtrapAtEachLifetime", "l", "Refit limit at each lifetime point to take into account differing efficiencies at A, B, C and D", Arg::Is::Optional),
		Arg("RescaleSignal", "r", "Rescale the expected signal in region A to this number during limit setting", Arg::Is::Optional),
	});

	// Make sure we got all the command line arguments we need
	if (argc == 1 || !args.Parse(argc, argv)) {
		cout << args.Usage("PlotSingleLimit") << endl;
		throw runtime_error("Bad command line arguments - exiting");
	}

	// Fill up the config guy and return it.
	config result;
	result.extrapolate_filename = args.Get("extrapFile");

	result.observed_data.A = args.GetAsFloat("nA");
	result.observed_data.B = args.GetAsFloat("nB");
	result.observed_data.C = args.GetAsFloat("nC");
	result.observed_data.D = args.GetAsFloat("nD");

	result.limit_settings.useToys = !args.IsSet("UseAsym");
	result.limit_settings.scaleLimitByEfficiency = !args.IsSet("ExtrapAtEachLifetime");


	auto baseFilename = result.extrapolate_filename.rfind("\\") == string::npos
		? result.extrapolate_filename
		: result.extrapolate_filename.substr(result.extrapolate_filename.rfind("\\") + 1);
	result.limit_settings.fileName = args.IsSet("OutputFile")
		? args.Get("OutputFile")
		: string("limit_result_") + baseFilename;
	result.limit_settings.rescaleSignalTo = args.IsSet("RescaleSignal")
		? args.GetAsFloat("RescaleSignal")
		: 0.0;

	return result;
}