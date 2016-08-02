//
// Given an extrapolation file as input, derive the limits as a function of
// lifetime for output.
//

#include "Wild/CommandLine.h"
#include "TApplication.h"
#include "TROOT.h"

#include <iostream>
#include <memory>

using namespace std;
using namespace Wild::CommandLine;

// Helper functions, etc.
struct config {

};
config parse_command_line(int argc, char **argv);

// Main entry point
int main(int argc, char **argv)
{
	// Get ROOT setup so we don't have funny problems with it later.
	int dummy_argc = 0;
	auto a = make_unique<TApplication>("extrapolate_betaw", &dummy_argc, argv);
	gROOT->SetBatch(true);

	// Protect against anything weird going wrong so we get a sensible error message.
	try {
		cout << "hi" << endl;
		auto c = parse_command_line(argc, argv);
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
		// Arguments that must be present on the command line, one after the other
		Arg("extrapFile", "e", "The extrapolation root file", Arg::Is::Required),

		// Options
		Flag("UseAsym", "a", "Do asymtotic fit rather than using toys (toys are slow!)", Arg::Is::Optional)
	});

	// Make sure we got all the command line arguments we need
	if (argc == 1 || !args.Parse(argc, argv)) {
		cout << args.Usage("PlotSingleLimit") << endl;
		throw runtime_error("Bad command line arguments - exiting");
	}
}