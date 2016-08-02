//
// Code to run the lifetime limit setting extrapolations
//
#include "limitSetting.h"
#include "extrap_file_wrapper.h"
#include "limit_output_file.h"

#include "SimulABCD.h"

#include <algorithm>

using namespace std;

// Helper functions
limit_result rescale_limit_by_efficiency(const limit_result &original, const vector<double> &original_lifetime, const vector<double> &new_lifetime);

// Extrapolate vs lifetime by:
//  1. Calculate the limit where the sample was generated
//  2. Use the efficiency evolution as a function of lifetime to
//     determine the limit at other points.
void extrapolate_limit_to_lifetime_by_efficency(const extrap_file_wrapper &input,
	const ABCD &dataObserved,
	const abcd_limit_config &limit_params)
{
	// Get the limit at the generated lifetime.
	auto generated = input.generated_lifetime();

	// Now, get the limit for it.
	auto generated_limit = do_abcd_limit(dataObserved, generated.signalEvents, limit_params);

	// For each lifetime in the input, rescale the limit till we have a limit number for all
	// the outputs.
	auto lifetimes = input.list_of_lifetimes();
	vector<limit_result> results;
	transform(lifetimes.begin(), lifetimes.end(), back_inserter(results),
		[&generated_limit, &generated, &input](double ctau)
		{
			return rescale_limit_by_efficiency(generated_limit, generated.efficiency, input.lifetime(ctau).efficiency);
		}
	);

	// Great. Next, we have to build the output file and write this!
	write_limit_output_file("limits.root", results);
}

// Given a limit that was generated at a particular lifetime, rescale for the new lifetime.
// The beauty about doing this this way, is that signal strength, which is what the limit it, does not
// change. Just the expected signal.
limit_result rescale_limit_by_efficiency(const limit_result &original,
	const vector<double> &original_efficiency,
	const vector<double> &new_efficiency)
{
	auto result = original;
	result.expected_signal.A *= new_efficiency[A] / original_efficiency[A];
	result.expected_signal.B *= new_efficiency[B] / original_efficiency[B];
	result.expected_signal.C *= new_efficiency[C] / original_efficiency[C];
	result.expected_signal.D *= new_efficiency[D] / original_efficiency[D];

	return result;
}
