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
limit_result rescale_limit_by_efficiency(const limit_result &original, const signal_lifetime &original_lifetime, const signal_lifetime &new_lifetime);

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
	auto generated_limit = do_abcd_limit(dataObserved, generated, limit_params);

	// For each lifetime in the input, rescale the limit till we have a limit number for all
	// the outputs.
	auto lifetimes = input.list_of_lifetimes();
	vector<limit_result> results;
	transform(lifetimes.begin(), lifetimes.end(), back_inserter(results),
		[&generated_limit, &generated, &input](double ctau)
		{
			return rescale_limit_by_efficiency(generated_limit, generated, input.lifetime(ctau));
		}
	);

	// Great. Next, we have to build the output file and write this!
	write_limit_output_file(limit_params, results, input.get_ctau_binning());
}

// Extrapolate vs lifetime by:
//  1. Calculate the limit at each lifetime
//  2. Use the efficiency evolution as a function of lifetime to
//     determine the limit at other points.
void extrapolate_limit_to_lifetime(const extrap_file_wrapper &input,
	const ABCD &dataObserved,
	const abcd_limit_config &limit_params)
{
	// For each lifetime in the input, rescale the limit till we have a limit number for all
	// the outputs.
	auto lifetimes = input.list_of_lifetimes();
	vector<limit_result> results;
	transform(lifetimes.begin(), lifetimes.end(), back_inserter(results),
		[&input, &dataObserved, &limit_params](double ctau)
	{
		// Get the signal info and do the lmit.
		auto signal = input.lifetime(ctau);
		return do_abcd_limit(dataObserved, signal, limit_params);
	}
	);

	// Great. Next, we have to build the output file and write this!
	write_limit_output_file(limit_params, results, input.get_ctau_binning());
}

// Given a limit that was generated at a particular lifetime, rescale for the new lifetime.
limit_result rescale_limit_by_efficiency(const limit_result &original,
	const signal_lifetime &original_efficiency,
	const signal_lifetime &new_efficiency)
{
	// Update the number of events expected
	auto result = original;
	result.signal = new_efficiency;

	// Next we have to update the signal strength. Mu represents the number of events
	// from signal that are compatible with observed data and errors. So that number
	// of events is what will remain constant. Rescale mu keeping that in mind.
	result.cl_95 *= original.signal.signalEvents.A / result.signal.signalEvents.A;
	result.cl_p1sigma *= original.signal.signalEvents.A / result.signal.signalEvents.A;
	result.cl_p2sigma *= original.signal.signalEvents.A / result.signal.signalEvents.A;
	result.cl_n1sigma *= original.signal.signalEvents.A / result.signal.signalEvents.A;
	result.cl_n2sigma *= original.signal.signalEvents.A / result.signal.signalEvents.A;

	return result;
}
