// Some dirt simple data structures that need to be shared, and so the defins are pulled out.
#ifndef __limit_datastructures__
#define __limit_datastructures__

#include <vector>

// Number of events, signal, etc., in the ABCD regions.
struct ABCD {
	double A, B, C, D;
};

// Allow use to refer to regions directly instead of by numbers.
enum Region {
	A = 0, B = 1, C = 2, D = 3
};

// For signal, the number of events is also associated with a lifetime.
struct signal_lifetime {
	ABCD signalEvents;
	double lifetime;
	std::vector<double> efficiency; // A, B, C, D efficiencies for the sample
};

// Control how the limit is actually run
struct abcd_limit_config {
	bool useToys; // True if we should run toys, otherwise run asym fit.
};

// Result (and input parameters) from a limit.
struct limit_result {
	ABCD observed_data;
	ABCD expected_signal;
	double cl_95;
	double cl_1sigma;
	double cl_2sigma;
};

#endif
