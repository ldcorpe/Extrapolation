// Contains the headers to run the various types of limit setting

#ifndef __limitSetting__
#define __limitSetting__

#include "limit_datastructures.h"

class extrap_file_wrapper;

// Run the limit where the MC file was generated, and then scale it.
void extrapolate_limit_to_lifetime_by_efficency(const extrap_file_wrapper &input,
	const ABCD &dataObserved,
	const abcd_limit_config &limit_params);

// Run the limit at each point in lifetime.
void extrapolate_limit_to_lifetime(const extrap_file_wrapper &input,
	const ABCD &dataObserved,
	const abcd_limit_config &limit_params);

#endif

