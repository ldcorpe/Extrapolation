// Make access to items in an extrapolation file easy so we don't have to pepper
// our code with weird ROOT calls.

#ifndef __extrap_file_wrapper__
#define __extrap_file_wrapper__

#include "limit_datastructures.h"

#include <TFile.h>

#include <string>
#include <memory>
#include <vector>

class extrap_file_wrapper {
public:
	inline extrap_file_wrapper(const std::string &filename)
	{
		_file = std::unique_ptr<TFile>(TFile::Open(filename.c_str(), "READ"));
		if (!_file->IsOpen()) {
			throw std::runtime_error((std::string("Unable to locate file ") + filename + " - can't continue").c_str());
		}
	}

	// Return the generated lifetime and A, B, C, and D for this MC input.
	inline signal_lifetime generated_lifetime() const;

	// Return a list of the lifetimes that we have efficiencies and #'s of events for
	inline std::vector<double> list_of_lifetimes() const;

	// Return what the signal sample woudl look like at a particular point
	inline signal_lifetime lifetime(double ctau) const;
private:
	std::unique_ptr<TFile> _file;
};

#endif
