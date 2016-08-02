// Contains methods for moveing between the A, B, C, D regions in CalR and LJ,
// which are defined slightly differently and thus the ratios are slightly off.
// it was thought this was easier than modifying the workspace they use
// (and that we copy).

#ifndef __CalRLJConverter__
#define __CalRLJConverter__

#include <vector>

inline std::vector<double> ConvertFromCalRToLJ(const double calRRegions[4])
{
	std::vector<double> result;

	// regions as defined in our note: A = B*C/D. Here it's defined A = B*D/C so changing C <--> D as def in the note
	result.push_back(calRRegions[0]); // A
	result.push_back(calRRegions[1]); // B
	result.push_back(calRRegions[3]); // D
	result.push_back(calRRegions[2]); // C

	return result;
}

#endif