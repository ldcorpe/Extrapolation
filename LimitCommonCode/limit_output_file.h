// Write out a limit file
#ifndef __limit_output_file__
#define __limit_output_file__

#include "limit_datastructures.h"

#include <string>
#include <vector>

void write_limit_output_file(const std::string &fname, const std::vector<limit_result> &results);

#endif
