#ifndef SIM_REFL_TRANS_H_INCLUDED
#define SIM_REFL_TRANS_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>

#include "rand.h"
#include "media_structure.h"

extern const double PI;

bool simRTFixedWeightIso(std::string output_file);
bool simRTVarWeightAnIso(std::string output_file);

#endif // SIM_REFL_TRANS_H_INCLUDED
