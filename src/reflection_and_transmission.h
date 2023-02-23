#ifndef SIM_REFL_TRANS_H_INCLUDED
#define SIM_REFL_TRANS_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>

#include "random_number.h"

extern const double PI;

bool simFixedWeightIsoScatter(std::string output_file);
bool simVarWeightAnisoScatter(std::string output_file);

#endif // SIM_REFL_TRANS_H_INCLUDED
