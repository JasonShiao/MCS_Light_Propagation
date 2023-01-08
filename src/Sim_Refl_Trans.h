#ifndef SIM_REFL_TRANS_H_INCLUDED
#define SIM_REFL_TRANS_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>

#include "Rand.h"
#include "MediaStructure.h"

extern const double PI;

bool sim_rt_fixed_weight_iso(std::string output_file);
bool sim_rt_var_weight_aniso(std::string output_file);

#endif // SIM_REFL_TRANS_H_INCLUDED
