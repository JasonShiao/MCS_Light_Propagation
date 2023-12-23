#ifndef STORE_GRID_INFO_H_INCLUDED
#define STORE_GRID_INFO_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

#include "random_number.h"
#include "medium_structure.h"

extern const double PI;

enum class InputBeamType {
    kInifiniteNarrowBeam,
    kUniformDistribBeam,
    kGaussianDistribBeam
};

bool simAbsorptionGridInfo(std::string output_file, InputBeamType input_beam_type);


#endif // STORE_GRID_INFO_H_INCLUDED
