#include <iostream>
#include <ctime>       /* time */
#include <stdlib.h>
#include <string>

#include "random_number.h"
#include "absorb_only.h"
#include "reflection_and_transmission.h"
#include "mismatch_bound.h"
#include "store_grid_info.h"

#include <argparse/argparse.hpp>

using namespace std;

const double PI = 3.14159265;

int main(int argc, char *argv[])
{
  // 1. Parsing input arguments
  argparse::ArgumentParser program("MCSLightPropagation");
  program.add_argument("-o", "--output")
    .required()
    .help("Filepath and name for the output data file");

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  auto output_file = program.get<std::string>("-o");
  std::cout << "Output filepath: " << output_file << std::endl;

  // 2. Setup
  // 2.5 Test (Verify) random number
  //random_number.test(10000, output_file, output_file.append("_classify"));

  // 3. Start simulation
  // Simulate light propagation in absorption-only media
  //absorbOnly(output_file); // "data/Distr_Out.csv"

  // Simulate Reflectance and Transmittance of Single-layer media
  //simRTFixedWeightIso(output_file); // "data/RT_fixed_weight_Out.txt"
  
  // Simulate Reflectance and Transmittance of Single-layer media
  //simRTVarWeightAnIso(output_file); // "data/RT_variable_weight_Out.txt"

  // Simulate Mismatch Boundary: Air-Tissue
  //simMismatchBoundary(output_file); // "data/Refl_Mismatch_Boundary_Out.txt"

  // Simulate and Store Absorption in grids
  // Infinite narrow beam
  //simGridInfo(output_file); // "data/Grid_Out.csv"
  // TODO: Gaussian beam

  // TODO: Uniform beam

  return 0;
}
