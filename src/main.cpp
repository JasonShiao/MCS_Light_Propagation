#include <iostream>
#include <ctime>       /* time */
#include <stdlib.h>
#include <string>

#include "Rand.h"
#include "AbsorbOnly.h"
#include "Sim_Refl_Trans.h"
#include "Mismatch_Bound.h"
#include "Store_grid_info.h"

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
  srand(time(NULL));
  //testRand(output_file, "data/interval.txt");

  // 3. Start simulation
  // Simulate light propagation in absorption-only media
  //absorb_only(output_file); // "data/Distr_Out.csv"

  // Simulate Reflectance and Transmittance of Single-layer media
  //sim_rt_fixed_weight_iso(output_file); // "data/RT_fixed_weight_Out.txt"
  
  // Simulate Reflectance and Transmittance of Single-layer media
  //sim_rt_var_weight_aniso(output_file); // "data/RT_variable_weight_Out.txt"

  // Simulate Mismatch Boundary: Air-Tissue
  //sim_mismatch_boundary(output_file); // "data/Refl_Mismatch_Boundary_Out.txt"

  // Simulate and Store Absorption in grids
  // Infinite narrow beam
  Sim_gridInfo(output_file); // "data/Grid_Out.csv"
  // TODO: Gaussian beam

  // TODO: Uniform beam

  return 0;
}
