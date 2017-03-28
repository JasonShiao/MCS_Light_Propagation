#include <iostream>
#include <ctime>       /* time */
#include <stdlib.h>

#include "Rand.h"
#include "AbsorbOnly.h"
#include "Sim_Refl_Trans.h"
#include "Mismatch_Bound.h"
#include "Store_grid_info.h"

using namespace std;

const double PI = 3.14159265;

int main()
{
   srand(time(NULL));
   //testRand();

   // Simulate light propagation in absorption-only media
   //absorb_only();

   // Simulate Reflectance and Transmittance of Single-layer media
   //sim_rt_fixed_weight_iso();

   // Simulate Reflectance and Transmittance of Single-layer media
   //sim_rt_var_weight_aniso();

   // Simulate Mismatch Boundary: Air-Tissue
   //sim_mismatch_boundary();

   // Simulate and Store Absorption in grids
   // Infinite narrow beam
   Sim_gridInfo();
   // Gaussian beam

   // Uniform beam

   return 0;
}
