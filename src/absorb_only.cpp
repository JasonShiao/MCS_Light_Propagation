#include "absorb_only.h"
#include "random_number.h"
#include "medium_structure.h"

using namespace std;

/**
 * @brief Simulate a collimated (parallel) beam propagating along the +z axis 
 *        in an absorbing medium that is 1.0 cm thick
*/
bool absorbOnly(std::string output_file){

   ofstream distrib_output;
   distrib_output.open(output_file.c_str());

   // NOTE: Direction of +z-axis is downward
   const int OBSERVE_BOUND_LEFT = -9999; // cm
   const int OBSERVE_BOUND_RIGHT = 9999; // cm
   const int OBSERVE_BOUND_TOP = -9999;  // cm
   const int OBSERVE_BOUND_BOTTOM = 9999;// cm

   Medium plate_struct(10, 0, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 1);
   const int kTotalNumPhoton = 10000;
   //const double kStepSize = 0.025;    // The simulation result with this value will not be very close to Beer-Lambert Law)
   const double kStepSize = 0.001;      // cm (This value will heavily affect the result, smaller value -> closer to Beer-Lambert Law)
   const double kDepthInterval = kStepSize; // cm, Q: Why = kStepSize? 

   int num_depth_interval = int(plate_struct.getLength() / kDepthInterval + 0.5);

   vector<int> absorb_distrib(num_depth_interval, 0); // store the number of photon in each depth interval

   for (int i = 0; i < kTotalNumPhoton; i++) {

      int step_index = 0;
      double pos_z = 0;

      // Propagation of One photon
      while (1) {
         // 1. Move
         step_index++;
         pos_z += kStepSize; // Since there is no scattering in absorb only medium, step size = move distance in +z-axis
         // 2. Check if random number < absorption probability, if so, absorb, otherwise, continue
         if (random_number.gen() < plate_struct.getMuA() * kStepSize ){
            absorb_distrib[int(step_index * kStepSize / kDepthInterval - 1)] += 1;
            break;
         }
      }

   }

   for (auto & num_photon_in_interval: absorb_distrib) {
      distrib_output << num_photon_in_interval << ", " << (double)(num_photon_in_interval) / kTotalNumPhoton << endl;
   }

   // TODO: Compare with Beer-Lambert Law: I = I_0 * e^(-u_s * dist)

   return true;
}
