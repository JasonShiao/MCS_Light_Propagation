#include "reflection_and_transmission.h"
#include "random_number.h"
#include "medium_structure.h"

using namespace std;

/**
 * @brief Simulate a colliminated light propagating along the +z-axis
 *        in a tissue slab medium with u_s = 90 cm^-1, u_a = 10 cm^-1, thickness = 0.02 cm
 *        
 *        Simulate with "isotropic scattering" and 10,000 "fixed weight" photons
 * 
 *        Assume the index matched with the outside medium ( = no direct reflection in the boundary)
 *        Calculate the transmission fraction T and reflection fraction R,
 *        The solution of adding-doubling (by Van de Hulst) gives reflectance R = 0.3616 and transmittance T = 0.3565.
*/
bool simFixedWeightIsoScatter(std::string output_file){

   ofstream RT_Out;
   RT_Out.open(output_file.c_str(), ofstream::app);

   // Direction of Z-axis is Downward
   const int OBSERVE_BOUND_LEFT = -9999; // cm
   const int OBSERVE_BOUND_RIGHT = 9999; // cm
   const int OBSERVE_BOUND_TOP = -9999; // cm
   const int OBSERVE_BOUND_BOTTOM = 9999; // cm

   Medium tissue_slab_medium(10, 90, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.02);
   const int kTotalNumPhoton = 10000;

   int w = 1;  // fixed weight
   int reflection_fraction = 0, transmission_fraction = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;   // position
   double theta, phi;  // in radians
   double step_size;

   for (int i = 0; i < kTotalNumPhoton; i++) {
      int step_index = 0;
      cx = 0; cy = 0; cz = 1;
      x = 0; y = 0; z = 0;

      // Propagation of One photon
      while (1) {
         //get step size
         step_size = -log(random_number.gen()) / tissue_slab_medium.getMuT();

         // Move
         step_index++;
         //move photon
         x += step_size * cx;
         y += step_size * cy;
         z += step_size * cz;

         if (z <= tissue_slab_medium.getBound('T')) {
            reflection_fraction += 1.0;   // fixed weight photon
            break;
         } else if (z >= tissue_slab_medium.getBound('B')) {
            transmission_fraction += 1.0; // fixed weight photon
            break;
         } else {
            if (random_number.gen() < tissue_slab_medium.getMuA() / tissue_slab_medium.getMuT()){  // Absorbed in Medium
               break;
            } else { // Scattered
               // get theta, phi
               phi = 2 * PI * (random_number.gen());
               theta = acos(2 * random_number.gen() - 1);
               // update cx, cy, cz
               double cx_t, cy_t, cz_t;
               if (cz > 0.99999 || cz < -0.99999) {
                  cx_t = sin(theta) * cos(phi);
                  cy_t = sin(theta) * sin(phi);
                  cz_t = cz >= 0 ? cos(theta) : -cos(theta);
               } else {
                  cx_t = sin(theta) * (cx*cz*cos(phi) - cy*sin(phi)) / sqrt(1 - cz * cz) + cx * cos(theta);
                  cy_t = sin(theta) * (cy*cz*cos(phi) + cx*sin(phi)) / sqrt(1 - cz * cz) + cy * cos(theta);
                  cz_t = -sin(theta) * cos(phi) * sqrt(1 - cz * cz) + cz * cos(theta);
               }
               cx = cx_t;
               cy = cy_t;
               cz = cz_t;
            }
         }

      }

   }

   cout << "R = " << (double)reflection_fraction / kTotalNumPhoton << ", T = " << (double)transmission_fraction/kTotalNumPhoton << endl;
   RT_Out << "R = " << (double)reflection_fraction / kTotalNumPhoton << ", T = " << (double)transmission_fraction/kTotalNumPhoton << endl;
   RT_Out.close();

   return true;
}

/**
 * @brief Simulate a colliminated light propagating along the +z-axis
 *        in a tissue slab medium with u_s = 90 cm^-1, u_a = 10 cm^-1, thickness = 0.02 cm
 *        
 *        Simulate with "anisotropic scattering" and 10,000 "variable weight" photons
 *        Use the Henyey-Greenstein phase function with Anisotropy factor g = 0.75
 * 
 *        Assume the index matched with the outside medium ( = no direct reflection in the boundary)
 *        Calculate the transmission fraction T and reflection fraction R,
 *        The solution of adding-doubling (by Van de Hulst) gives reflectance R = 0.09739 and transmittance T = 0.66096.
*/
bool simVarWeightAnisoScatter(std::string output_file) {

   ofstream RT_Out;
   RT_Out.open(output_file.c_str(), ofstream::app);

   const double kAnisotropic = 0.75;

   // Direction of Z-axis is Downward
   const int OBSERVE_BOUND_LEFT = -9999; // cm
   const int OBSERVE_BOUND_RIGHT = 9999; // cm
   const int OBSERVE_BOUND_TOP = -9999; // cm
   const int OBSERVE_BOUND_BOTTOM = 9999; // cm

   Medium tissue_slab_medium(10, 90, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.02);
   const int kTotalNumPhoton = 10000;

   double w;  // fixed weight
   double reflection_fraction = 0, transmission_fraction = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;     // position
   double theta, phi;  // in radians
   double step_size;

   for (int i = 0; i < kTotalNumPhoton; i++) {
      w = 1 ;
      int step_index = 0;
      cx = 0; cy = 0; cz = 1;
      x = 0; y = 0; z = 0;

      // Propagation of One photon
      while(1){
         // 1. get step size
         step_size = -log(random_number.gen()) / tissue_slab_medium.getMuT();
         // 2. move photon
         x += step_size * cx;
         y += step_size * cy;
         z += step_size * cz;
         // 3. Check condition
         if (z <= tissue_slab_medium.getBound('T')) {
            reflection_fraction += w;
            break;
         } else if(z >= tissue_slab_medium.getBound('B')) {
            transmission_fraction += w;
            break;
         } else {
            // Update photon weight (decay in medium)
            w = w * tissue_slab_medium.getMuS() / tissue_slab_medium.getMuT();
            if (w < 0.001) { // weight too small
               // play roulette to determine whether the photon can go further
               if (random_number.gen() > 0.05) { // Absorbed
                  break;
               } else { // continue but with much smaller weight?
                  w *= 0.05;
                  continue;
               }
            } else {
               // get theta, phi
               phi = 2 * PI * random_number.gen();
               theta = acos(1 / (2 * kAnisotropic) * (1 + pow(kAnisotropic, 2.0) - pow((1 - pow(kAnisotropic, 2.0)) / (1 - kAnisotropic + 2 * kAnisotropic * random_number.gen()), 2.0)));
               // update cx, cy, cz
               double cx_t, cy_t, cz_t;
               if (cz > 0.99999 || cz < -0.99999) {
                  cx_t = sin(theta) * cos(phi);
                  cy_t = sin(theta) * sin(phi);
                  cz_t = cz >= 0 ? cos(theta) : -cos(theta);
               } else {
                  cx_t = sin(theta) * (cx * cz * cos(phi)-cy * sin(phi)) / sqrt(1 - cz * cz) + cx * cos(theta);
                  cy_t = sin(theta) * (cy * cz * cos(phi)+cx * sin(phi)) / sqrt(1 - cz * cz) + cy * cos(theta);
                  cz_t = - sin(theta) * cos(phi) * sqrt(1 - cz * cz) + cz * cos(theta);
               }
               cx = cx_t;
               cy = cy_t;
               cz = cz_t;
            }
         }
      } // end of while loop for propogation of a single photon
   } // end of for loop for iterating all photons

   cout << "R = " << reflection_fraction / kTotalNumPhoton << ", T = " << transmission_fraction/kTotalNumPhoton << endl;
   RT_Out << "R = " << reflection_fraction / kTotalNumPhoton << ", T = " << transmission_fraction/kTotalNumPhoton << endl;
   RT_Out.close();

   return true;

}
