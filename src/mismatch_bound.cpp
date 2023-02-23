#include "mismatch_bound.h"
#include "random_number.h"
#include "medium_structure.h"


using namespace std;

/**
 * @brief Simulate with mismatched boundary between the tissue and the outside medium. 
 *        Assume a semi-infinite slab (make d at least 20 optical depths (>0.2cm?)) with 
 *        mu_a = 10 cm-1, mu_s = 90 cm-1 and "isotropic scattering". 
 *        Collinimate beam? or unpolarized beam?
 *        Let n1 (air) = 1, n2 (tissue) = 1.5.
 * 
 *        Determine R using 10,000 "variable weight" photons.
 *        Giovanelli has reported a theoretical value of R = 0.2600.
 * 
 *        Fresnel equations: describe the reflection and transmission of light (or electromagnetic radiation 
 *                           in general) when incident on an interface between different optical media.
*/
bool simMismatchBoundary(std::string output_file){

   ofstream result_output;
   result_output.open(output_file.c_str(), ofstream::app);

   // Direction of Z-axis is Downward
   const double OBSERVE_BOUND_LEFT = -9999; // cm
   const double OBSERVE_BOUND_RIGHT = 9999; // cm
   const double OBSERVE_BOUND_TOP = -9999; // cm
   const double OBSERVE_BOUND_BOTTOM = 9999; // cm

   Medium air_medium(0, 0, 1, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, OBSERVE_BOUND_TOP, 0);
   Medium tissue_medium(10, 90, 1.5, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.2);

   const int kTotalNumPhoton = 10000;

   double w;
   double reflection_amount = 0, transimission_amount = 0, absorption_amount = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;     // position
   double theta, phi;  // in radians
   double in_angle;    // in radians
   double theta_c = asin(air_medium.getN() / tissue_medium.getN());
   double step_size;

   for (int i = 0; i < kTotalNumPhoton; i++){
      // 1. Initialize photon
      int step_index = 0;
      cx = 0; cy = 0; cz = 1;   // cx^2+cy^2+cz^2 = 1
      x = 0; y = 0; z = 0;
      w = 1;

      // 2. Incident from outside (air) to tissue
      // Approach 1: All reflect or all transmit with probability that follows Fresnel equation, 
      if (random_number.gen() < pow((tissue_medium.getN() - air_medium.getN()) / (air_medium.getN() + tissue_medium.getN()), 2 )) { // Fresnel equations for normal incident
         reflection_amount += w;
         continue;
      }
      // Approach 2: Partially reflect with ratio that follows Fresnel equation
      //reflection_amount += w * pow((tissue_medium.getN() - air_medium.getN()) / (air_medium.getN() + tissue_medium.getN()), 2);
      //w = 1 - w * pow((tissue_medium.getN() - air_medium.getN()) / (air_medium.getN() + tissue_medium.getN()), 2);

      // Propagation of One photon
      while (1) {
         // 1. get step size
         step_size = -log(random_number.gen()) / tissue_medium.getMuT();

         // 2. move photon
         double r, t, Rs, Rp;
         step_index++;
         x += step_size * cx;
         y += step_size * cy;
         z += step_size * cz;

         // Check condition (whether still in tissue)
         if (z < tissue_medium.getBound('T')) { // From tissue to air through top boundary
            //in_angle = -atan(sqrt(1 - cz * cz) / cz);
            in_angle = acos(-cz);
            // Check Total Internal Reflection (TIR) condition
            if (in_angle > theta_c) {
               r = 1;
               t = 0;
            } else { // Partial reflection
               double cos_t = sqrt( 1 - pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2));
               if (pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2) > 0.999999)
                  cos_t = 0;
               Rs = pow( ( ( tissue_medium.getN() * cos(in_angle) - air_medium.getN() * cos_t ) / ( tissue_medium.getN() * cos(in_angle) + air_medium.getN() * cos_t ) ), 2);
               Rp = pow( ( ( tissue_medium.getN() * cos_t - air_medium.getN() * cos(in_angle) ) / ( tissue_medium.getN() * cos_t + air_medium.getN() * cos(in_angle) ) ), 2);
               r = (Rs + Rp) / 2;
               t = 1 - r;
            }

            if (random_number.gen() < t) {
               reflection_amount = reflection_amount + w;
               break;
            }

            z = tissue_medium.getBound('T') - z;
            cz = -cz; // bounce back from the boundary
         } else if (z > tissue_medium.getBound('B')) { // From tissue to air through bottom boundary
            //in_angle = atan(sqrt(1 - cz * cz) / cz);
            in_angle = acos(cz);
            // Check Total Internal Reflection (TIR) condition
            if (in_angle > theta_c) {
               r = 1;
               t = 0;
            } else {  // Partial reflection
               double cos_t = sqrt(1 - pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2) );
               if (pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2) > 0.999999)
                  cos_t = 0;
               Rs = pow(
                        ( (tissue_medium.getN() * cos(in_angle) - air_medium.getN() * cos_t) / (tissue_medium.getN() * cos(in_angle) + air_medium.getN() * cos_t) ), 
                        2
                     );
               Rp = pow( 
                        ( (tissue_medium.getN() * cos_t - air_medium.getN() * cos(in_angle)) / (tissue_medium.getN() * cos_t + air_medium.getN() * cos(in_angle)) ), 
                        2
                     );
               r = (Rs + Rp) / 2;
               t = 1 - r;
            }

            if (random_number.gen() < t) {
               transimission_amount = transimission_amount + w;
               break;
            }

            z = 2 * tissue_medium.getBound('B') - z;
            cz = -cz; // bounce back from the boundary
         }

         // absorption
         absorption_amount += w * tissue_medium.getMuA() / tissue_medium.getMuT();
         // update weight
         w = w * tissue_medium.getMuS() / tissue_medium.getMuT();
         if (w < 0.001) {
            // play roulette
            if (random_number.gen() > 0.05) {
               //absorption_amount+=w;
               break;
            } else{
               w *= 20;
               continue;
            }
         } else {
            // get theta, phi
            phi = 2 * PI * (random_number.gen());
            theta = acos(2 * random_number.gen() - 1);
            // update cx, cy, cz
            double cx_t, cy_t, cz_t;
            if (cz > 0.99999 || cz < -0.99999) {
               cx_t = sin(theta) * cos(phi);
               cy_t = cz >= 0 ? sin(theta) * sin(phi) : -sin(theta) * sin(phi);
               cz_t = cz >= 0 ? cos(theta) : -cos(theta);
            } else {
               cx_t = sin(theta) * (cx * cz * cos(phi) - cy * sin(phi)) / sqrt(1 - cz * cz) + cx * cos(theta);
               cy_t = sin(theta) * (cy * cz * cos(phi) + cx * sin(phi)) / sqrt(1 - cz * cz) + cy * cos(theta);
               cz_t = - sin(theta) * cos(phi) * sqrt(1 - cz * cz) + cz * cos(theta);
            }
            cx = cx_t;
            cy = cy_t;
            cz = cz_t;
         }

      } // end of while loop for propagation of a single photon
   } // end of for loop iteration for all photons

   cout <<  "R = " << reflection_amount / kTotalNumPhoton << ", T = " << transimission_amount / kTotalNumPhoton << ", A = " << absorption_amount/kTotalNumPhoton << endl;
   cout << "R + T + A = " << (reflection_amount + transimission_amount+absorption_amount) / kTotalNumPhoton << endl;
   result_output << "R = " << reflection_amount / kTotalNumPhoton << ", T = "<< transimission_amount / kTotalNumPhoton << endl;
   result_output.close();

   return true;

}
