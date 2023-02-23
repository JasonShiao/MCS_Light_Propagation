#include "store_grid_info.h"
#include "medium_structure.h"
#include "random_number.h"

#define TISSUE_WIDTH 0.3      // in cm
#define GRID_SIZE_R 0.01      // in cm
#define GRID_SIZE_Z 0.01      // in cm

using namespace std;

/**
 * @brief Simulate with mu_a = 6 cm-1, mu_s = 414 cm-1, g = 0.91, Henyey-Greenstein phase function, 
 *        n1(air) = 1, n2(tissue) = 1.37, tissue thickness = 1.5 mm.
 *        Develop a grid (3 mm wide, 1.5 mm deep). 
 * 
 *        Let delta_r = delta_z = 0.1 mm. 
 *        Use an infinitely narrow beam with normal incidence at the origin
*/
bool simAbsorptionGridInfo(std::string output_file){

   ofstream absorb_grid_output;
   absorb_grid_output.open(output_file.c_str());

   const double kAnisotropic = 0.91;
   // Direction of Z-axis is Downward
   const double OBSERVE_BOUND_LEFT = -9999; // cm
   const double OBSERVE_BOUND_RIGHT = 9999; // cm
   const double OBSERVE_BOUND_TOP = -9999; // cm
   const double OBSERVE_BOUND_BOTTOM = 9999; // cm

   Medium air_medium(0, 0, 1, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, OBSERVE_BOUND_TOP, 0);
   Medium tissue_medium(6, 414, 1.37, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.15);

   const int kTotalNumPhoton = 10000;

   double reflection = 0, transmission = 0, absorption = 0;

   // Initialize storage grid
   int grid_row_count = int((tissue_medium.getBound('B') - tissue_medium.getBound('T')) / GRID_SIZE_Z + 0.5); // round up
   int grid_col_count = int(TISSUE_WIDTH/GRID_SIZE_R + 0.5);   // round up
   double** absorption_grid;
   absorption_grid = new double*[grid_row_count];
   for (int i = 0; i < grid_row_count; ++i) {
      absorption_grid[i] = new double[grid_col_count];
   }
   for (int i = 0; i < grid_row_count; ++i) {
      for(int j = 0; j < grid_col_count; ++j)
         absorption_grid[i][j] = 0;
   }
   //double absorption_grid[(int)((tissue_medium.getBound('B')-tissue_medium.getBound('T'))/GRID_SIZE_Z)][(int)(TISSUE_WIDTH/GRID_SIZE_R)] = {0};   // (z,r)

   cout << "Simulating ..." << endl;
   for (int i = 0; i < kTotalNumPhoton; i++) {
      double w;
      double x, y, z;     // position
      double cx, cy, cz;  // velocity direction
      double theta, phi;  // scattering, in unit of radians
      double theta_c = asin(air_medium.getN() / tissue_medium.getN());
      double step_size;

      if (i % 200 == 0) {
         cout << setw(2) << setprecision(2) << (int)((float)i / kTotalNumPhoton * 100) << "%" << "\r";
      }

      // Single photon initialization
      int step_index = 0;
      cx = 0; cy = 0; cz = 1;   // cx^2+cy^2+cz^2 = 1
      // ============== Infinite Narrow beam =================
      //x = 0;y = 0;z = 0;
      z = 0;
      // ================== Uniform distrib beam =====================
      /*double rBeam = 0.05*sqrt(random_number.gen());
      double thetaBeam = random_number.gen();
      x = rBeam*cos(2*PI*thetaBeam);
      y = rBeam*sin(2*PI*thetaBeam);*/
      // ================== Gaussian distrib beam ====================
      double rBeam = 0.05*sqrt(-log(1-random_number.gen())/2);
      double thetaBeam = random_number.gen();
      x = rBeam*cos(2*PI*thetaBeam);
      y = rBeam*sin(2*PI*thetaBeam);

      w = 1;

      //incident from outside to tissue
      //if(random_number.gen() < pow((tissue_medium.getN()-air_medium.getN())/(air_medium.getN()+tissue_medium.getN()), 2 ) ){
      //   reflection += w;
      //   continue;
      //}
      reflection += w*pow((tissue_medium.getN()-air_medium.getN())/(air_medium.getN()+tissue_medium.getN()), 2 );
      w = 1 - w*pow((tissue_medium.getN()-air_medium.getN())/(air_medium.getN()+tissue_medium.getN()), 2 );

      // Propagation of One photon
      while (1) {
         // 1. Get step size
         step_size = - log(random_number.gen())/tissue_medium.getMuT();

         // 2. Move photon
         step_index++;
         x += step_size * cx;
         y += step_size * cy;
         z += step_size * cz;

         // 3. Check whether still in tissue, handle reflection and transmission at boundary
         double fresnel_reflect, fresnel_transmit, fresnell_reflect_s_polar, fresnell_reflect_p_polar;
         double in_angle;    // in radians
         if (z < tissue_medium.getBound('T')) {
            //in_angle = acos(-cz);
            in_angle = -atan(sqrt(1 - cz * cz) / cz);
            // Check whether Total Internal Reflection (TIR)
            if (in_angle > theta_c) {
               fresnel_reflect = 1;
               fresnel_transmit = 0;
            } else {
               double cos_t = sqrt( 1 - pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2) );
               if (pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2) > 0.999999) {
                  cos_t = 0;
               }
               fresnell_reflect_s_polar = pow( ((tissue_medium.getN() * cos(in_angle) - air_medium.getN() * cos_t) / (tissue_medium.getN() * cos(in_angle) + air_medium.getN() * cos_t)), 2);
               fresnell_reflect_p_polar = pow( ((tissue_medium.getN() * cos_t - air_medium.getN() * cos(in_angle)) / (tissue_medium.getN() * cos_t + air_medium.getN() * cos(in_angle))), 2);
               fresnel_reflect = (fresnell_reflect_s_polar + fresnell_reflect_p_polar) / 2;
               fresnel_transmit = 1 - fresnel_reflect;
            }

            if (random_number.gen() < fresnel_transmit) {
               reflection = reflection + w;
               break;
            }

            z = tissue_medium.getBound('T') - z;
            cz = -cz;
         } else if (z > tissue_medium.getBound('B')) {
            //in_angle = acos(cz);
            // TIR ?
            if (in_angle > theta_c ) {
               fresnel_reflect = 1;
               fresnel_transmit = 0;
            } else {
               double cos_t = sqrt(1 - pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2));
               if (pow(tissue_medium.getN() * sin(in_angle) / air_medium.getN(), 2) > 0.999999) {
                  cos_t = 0;
               }
               fresnell_reflect_s_polar = pow( ( ( tissue_medium.getN()*cos(in_angle)-air_medium.getN()*cos_t ) / ( tissue_medium.getN()*cos(in_angle)+air_medium.getN()*cos_t ) ), 2);
               fresnell_reflect_p_polar = pow( ( ( tissue_medium.getN()*cos_t-air_medium.getN()*cos(in_angle) ) / ( tissue_medium.getN()*cos_t+air_medium.getN()*cos(in_angle) ) ), 2);
               fresnel_reflect = (fresnell_reflect_s_polar+fresnell_reflect_p_polar) / 2;
               fresnel_transmit = 1 - fresnel_reflect;
            }

            if (random_number.gen() < fresnel_transmit) {
               transmission = transmission + w;
               break;
            }

            z = 2 * tissue_medium.getBound('B') - z;
            cz = -cz;
         }

         //absorption += w*tissue_medium.getMuA()/tissue_medium.getMuT();
         // 4. Absorption, update absorption grids
         int grid_z;
         int grid_r;
         grid_z = (int)(z / GRID_SIZE_Z);
         grid_r = (int)(sqrt(x * x + y * y) / GRID_SIZE_R);
         if (grid_z != 0 && z / GRID_SIZE_Z - grid_z == 0) // When between 2 z-grids
            grid_z--;
         if (grid_r != 0 && sqrt(x * x + y * y) / GRID_SIZE_R - grid_r == 0) // When between 2 r-grids
            grid_r--;
         if (grid_z < grid_row_count && grid_r < grid_col_count)
            absorption_grid[grid_z][grid_r] += w * tissue_medium.getMuA() / tissue_medium.getMuT();

         // 5. Update weight after absorption (remaining scattering part)
         w = w * tissue_medium.getMuS() / tissue_medium.getMuT();
         if (w < 0.001) {
            // play roulette
            if (random_number.gen() > 0.05) {
               //absorption+=w;
               break;
            } else {
               w *= 20;
               continue;
            }
         }

         // 6. Anisotropic Scattering: get theta & phi and update cx, cy, cz
         phi = 2 * PI * random_number.gen();
         theta = acos(1 / (2 * kAnisotropic) * (1 + pow(kAnisotropic, 2.0) - pow((1 - pow(kAnisotropic, 2.0)) / (1 - kAnisotropic + 2 * kAnisotropic * random_number.gen()), 2.0)));
         double cx_t, cy_t, cz_t;
         if (cz > 0.99999 || cz < -0.99999) {
            cx_t = sin(theta) * cos(phi);
            cy_t = sin(theta) * sin(phi);
            //cy_t = cz >= 0 ? sin(theta)*sin(phi):-sin(theta)*sin(phi);
            cz_t = cz >= 0 ? cos(theta) : -cos(theta);
         } else {
            cx_t = sin(theta) * (cx * cz * cos(phi) - cy * sin(phi)) / sqrt(1 - cz * cz) + cx * cos(theta);
            cy_t = sin(theta) * (cy * cz * cos(phi) + cx * sin(phi)) / sqrt(1 - cz * cz) + cy * cos(theta);
            cz_t = -sin(theta) * cos(phi) * sqrt(1 - cz * cz) + cz * cos(theta);
         }
         cx = cx_t;
         cy = cy_t;
         cz = cz_t;
      } // end of while loop of propagation of a single photon
   } // end of for loop of iteration for all photons

   // Display result
   cout << "R = " << (double)reflection / kTotalNumPhoton << ", T = " << (double)transmission / kTotalNumPhoton << endl;
   //absorb_grid_output << "R = "<< (double)reflection/kTotalNumPhoton <<", T = "<< (double)transmission/kTotalNumPhoton <<endl;
   for (int i = 0; i < grid_row_count; i++) {
      for (int j = 0; j < grid_col_count; j++) {
         absorb_grid_output << absorption_grid[i][j] << ", ";
      }
      absorb_grid_output << endl;
   }

   // Clean up
   for (int i = 0; i < grid_row_count; ++i) {
      delete [] absorption_grid[i];
   }
   delete [] absorption_grid;
   absorb_grid_output.close();

   return true;
}
