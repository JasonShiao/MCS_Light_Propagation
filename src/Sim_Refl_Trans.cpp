#include "Sim_Refl_Trans.h"

using namespace std;

bool sim_rt_fixed_weight_iso(){

   ofstream RT_Out;
   RT_Out.open("Data/RT_fixed_weight_Out.txt", ofstream::app);

   // Seed with generator
   mt19937 generator(time(0));
   uniform_real_distribution<double> distribution(0.0, 1.0);
   // Direction of Z-axis is Downward
   const int OBSERVE_BOUND_LEFT = -9999; // cm
   const int OBSERVE_BOUND_RIGHT = 9999; // cm
   const int OBSERVE_BOUND_TOP = -9999; // cm
   const int OBSERVE_BOUND_BOTTOM = 9999; // cm

   Media SimpleStruct(10, 90, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.02);
   const int PHOTON_NUM = 10000;

   int w = 1;  // fixed weight
   int Refl = 0, Trans = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;   // position
   double theta, phi;  // in radians
   double StepSize;

   for(int i=0;i<PHOTON_NUM;i++){

      int step_index = 0;
      cx = 0;cy = 0;cz = 1;
      x = 0;y = 0;z = 0;

      // Propagation of One photon
      while(1){
         //get step size
         StepSize = - log(distribution(generator))/SimpleStruct.get_mu_t();

         // Move
         step_index++;
         //move photon
         x=x+StepSize*cx;
         y=y+StepSize*cy;
         z=z+StepSize*cz;

         if( z <= SimpleStruct.get_bound('T')){
            Refl++;
            break;
         }else if( z >= SimpleStruct.get_bound('B')){
            Trans++;
            break;
         }else{
            if(distribution(generator) < SimpleStruct.get_mu_a()/SimpleStruct.get_mu_t()){  // Absorbed in Media
               break;
            }else {     // Scattered
               // get theta, phi
               phi = 2*PI*(distribution(generator));
               theta = acos(2*distribution(generator)-1);
               // update cx, cy, cz
               double cx_t, cy_t, cz_t;
               if(cz > 0.99999 || cz < -0.99999){
                  cx_t = sin(theta)*cos(phi);
                  cy_t = sin(theta)*sin(phi);
                  cz_t = cz >= 0 ? cos(theta):-cos(theta);
               }else{
                  cx_t = sin(theta)*(cx*cz*cos(phi)-cy*sin(phi))/sqrt(1-cz*cz) + cx*cos(theta);
                  cy_t = sin(theta)*(cy*cz*cos(phi)+cx*sin(phi))/sqrt(1-cz*cz) + cy*cos(theta);
                  cz_t = - sin(theta)*cos(phi)*sqrt(1-cz*cz) + cz*cos(theta);
               }
               cx = cx_t;
               cy = cy_t;
               cz = cz_t;
            }
         }

      }

   }

   cout << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM <<endl;
   RT_Out << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM <<endl;
   RT_Out.close();

   return true;


}


bool sim_rt_var_weight_aniso(){

   ofstream RT_Out;
   RT_Out.open("Data/RT_variable_weight_Out.txt", ofstream::app);

   const double ANISO = 0.75;

   // Seed with generator
   mt19937 generator(time(0));
   uniform_real_distribution<double> distribution(0.0, 1.0);
   // Direction of Z-axis is Downward
   const int OBSERVE_BOUND_LEFT = -9999; // cm
   const int OBSERVE_BOUND_RIGHT = 9999; // cm
   const int OBSERVE_BOUND_TOP = -9999; // cm
   const int OBSERVE_BOUND_BOTTOM = 9999; // cm

   Media SimpleStruct(10, 90, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.02);
   const int PHOTON_NUM = 10000;

   double w;  // fixed weight
   double Refl = 0, Trans = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;   // position
   double theta, phi;  // in radians
   double StepSize;

   for(int i=0;i<PHOTON_NUM;i++){

      w = 1 ;
      int step_index = 0;
      cx = 0;cy = 0;cz = 1;
      x = 0;y = 0;z = 0;

      // Propagation of One photon
      while(1){
         //get step size
         StepSize = - log(distribution(generator))/SimpleStruct.get_mu_t();
         //move photon
         x=x+StepSize*cx;
         y=y+StepSize*cy;
         z=z+StepSize*cz;
         if( z <= SimpleStruct.get_bound('T')){
            Refl += w;
            break;
         }else if( z >= SimpleStruct.get_bound('B') ){
            Trans += w;
            break;
         }else{
            // update weight
            w = w*SimpleStruct.get_mu_s()/SimpleStruct.get_mu_t();
            if( w < 0.001){
               // play roulette
               if(distribution(generator) > 0.05)
                  break;
               else{
                  w*=0.05;
                  continue;
               }
            }else {
               // get theta, phi
               phi = 2*PI*distribution(generator);
               theta = acos(1/(2*ANISO)*(1+pow(ANISO, 2.0)-pow((1-pow(ANISO, 2.0))/(1-ANISO+2*ANISO*distribution(generator)), 2.0)));
               // update cx, cy, cz
               double cx_t, cy_t, cz_t;
               if(cz > 0.99999 || cz < -0.99999){
                  cx_t = sin(theta)*cos(phi);
                  cy_t = sin(theta)*sin(phi);
                  cz_t = cz >= 0 ? cos(theta):-cos(theta);
               }else{
                  cx_t = sin(theta)*(cx*cz*cos(phi)-cy*sin(phi))/sqrt(1-cz*cz) + cx*cos(theta);
                  cy_t = sin(theta)*(cy*cz*cos(phi)+cx*sin(phi))/sqrt(1-cz*cz) + cy*cos(theta);
                  cz_t = - sin(theta)*cos(phi)*sqrt(1-cz*cz) + cz*cos(theta);
               }
               cx = cx_t;
               cy = cy_t;
               cz = cz_t;
            }
         }
      }

   }

   cout << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM <<endl;
   RT_Out << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM <<endl;
   RT_Out.close();

   return true;

}
