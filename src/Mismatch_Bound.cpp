#include "Mismatch_Bound.h"


using namespace std;

bool sim_mismatch_boundary(std::string output_file){

   ofstream Refl_Out;
   Refl_Out.open(output_file.c_str(), ofstream::app);

   // Seed with generator
   mt19937 generator(time(0));
   uniform_real_distribution<double> distribution(0.0, 1.0);

   // Direction of Z-axis is Downward
   const double OBSERVE_BOUND_LEFT = -9999; // cm
   const double OBSERVE_BOUND_RIGHT = 9999; // cm
   const double OBSERVE_BOUND_TOP = -9999; // cm
   const double OBSERVE_BOUND_BOTTOM = 9999; // cm

   Media Air(0, 0, 1, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, OBSERVE_BOUND_TOP, 0);
   Media Tissue(10, 90, 1.5, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.2);

   const int PHOTON_NUM = 10000;

   double w;
   double Refl = 0, Trans = 0, Absorb = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;   // position
   double theta, phi;  // in radians
   double in_angle;    // in radians
   double theta_c = asin(Air.get_n()/Tissue.get_n());
   double StepSize;

   for(int i=0;i<PHOTON_NUM;i++){

      int step_index = 0;
      cx = 0;cy = 0;cz = 1;   // cx^2+cy^2+cz^2 = 1
      x = 0;y = 0;z = 0;
      w = 1;

      //incident from outside to tissue
      if(distribution(generator) < pow((Tissue.get_n()-Air.get_n())/(Air.get_n()+Tissue.get_n()), 2 ) ){
         Refl += w;
         continue;
      }
      //Refl += w*pow((Tissue.get_n()-Air.get_n())/(Air.get_n()+Tissue.get_n()), 2 );
      //w = 1 - w*pow((Tissue.get_n()-Air.get_n())/(Air.get_n()+Tissue.get_n()), 2 );

      // Propagation of One photon
      while(1){
         //get step size
         StepSize = - log(distribution(generator))/Tissue.get_mu_t();

         double r, t, Rs, Rp;
         //move photon
         step_index++;
         x=x+StepSize*cx;
         y=y+StepSize*cy;
         z=z+StepSize*cz;

         // In tissue ?
         if( z < Tissue.get_bound('T') ){
            //in_angle = -atan(sqrt(1-cz*cz)/cz);
            in_angle = acos(-cz);
            // TIR ?
            if( in_angle > theta_c ){
               r=1;
               t=0;
            }else{
               double cos_t = sqrt( 1 - pow(Tissue.get_n()*sin(in_angle)/Air.get_n(),2) );
               if ( pow( Tissue.get_n()*sin(in_angle)/Air.get_n(), 2) > 0.999999)
                  cos_t = 0;
               Rs = pow( ( ( Tissue.get_n()*cos(in_angle)-Air.get_n()*cos_t ) / ( Tissue.get_n()*cos(in_angle)+Air.get_n()*cos_t ) ), 2);
               Rp = pow( ( ( Tissue.get_n()*cos_t-Air.get_n()*cos(in_angle) ) / ( Tissue.get_n()*cos_t+Air.get_n()*cos(in_angle) ) ), 2);
               r = (Rs+Rp)/2;
               t = 1-r;
            }

            if(distribution(generator) < t){
               Refl = Refl + w;
               break;
            }

            z = Tissue.get_bound('T') - z;
            cz = -cz;

         }else if( z > Tissue.get_bound('B') ){
            //in_angle = atan(sqrt(1-cz*cz)/cz);
            in_angle = acos(cz);
            // TIR ?
            if( in_angle > theta_c ){
               r=1;
               t=0;
            }else{
               double cos_t = sqrt( 1 - pow(Tissue.get_n()*sin(in_angle)/Air.get_n(),2) );
               if ( pow( Tissue.get_n()*sin(in_angle)/Air.get_n(), 2) > 0.999999)
                  cos_t = 0;
               Rs = pow( ( ( Tissue.get_n()*cos(in_angle)-Air.get_n()*cos_t ) / ( Tissue.get_n()*cos(in_angle)+Air.get_n()*cos_t ) ), 2);
               Rp = pow( ( ( Tissue.get_n()*cos_t-Air.get_n()*cos(in_angle) ) / ( Tissue.get_n()*cos_t+Air.get_n()*cos(in_angle) ) ), 2);
               r = (Rs+Rp)/2;
               t = 1-r;
            }

            if(distribution(generator) < t){
               Trans = Trans + w;
               break;
            }

            z = 2*Tissue.get_bound('B') - z;
            cz = -cz;

         }

         // absorption
         Absorb += w*Tissue.get_mu_a()/Tissue.get_mu_t();
         // update weight
         w = w*Tissue.get_mu_s()/Tissue.get_mu_t();
         if( w < 0.001){
            // play roulette
            if(distribution(generator) > 0.05){
               //Absorb+=w;
               break;
            }else{
               w*=20;
               continue;
            }
         }else {
            // get theta, phi
            phi = 2*PI*(distribution(generator));
            theta = acos(2*distribution(generator)-1);
            // update cx, cy, cz
            double cx_t, cy_t, cz_t;
            if(cz > 0.99999 || cz < -0.99999){
               cx_t = sin(theta)*cos(phi);
               cy_t = cz >= 0 ? sin(theta)*sin(phi):-sin(theta)*sin(phi);
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

   cout << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM << ", A = " << (double)Absorb/PHOTON_NUM << endl;
   cout<<"R + T + A = "<<(double)(Refl+Trans+Absorb)/PHOTON_NUM<<endl;
   Refl_Out << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM <<endl;
   Refl_Out.close();

   return true;

}
