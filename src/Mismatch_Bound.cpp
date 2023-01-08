#include "mismatch_bound.h"


using namespace std;

bool simMismatchBoundary(std::string output_file){

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
   double theta_c = asin(Air.getN()/Tissue.getN());
   double StepSize;

   for(int i=0;i<PHOTON_NUM;i++){

      int step_index = 0;
      cx = 0;cy = 0;cz = 1;   // cx^2+cy^2+cz^2 = 1
      x = 0;y = 0;z = 0;
      w = 1;

      //incident from outside to tissue
      if(distribution(generator) < pow((Tissue.getN()-Air.getN())/(Air.getN()+Tissue.getN()), 2 ) ){
         Refl += w;
         continue;
      }
      //Refl += w*pow((Tissue.getN()-Air.getN())/(Air.getN()+Tissue.getN()), 2 );
      //w = 1 - w*pow((Tissue.getN()-Air.getN())/(Air.getN()+Tissue.getN()), 2 );

      // Propagation of One photon
      while(1){
         //get step size
         StepSize = - log(distribution(generator))/Tissue.getMuT();

         double r, t, Rs, Rp;
         //move photon
         step_index++;
         x=x+StepSize*cx;
         y=y+StepSize*cy;
         z=z+StepSize*cz;

         // In tissue ?
         if( z < Tissue.getBound('T') ){
            //in_angle = -atan(sqrt(1-cz*cz)/cz);
            in_angle = acos(-cz);
            // TIR ?
            if( in_angle > theta_c ){
               r=1;
               t=0;
            }else{
               double cos_t = sqrt( 1 - pow(Tissue.getN()*sin(in_angle)/Air.getN(),2) );
               if ( pow( Tissue.getN()*sin(in_angle)/Air.getN(), 2) > 0.999999)
                  cos_t = 0;
               Rs = pow( ( ( Tissue.getN()*cos(in_angle)-Air.getN()*cos_t ) / ( Tissue.getN()*cos(in_angle)+Air.getN()*cos_t ) ), 2);
               Rp = pow( ( ( Tissue.getN()*cos_t-Air.getN()*cos(in_angle) ) / ( Tissue.getN()*cos_t+Air.getN()*cos(in_angle) ) ), 2);
               r = (Rs+Rp)/2;
               t = 1-r;
            }

            if(distribution(generator) < t){
               Refl = Refl + w;
               break;
            }

            z = Tissue.getBound('T') - z;
            cz = -cz;

         }else if( z > Tissue.getBound('B') ){
            //in_angle = atan(sqrt(1-cz*cz)/cz);
            in_angle = acos(cz);
            // TIR ?
            if( in_angle > theta_c ){
               r=1;
               t=0;
            }else{
               double cos_t = sqrt( 1 - pow(Tissue.getN()*sin(in_angle)/Air.getN(),2) );
               if ( pow( Tissue.getN()*sin(in_angle)/Air.getN(), 2) > 0.999999)
                  cos_t = 0;
               Rs = pow( ( ( Tissue.getN()*cos(in_angle)-Air.getN()*cos_t ) / ( Tissue.getN()*cos(in_angle)+Air.getN()*cos_t ) ), 2);
               Rp = pow( ( ( Tissue.getN()*cos_t-Air.getN()*cos(in_angle) ) / ( Tissue.getN()*cos_t+Air.getN()*cos(in_angle) ) ), 2);
               r = (Rs+Rp)/2;
               t = 1-r;
            }

            if(distribution(generator) < t){
               Trans = Trans + w;
               break;
            }

            z = 2*Tissue.getBound('B') - z;
            cz = -cz;

         }

         // absorption
         Absorb += w*Tissue.getMuA()/Tissue.getMuT();
         // update weight
         w = w*Tissue.getMuS()/Tissue.getMuT();
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
