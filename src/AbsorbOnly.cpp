#include "AbsorbOnly.h"

using namespace std;

bool absorb_only(std::string output_file){

   ofstream Distr_Out;
   Distr_Out.open(output_file.c_str());

   // Seed with generator
   mt19937 generator(time(0));
   uniform_real_distribution<double> distribution(0.0, 1.0);
   // Direction of Z-axis is Downward
   const int OBSERVE_BOUND_LEFT = -9999; // cm
   const int OBSERVE_BOUND_RIGHT = 9999; // cm
   const int OBSERVE_BOUND_TOP = -9999; // cm
   const int OBSERVE_BOUND_BOTTOM = 9999; // cm

   Media SimpleStruct(10, 0, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 1);
   const int NumPhoton = 10000;
   const double step_size = 0.025; // cm

   int NumInterval = int(SimpleStruct.get_length()/step_size+0.5);

   vector<int> absorb_distr( NumInterval, 0);

   for(int i=0;i<NumPhoton;i++){

      int step_index = 0;
      double pos_z = 0;

      // Propagation of One photon
      while(1){
         // Move
         step_index++;
         pos_z += step_size;
         // if random number < absorption probability -> Absorb
         if( distribution(generator) < SimpleStruct.get_mu_a()*step_size ){
            absorb_distr[step_index-1] =  absorb_distr[step_index-1]+1;
            break;
         }

      }

   }

   for(vector<int>::iterator it = absorb_distr.begin() ; it != absorb_distr.end(); ++it){
      Distr_Out << *it << ", " << (double)(*it)/NumPhoton<<endl;
   }

   return true;
}
