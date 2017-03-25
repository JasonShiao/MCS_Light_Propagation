#include "AbsorbOnly.h"

using namespace std;

bool absorb_only(){

   ofstream Distr_Out;
   Distr_Out.open("Data/Distr_Out.csv");

   Media SimpleStruct;
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
         if( genRand() < SimpleStruct.get_mu_a(pos_z)*step_size ){
            absorb_distr[step_index-1] =  absorb_distr[step_index-1]+1;
            break;
         }

      }

   }

   for(vector<int>::iterator it = absorb_distr.begin() ; it != absorb_distr.end(); ++it){
      Distr_Out << *it << ", " << (double)(*it)/NumPhoton<<endl;
   }

}
