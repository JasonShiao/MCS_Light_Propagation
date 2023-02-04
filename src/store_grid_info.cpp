#include "store_grid_info.h"

#define TISSUE_WIDTH 0.3      // in cm
#define GRID_SIZE_R 0.01      // in cm
#define GRID_SIZE_Z 0.01      // in cm

using namespace std;

bool simGridInfo(std::string output_file){

   ofstream Grid_Out;
   Grid_Out.open(output_file.c_str());

   // Seed with generator
   mt19937 generator(time(0));
   uniform_real_distribution<double> distribution(0.0, 1.0);

   const double ANISO = 0.91;
   // Direction of Z-axis is Downward
   const double OBSERVE_BOUND_LEFT = -9999; // cm
   const double OBSERVE_BOUND_RIGHT = 9999; // cm
   const double OBSERVE_BOUND_TOP = -9999; // cm
   const double OBSERVE_BOUND_BOTTOM = 9999; // cm

   Media Air(0, 0, 1, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, OBSERVE_BOUND_TOP, 0);
   Media Tissue(6, 414, 1.37, OBSERVE_BOUND_LEFT, OBSERVE_BOUND_RIGHT, 0, 0.15);

   const int PHOTON_NUM = 10000;

   double w;
   double Refl = 0, Trans = 0, Absorb = 0;
   double cx, cy, cz;  // velocity direction
   double x, y, z;   // position
   double theta, phi;  // in radians
   double in_angle;    // in radians
   double theta_c = asin(Air.getN()/Tissue.getN());
   double StepSize;

   int Grid_rowCount = int((Tissue.getBound('B')-Tissue.getBound('T'))/GRID_SIZE_Z + 0.5); // round up
   int Grid_colCount = int(TISSUE_WIDTH/GRID_SIZE_R + 0.5);   // round up
   double** absorption_grid;
   absorption_grid = new double*[Grid_rowCount];
   for(int i = 0; i < Grid_rowCount; ++i){
      absorption_grid[i] = new double[Grid_colCount];
   }
   for(int i = 0; i < Grid_rowCount; ++i){
      for(int j = 0; j < Grid_colCount; ++j)
         absorption_grid[i][j] = 0;
   }
   //double absorption_grid[(int)((Tissue.getBound('B')-Tissue.getBound('T'))/GRID_SIZE_Z)][(int)(TISSUE_WIDTH/GRID_SIZE_R)] = {0};   // (z,r)


   cout<<"Simulating ..."<<endl;

   for(int i=0;i<PHOTON_NUM;i++){
      if( i % 200 == 0)
         cout << setw(2)<<setprecision(2)<< (int)((float)i/PHOTON_NUM*100) <<"%"<<"\r";
      //initialization
      int step_index = 0;
      cx = 0;cy = 0;cz = 1;   // cx^2+cy^2+cz^2 = 1
      // ============== Infinite Narrow =================
      //x = 0;y = 0;z = 0;
      z = 0;
      // ================== Uniform =====================
      /*double rBeam = 0.05*sqrt(distribution(generator));
      double thetaBeam = distribution(generator);
      x = rBeam*cos(2*PI*thetaBeam);
      y = rBeam*sin(2*PI*thetaBeam);*/
      // ================== Gaussian ====================
      double rBeam = 0.05*sqrt(-log(1-distribution(generator))/2);
      double thetaBeam = distribution(generator);
      x = rBeam*cos(2*PI*thetaBeam);
      y = rBeam*sin(2*PI*thetaBeam);

      w = 1;

      //incident from outside to tissue
      //if(distribution(generator) < pow((Tissue.getN()-Air.getN())/(Air.getN()+Tissue.getN()), 2 ) ){
      //   Refl += w;
      //   continue;
      //}
      Refl += w*pow((Tissue.getN()-Air.getN())/(Air.getN()+Tissue.getN()), 2 );
      w = 1 - w*pow((Tissue.getN()-Air.getN())/(Air.getN()+Tissue.getN()), 2 );

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
            //in_angle = acos(-cz);
            in_angle = -atan(sqrt(1-cz*cz)/cz);
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
            //in_angle = acos(cz);
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

         //Absorb += w*Tissue.getMuA()/Tissue.getMuT();
         // update Absorption grids
         int grid_z;
         int grid_r;
         grid_z = (int)(z/GRID_SIZE_Z);
         grid_r = (int)(sqrt(x*x+y*y)/GRID_SIZE_R);
         if( grid_z != 0 && z/GRID_SIZE_Z - grid_z == 0) // When between 2 z-grids
            grid_z--;
         if( grid_r != 0 && sqrt(x*x+y*y)/GRID_SIZE_R - grid_r == 0) // When between 2 r-grids
            grid_r--;
         if(grid_z < Grid_rowCount && grid_r < Grid_colCount)
            absorption_grid[grid_z][grid_r] += w*Tissue.getMuA()/Tissue.getMuT();

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
         }

         // get theta, phi
         phi = 2*PI*distribution(generator);
         theta = acos(1/(2*ANISO)*(1+pow(ANISO, 2.0)-pow((1-pow(ANISO, 2.0))/(1-ANISO+2*ANISO*distribution(generator)), 2.0)));
         // update cx, cy, cz
         double cx_t, cy_t, cz_t;
         if(cz > 0.99999 || cz < -0.99999){
            cx_t = sin(theta)*cos(phi);
            cy_t = sin(theta)*sin(phi);
            //cy_t = cz >= 0 ? sin(theta)*sin(phi):-sin(theta)*sin(phi);
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

   cout << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM << endl;
   //Grid_Out << "R = "<< (double)Refl/PHOTON_NUM <<", T = "<< (double)Trans/PHOTON_NUM <<endl;
   for(int i=0;i<Grid_rowCount;i++){
      for(int j=0;j<Grid_colCount;j++){
         Grid_Out << absorption_grid[i][j] << ", ";
      }
      Grid_Out << endl;
   }

   for(int i = 0; i < Grid_rowCount; ++i)
      delete [] absorption_grid[i];
   delete [] absorption_grid;
   Grid_Out.close();

   return true;

}
