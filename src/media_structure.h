#ifndef MEDIASTRUCTURE_H_INCLUDED
#define MEDIASTRUCTURE_H_INCLUDED

#include <iostream>
using namespace std;

class Media{
   public:
      Media();
      Media(double, double, double, double, double, double);
      Media(double, double, double, double, double, double, double);

      double getMuA(){ return mu_a; };
      double getMuS(){ return mu_s; };
      double getMuT(){ return mu_a+mu_s; };
      double getN(){ return diffrac_coeff; };

      double getLength(){ return bottom_boundary - top_boundary; };
      double getBound( char side) {
         switch(side){
            case 'L': return left_boundary;
            case 'R': return right_boundary;
            case 'T': return top_boundary;
            case 'B': return bottom_boundary;
            default: cout << side << " Boundary Not Defined" << endl ;return 0;
         }
      };

      virtual ~Media();
   protected:
   private:
      double mu_a;
      double mu_s;
      // Z-axis Downward
      double left_boundary;
      double right_boundary;
      double top_boundary;
      double bottom_boundary;
      double diffrac_coeff;

};

#endif // MEDIASTRUCTURE_H_INCLUDED
