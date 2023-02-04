#ifndef MEDIUMSTRUCTURE_H_INCLUDED
#define MEDIUMSTRUCTURE_H_INCLUDED

#include <iostream>
using namespace std;

class Medium{
   public:
      Medium();
      Medium(double mu_absorb, double mu_scatter, double b_left, double b_right, double b_top, double b_bottom);
      Medium(double mu_absorb, double mu_scatter, double n, double b_left, double b_right, double b_top, double b_bottom);

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

      virtual ~Medium();
   protected:
   private:
      double mu_a;           // in unit of cm^(-1)
      double mu_s;           // no unit
      // Z-axis Downward
      double left_boundary;  // in unit of cm
      double right_boundary; // in unit of cm
      double top_boundary;   // in unit of cm
      double bottom_boundary;// in unit of cm
      double diffrac_coeff;

};

#endif // MEDIUMSTRUCTURE_H_INCLUDED
