#ifndef MEDIASTRUCTURE_H_INCLUDED
#define MEDIASTRUCTURE_H_INCLUDED

#include <iostream>
using namespace std;

class Media{
   public:
      Media();
      Media(double, double, double, double, double, double);
      Media(double, double, double, double, double, double, double);

      double get_mu_a(){ return mu_a; };
      double get_mu_s(){ return mu_s; };
      double get_mu_t(){ return mu_a+mu_s; };
      double get_n(){ return diffrac_coeff; };

      double get_length(){ return bottom_boundary - top_boundary; };
      double get_bound( char side) {
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
