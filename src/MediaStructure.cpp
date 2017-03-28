#include "MediaStructure.h"

Media::Media(){
   mu_a = 0;
   mu_s = 0;
   left_boundary = 0;
   right_boundary = 0;
   top_boundary = 0;
   bottom_boundary = 0;
   diffrac_coeff = 1;
}

Media::Media(double mu_absorb, double mu_scatter, double b_left, double b_right, double b_top, double b_bottom){
   mu_a = mu_absorb;
   mu_s = mu_scatter;
   left_boundary = b_left;
   right_boundary = b_right;
   top_boundary = b_top;
   bottom_boundary = b_bottom;
   diffrac_coeff = 1;
}

Media::Media(double mu_absorb, double mu_scatter, double n, double b_left, double b_right, double b_top, double b_bottom){
   mu_a = mu_absorb;
   mu_s = mu_scatter;
   left_boundary = b_left;
   right_boundary = b_right;
   top_boundary = b_top;
   bottom_boundary = b_bottom;
   diffrac_coeff = n;
}

Media::~Media()
{

}
