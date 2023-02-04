#ifndef RAND_H_INCLUDED
#define RAND_H_INCLUDED

#include <iostream>
#include <stdlib.h>
#include <random>
#include <string>

class RandomNumber {
 public:
   RandomNumber() = delete;
   RandomNumber(double min, double max);

   double gen();
   void resetSeed();
   void test(int num_test, std::string output_file_raw, std::string output_file_interval);

 private:
   std::random_device rd_;
   std::mt19937 mt_;
   std::uniform_int_distribution<> dist_; // uniformly distributed int number on the interval  [min,max]
   // Q: Is uniform_real_distribution not spread evenly enough? 
   // Q: Is uniform_int_distribution spread more evenly?
   //std::uniform_real_distribution<double> dist_; // uniformly distributed double number on the interval  [min,max)

   double min_;
   double max_;
};

extern RandomNumber random_number;

//double genRand();
//bool testRand(std::string output_file_raw, std::string output_file_interval);

#endif // RAND_H_INCLUDED
