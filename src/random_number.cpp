#include "random_number.h"

#include <fstream>
#include <cmath>
#include <vector>
//#include <ctime>       /* time */

using namespace std;

RandomNumber random_number{0.0, 1.0};

RandomNumber::RandomNumber(double min, double max): mt_(rd_()), min_{min}, max_{max} {
   //dist_ = std::uniform_real_distribution<>(min, max);
   dist_ = std::uniform_int_distribution<>(int(min * 10000), int(max * 10000));
}

void RandomNumber::resetSeed() {
   mt_ = std::mt19937(rd_());
}

double RandomNumber::gen(){
   // return dist_(mt_); // for real_uniform_distribution
   return (double)(dist_(mt_)) / 10000;
}

void RandomNumber::test(int num_test, std::string output_file_raw, std::string output_file_interval){

   ofstream classify, rawdata;
   classify.open(output_file_interval.c_str());
   rawdata.open(output_file_raw.c_str());

   //const int kNumTest = 10000;

   std::vector<double> rand_array;
   for(auto i = 0; i < num_test; ++i) {
      rand_array.push_back(this->gen());
   }
   
   const int kIntervalCnt = 20;
   int interval[kIntervalCnt] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};

   size_t idx = 0;
   for (const auto & num: rand_array) {
      ++interval[int((num - min_) * kIntervalCnt / (max_ - min_)) >= kIntervalCnt ? 
                     (kIntervalCnt - 1) : 
                     int((num - min_) * kIntervalCnt / (max_ - min_)) ];
      //num > 0.95 ? (interval[19]++) : num > 0.9 ? (interval[18]++) : num>0.85 ? (interval[17]++) : num>0.8 ? (interval[16]++):\
      //num > 0.75 ? (interval[15]++) : num > 0.7 ? (interval[14]++) : num>0.65 ? (interval[13]++) : num>0.6 ? (interval[12]++):\
      //num > 0.55 ? (interval[11]++) : num > 0.5 ? (interval[10]++) : num>0.45 ? (interval[9]++)  : num>0.4 ? (interval[8]++):\
      //num > 0.35 ? (interval[7]++)  : num > 0.3 ? (interval[6]++)  : num>0.25 ? (interval[5]++)  : num>0.2 ? (interval[4]++):\
      //num > 0.15 ? (interval[3]++)  : num > 0.1 ? (interval[2]++)  : num>0.05 ? (interval[1]++)  : num>0   ? (interval[0]++) : (interval[0]++);

      rawdata << num << ", ";
      idx++;
      if (idx % 5 == 0) rawdata << endl;
   }

   cout << "=========================================================" << endl;
   cout << "Statistic of one random sample with sample size of 10000:" << endl;

   for(int i = 0; i < kIntervalCnt; ++i){
      cout << "[" << min_ + i * (max_ - min_) / kIntervalCnt << ", " 
                  << min_ + (i + 1 ) * (max_ - min_) / kIntervalCnt <<  "]: " << interval[i] << endl;
      classify << "[" << min_ + i * (max_ - min_) / kIntervalCnt << ", "
                      << min_ + (i + 1) * (max_ - min_) / kIntervalCnt << "]: " << interval[i] << endl;
   }
   cout << endl;
}
