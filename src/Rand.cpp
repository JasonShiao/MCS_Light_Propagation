#include "rand.h"

using namespace std;

double genRand(){
   return (double)(rand()%1000)/1000;
}

bool testRand(std::string output_file_raw, std::string output_file_interval){

   default_random_engine generator;
   uniform_real_distribution<double> distribution(0.0, 1.0);

   const int NumTest = 10000;

   ofstream classify, rawdata;
   classify.open(output_file_interval.c_str());
   rawdata.open(output_file_raw.c_str());
   double RandArray[NumTest];

   for(int i=0;i<10000;++i) {
      RandArray[i] = distribution(generator);
   }

   int interval[20]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};

   for(int i=0;i<10000;++i){
      RandArray[i]>0.95?(interval[19]++):RandArray[i]>0.9?(interval[18]++):RandArray[i]>0.85?(interval[17]++):RandArray[i]>0.8?(interval[16]++):\
      RandArray[i]>0.75?(interval[15]++):RandArray[i]>0.7?(interval[14]++):RandArray[i]>0.65?(interval[13]++):RandArray[i]>0.6?(interval[12]++):\
      RandArray[i]>0.55?(interval[11]++):RandArray[i]>0.5?(interval[10]++):RandArray[i]>0.45?(interval[9]++):RandArray[i]>0.4?(interval[8]++):\
      RandArray[i]>0.35?(interval[7]++):RandArray[i]>0.3?(interval[6]++):RandArray[i]>0.25?(interval[5]++):RandArray[i]>0.2?(interval[4]++):\
      RandArray[i]>0.15?(interval[3]++):RandArray[i]>0.1?(interval[2]++):RandArray[i]>0.05?(interval[1]++):RandArray[i]>0?(interval[0]++):(interval[0]++);

      rawdata << RandArray[i]<<", ";
      if(i%5==4) rawdata<<endl;
   }

   cout << "=========================================================" << endl;
   cout << "Statistic of one random sample with sample size of 10000:" << endl;

   for(int i=0;i<20;++i){
      cout<<"["<< i*0.05<<", "<<(i+1)*0.05<<"]: "<<interval[i]<<endl;
      classify << "["<< i*0.05<<", "<<(i+1)*0.05<<"]: "<<interval[i]<<endl;
   }
   cout<<endl;

   int sum=0;
   int sqr_sum=0;
   for(int i=0;i<20;i++){
      sum += interval[i];
      sqr_sum += interval[i]*interval[i];
   }
   double mean = (double)(sum)/20;

   cout << "mean: " << mean << ", std = " << sqrt((sqr_sum-20*mean*mean)/19) << endl;
   cout << "=========================================================" << endl;

   classify << endl;
   classify << "mean: " << mean << ", std = " << sqrt((sqr_sum-20*mean*mean)/19) << endl;
   rawdata << endl;

   return true;
}
