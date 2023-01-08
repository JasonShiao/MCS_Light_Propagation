#ifndef RAND_H_INCLUDED
#define RAND_H_INCLUDED

#include <iostream>
#include <ctime>       /* time */
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <random>
#include <string>

double genRand();
bool testRand(std::string output_file_raw, std::string output_file_interval);

#endif // RAND_H_INCLUDED
