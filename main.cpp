#include <iostream>
#include <ctime>       /* time */
#include <stdlib.h>

#include "Rand.h"
#include "AbsorbOnly.h"

using namespace std;


int main()
{
   srand(time(NULL));
   //testRand();

   // Simulate light propagation in absorption-only media
   absorb_only();

   return 0;
}
