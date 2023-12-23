# Monte-Carlo Simulation for Light Propagation in media


## Run CMake
```cmake -S . -B build```

## Build
```cmake --build build```

## Notes
1. When we reduce n for rand()%n/n, 
the std of number of data between 20 intervals will decrease!

If we need a finer rand with acceptable std -> rand()*5%10000/10000

2. When step size decreases, the data will be closer to the curve generated from Beer-Lambert law
