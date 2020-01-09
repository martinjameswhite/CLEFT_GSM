#!/bin/bash
#
g++ --std=c++11 -funroll-loops -fopenmp \
    -o zeldovich zeldovich.cpp -lfftw3 -lfftw3_omp -lm
#
export OMP_NUM_THREADS=8
#
./zeldovich pklin_RunPB.txt > test_cpp.txt
#
rm -f zeldovich
#
