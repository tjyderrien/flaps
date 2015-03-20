#!/bin/bash

# export KMP_STACKSIZE=1048576000
export OMP_NUM_THREADS=4
ulimit -s unlimited
cp FermiDatasE.sav FermiDatasE.dat
cp FermiDatasH.sav FermiDatasH.dat

# wihout lapack
ifort main.f90 -traceback -heap-arrays -shared-intel -lm -openmp -I/usr/include -o Flaps.out
# ifort main.f90 -O3 -lm -I/usr/include -o Flaps.out

# with lapack
# ifort -lm -I/usr/include -lapack -openmp -o Flaps.out 

