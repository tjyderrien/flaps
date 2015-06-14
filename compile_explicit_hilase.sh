#!/bin/bash

# export KMP_STACKSIZE=1048576000
export OMP_STACKSIZE=10485760
export OMP_NUM_THREADS=16
# export MKL_NUM_THREADS=16
export MKL_DYNAMIC="FALSE"
export OMP_DYNAMIC="FALSE"
export OMP_SCHEDULE="DYNAMIC,800"
ulimit -s unlimited

KMP_AFFINITY=granularity=fine,compact

rm Flaps_explicit.out; rm *.o
cp FermiDatasE.sav FermiDatasE.dat
cp FermiDatasH.sav FermiDatasH.dat
# wihout OpenMP
# ifort main.f90 -traceback -mcmodel=large -shared-intel -lm -I/usr/include -o Flaps.out -axSSE4.2

# without OpenMP, with libraries
# ifort -w amos/*.f Bivariate.f zeroin.f matrixtools.f90 main_explicit_Mie.f90 -traceback -mcmodel=large -I/usr/include -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -shared-intel -llapack -lm -o Flaps_explicit.out

# with multithreaded lapack mkl
# ifort -w amos/*.f Bivariate.f zeroin.f matrixtools.f90 main_explicit_Mie.f90 -O2 -traceback -mcmodel=large -I/usr/include -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -R/opt/intel/mkl/lib/intel64 -shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -openmp -lpthread -W1 -lm -o Flaps_explicit.out

# using another compiler
cd amos/
gfortran -c *.f
cd ..
gfortran -c Bivariate.f zeroin.f
gfortran -c -g -w -fbacktrace gmsh/libmsh2vf.f90 -ffree-line-length-none -ffixed-line-length-none -finteger-4-integer-8
gfortran -c main_explicit_Mie.f90 -g -Wno-unused-variable -Wall -fbacktrace -fbounds-check -O2 -mcmodel=small -I/usr/include -llapack -fopenmp -ffree-line-length-none -ffixed-line-length-none -finteger-4-integer-8 -lm
gfortran main_explicit_Mie.o -fbacktrace -fopenmp -mcmodel=small amos/*.o Bivariate.o zeroin.o libmsh2vf.o -o Flaps_explicit.out 
./Flaps_explicit.out

# ifort main.o Bivariate.o -o Flaps.out
#  -axSSE4.2 
# ifort main.f90 --O3 -lm -I/usr/include -o Flaps.out
# ./Flaps_explicit.out
