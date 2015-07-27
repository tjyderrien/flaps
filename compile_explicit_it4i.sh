#!/bin/bash

# export KMP_STACKSIZE=104857600
# export OMP_STACKSIZE=104857600
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export MKL_DYNAMIC="FALSE"
export OMP_DYNAMIC="FALSE"
export OMP_SCHEDULE="DYNAMIC,800"
# ulimit -s unlimited

KMP_AFFINITY=granularity=fine,compact

rm Flaps_explicit.out
cp FermiDatasE.sav FermiDatasE.dat
cp FermiDatasH.sav FermiDatasH.dat
# wihout OpenMP
# ifort main.f90 -traceback -mcmodel=large -shared-intel -lm -I/usr/include -o Flaps.out -axSSE4.2

# without OpenMP, with libraries
# ifort -w amos/*.f Bivariate.f zeroin.f matrixtools.f90 main_explicit_Mie.f90 -traceback -mcmodel=large -I/usr/include -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -shared-intel -llapack -lm -o Flaps_explicit.out

# with multithreaded lapack mkl 
ifort -ipo -O3 -vec -xAVX -vec-report1 -w gmsh/libmsh2vf.f90 amos/*.f Bivariate.f zeroin.f control_file.f90 main_explicit_Mie.f90 -traceback -mcmodel=large \
-I/usr/include -I$MKL_INC_DIR -L$MKL_LIB_DIR \
-shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -liomp5 -lpthread -W1 -lm -o Flaps_explicit.out

# using another compiler
# gfortran -w amos/*.f Bivariate.f zeroin.f main_explicit_Mie.f90 -O2 -xT -mcmodel=large -I/usr/include -openmp -lm -o Flaps_explicit.out

# ifort main.o Bivariate.o -o Flaps.out
#  -axSSE4.2 
# ifort main.f90 --O3 -lm -I/usr/include -o Flaps.out
# ./Flaps_explicit.out
