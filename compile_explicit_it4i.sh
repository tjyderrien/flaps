#!/bin/bash

module load intel
module load mkl

export KMP_STACKSIZE=104857600
export OMP_STACKSIZE=104857600
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

# wihout OpenMP - without MKL - works
# ifort gmsh/libmsh2vf.f90 amos/*.f Bivariate.f zeroin.f control_file.f90 main_explicit_Mie.f90 \
# -traceback -mcmodel=large \
# -shared-intel -lm \
# -I/usr/include \
# -o Flaps_explicit.out 

# with iOMP5 - without MKL - runs, but on single core :-(
#ifort gmsh/libmsh2vf.f90 amos/*.f Bivariate.f zeroin.f control_file.f90 main_explicit_Mie.f90 \
#-traceback -mcmodel=large \
#-shared-intel -lm -liomp5 \
#-I/usr/include \
#-o Flaps_explicit.out 

# with OpenMP - without MKL - segmentation fault because of ulimit
ifort gmsh/libmsh2vf.f90 amos/*.f Bivariate.f zeroin.f control_file.f90 main_explicit_Mie.f90 \
-traceback -mcmodel=large \
-shared-intel -lm -openmp \
-I/usr/include \
-o Flaps_explicit.out 

# with IT4I optimizations, multithreaded lapack mkl - runs on one core only :-(
# ifort -ipo -O3 -vec -xAVX -vec-report1 -w gmsh/libmsh2vf.f90 amos/*.f Bivariate.f zeroin.f control_file.f90 main_explicit_Mie.f90 -traceback -mcmodel=large \
#-I/usr/include -I$MKL_INC_DIR -L$MKL_LIB_DIR \
#-shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -openmp -lpthread -W1 -lm -o Flaps_explicit.out

# Without optimizations, use INTEL compilers, use MKL library - problem in linking
# ifort -w gmsh/libmsh2vf.f90 amos/*.f Bivariate.f zeroin.f control_file.f90 main_explicit_Mie.f90 -traceback -shared-intel -mcmodel=large \
# -I/usr/include \
# -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm \
# -i8 -openmp -I${MKLROOT}/include \
# -o Flaps_explicit.out

# -I${MKL_INC_DIR} -L${MKL_LIB_DIR} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 \


#-I${MKL_INC_DIR} -L${MKL_LIB_DIR} \
#-lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -lpthread \




# using GNU compiler - not tested
# gfortran -w amos/*.f Bivariate.f zeroin.f main_explicit_Mie.f90 -O2 -xT -mcmodel=large -I/usr/include -openmp -lm -o Flaps_explicit.out

# ifort main.o Bivariate.o -o Flaps.out
#  -axSSE4.2 
# ifort main.f90 --O3 -lm -I/usr/include -o Flaps.out
# ./Flaps_explicit.out
