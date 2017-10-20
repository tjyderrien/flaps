#!/bin/bash
# DebugLine=-g -check all -warn all -backtrace 
DebugLine="" #-g -bounds-check -w -backtrace -Wno-unused-variable -Wall
Optimization="$(DebugLine) -O3 -I/usr/include -llapack -lm  -openmp" #-openmp #-mcmodel=large"

FC="ifort" #gfortran
FCFLAGS="${Optimization}"

make clean
autoreconf -i
./configure FC=${FC} FCLAGS=${FCFLAGS}
make 
# make install

cp src/flaps ./
