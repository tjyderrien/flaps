#!/bin/bash
# DebugLine=-g -check all -warn all -backtrace 
DebugLine="" #-g -bounds-check -w -backtrace -Wno-unused-variable -Wall
Optimization="$(DebugLine) -O3 -I/usr/include -llapack -lm  -openmp" #-openmp #-mcmodel=large"

C="ifort" #gfortran
CFLAGS="${Optimization}"

make clean
autoreconf -i
make 
make install



