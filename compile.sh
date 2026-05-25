#!/bin/bash
# DebugLine="-g -fcheck=bounds -fbacktrace -Wno-unused-variable -Wall"
Optimization="-O3 -fopenmp"

FC=gfortran
FCFLAGS="${Optimization}"

module load openblas

OPENBLAS_LIBS=$(pkg-config --libs openblas)
OPENBLAS_LIBDIR=$(pkg-config --variable=libdir openblas)

make clean
autoreconf -i
./configure FC=${FC} FCFLAGS="${FCFLAGS}" \
  LDFLAGS="-Wl,-rpath,${OPENBLAS_LIBDIR}" \
  --with-blas="${OPENBLAS_LIBS}" \
  --with-lapack="${OPENBLAS_LIBS}"
make
make install



