#!/bin/bash
echo "Loading FLAPS:OMP variables."
export OMP_STACKSIZE=10M
# export OMP_NUM_THREADS=8
# export MKL_NUM_THREADS=8
export MKL_DYNAMIC="FALSE"
export OMP_DYNAMIC="FALSE"
export OMP_SCHEDULE="DYNAMIC,1600"
ulimit -s unlimited

cp FermiDatasE.sav FermiDatasE.dat
cp FermiDatasH.sav FermiDatasH.dat

