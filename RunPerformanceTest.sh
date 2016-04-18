#!/bin/bash

export OMP_STACKSIZE=1048576
# export OMP_NUM_THREADS=4
# export MKL_NUM_THREADS=8
export MKL_DYNAMIC="FALSE"
export OMP_DYNAMIC="FALSE"
# export OMP_SCHEDULE="DYNAMIC,800"
export OMP_SCHEDULE="STATIC"
export KMP_AFFINITY=compact
ulimit -s unlimited

for i in 4 2 1
do
	export OMP_NUM_THREADS=$i
	OutputFile="M2001-N151-${i}cores.log"
#	OutputFile2="${OutputFile}.2.log"
	echo ${OutputFile}
	time ./Flaps_explicit.out >& `echo ${OutputFile}`
	echo $i > /dev/null #This is just to flush variables.
done
