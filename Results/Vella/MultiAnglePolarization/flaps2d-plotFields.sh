#!/bin/bash

# This file collects all the Field distributions calculated by FLAPS and return a chain which can be used to crate a Tar package file. 

# FlapsRoot='/home/thibault/Documents/Codes/Flaps/20150428-flaps-gmsh/flaps2d' #Budecska
FlapsRoot='/home/tderrien/Codes/flaps2d' #IT4I Salomon
Files=''
for i in `du --max-depth=1 | awk '{ print $2 }'`
do
  echo [Directory] $i
  cd $i
  pwd
  gnuplot $FlapsRoot/Scripts/plot/plotTipField.gnu
  Files="${Files} ${i}/Field.dat"
  FilesEps="${FilesEps} ${i}/Field.eps"
  cd -
done

echo ${Files} ${FilesEps}
tar -cvf 'FieldMaps.tar' ${FilesEps} ${Files}
bzip2 'FieldMaps.tar'
