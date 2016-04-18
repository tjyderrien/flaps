#!/bin/bash

# Converts MAP file of Depth.dat into POS file for GMSH

file=`ls -1 -t *.vid | sort -n | tail -n1`
awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$4}' ${file}  > Intensity.map
awk '{ if($5!="") print }' mesh.dat > mesh.tmp

# Now we have same number of elements in each files
wc -l Intensity.map 
wc -l mesh.tmp

awk '{ print $5, $1*1E6, $2*1E6, $2*0.0 }' mesh.tmp > Intensity.0
awk '{ print $5 }' mesh.tmp > Intensity.2
awk '{ print $5, 15,1, $5*1E6 }' mesh.tmp > Intensity.4
awk '{ print $3*1E4 }' Intensity.map > Intensity.3
paste Intensity.2 Intensity.3 > Intensity.1
rm Intensity.3

# Now files for construction of POS file are ready. 
wc -l Intensity.0
wc -l Intensity.1 
