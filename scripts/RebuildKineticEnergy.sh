#!/bin/bash
# extract the data 
# awk -f 'Scripts/video/BuildMaps.awk' Depth.dat 
# create a file with CellVols
awk '{ print $17 }' meshElements.dat > CellVols.tmp
numline=`wc -l CellVols.tmp | awk '{ print $1 }'` 
# let numline=numline-1 #6 lines were missed by BuildMaps.awk (??)
head -n${numline} CellVols.tmp > CellVols.dat
# for each created vid file, paste the CellVol data
for i in `ls *.vid -1 | sort -g`
do
	# filter data with data for energy conservation
	awk '{ print $1, $2, $3, $5, $6, $7, $8, $9 }' $i > $i.tmp
	paste $i.tmp CellVols.dat > $i.vol
	# t, x, y, Te, Th, Ts, Ne, Nh, CellVol
	# 1, 2, 3, 4 , 5 , 6 , 7 , 8 , 9
	#### calculate potential energy, kinetic energy of electrons, kinetic energy of lattice
	# t, pot-e, pot-h, kinetic e, kin-h, kin-s
	# echo $i.vol
	awk 'BEGIN {kb=1.38e-23; n0=5E28; ec=1.6e-19; gap=1.12*ec; pote=0e0; poth=0e0; kine=0e0; kinh=0e0; kins=0e0 } {
	time=$1
	kine=kine+1.5e0*kb*$7*$4*$9;
	kinh=kinh+1.5e0*kb*$8*$5*$9;
	pote=pote+gap*$7*$9
	poth=poth+gap*$8*$9
	kins=kins+kb*n0*$6*$9
	} END { print time, kine, kinh, kins, pote, poth }' $i.vol
done
