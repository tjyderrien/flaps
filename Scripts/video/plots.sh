#! /bin/bash

valeurX=2.40024e-09
# valeurX=1.5557e-09
valeurY=0

# 1.00026e-08
# 1.60011e-09; 
# 2.50017e-9

valeurT=50e-15

# ## MAKE file for T,X profiles
for i in `ls Bonse*.dat`
do
 	echo "Filtering XZ $i"
#	awk "{ if(\$4==${valeurY}) print }" $i > $i.2DX &
done

## MAKE file for T,X profiles
for i in `ls Main*.dat`
do
	echo "Filtering TX $i"
#	awk "{ if(\$4==${valeurX}) print }" $i > $i.2DX
done

## MAKE file for T, Z profiles
for i in `ls Main*.dat`
do 
	echo "Filtering TZ $i"
	awk '{if($3==0) print }' $i > $i.2DZ
done

############ BUILD zOfFluence.dat ##############
# build an unordered Fluences.dat file
set LC_ALL=C
rm Fluences.tmp
rm zOfFluence*.tmp 

for i in `ls Main*.dat -1`
do
	awk '{print $2}' $i | head -n1 >> Fluences.tmp
done



# extract unordrer datas from TIME*.dat
# rm Time0.dat
for i in `ls Time*.dat -1`
do
	awk -f '/home/thibault/Documents/LaAPT/Scripts/MeltingOfFluence-Enhanced.awk' $i >> zOfFluence.tmp
done
sort -g Fluences.tmp > Fluences.dat
awk '{ printf("%15.15f\t%15.15f\t%15.15f\t%15.15f\t%15.15f\t%g\t%g\t%g\t%g\n", $1, $2, $3, $4, $5, $6, $7, $8, $9)}' zOfFluence.tmp > zOfFluenceSorted.tmp
# sort -g zOfFluence.tmp > zOfFluenceSorted.tmp
paste Fluences.dat zOfFluenceSorted.tmp > zOfFluence.dat

rm Fluences.tmp
rm zOfFluence.tmp
rm zOfFluenceSorted.tmp

####### Build Video ###########

rm *.vid
# awk -f '/home/thibault/Documents/LaAPT/Scripts/BuildFilesXZ.awk' DepthVessel.dat
awk -f '/home/thibault/Documents/LaAPT/Scripts/BuildMaps.awk' Depth.dat
# awk -f '/home/thibault/Documents/LaAPT/Scripts/BuildFilesXZ.awk' DualDepth.dat
ls *.vid

file=`ls -1 -t *.vid | sort -n | tail -n1`

# filter for Depth.dat
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$21}' ${file}  > SourceE.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$4}' ${file}  > Intensity.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$5}' ${file}  > Te.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$6}' ${file}  > Th.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$7}' ${file}  > Ts.map
# awk '{ if($2<100e-6) print $2*1E6,$3*1E6,1E-9*($35**2+$36**2)**0.5 }' $file  > Field.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,($44**2)**0.5 }' $file  > Field.map
#  awk '{ if($2==0 && $3==0) print $1,$35,$36 }' Depth.dat > TimeApexField.dat
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$8}' ${file}  > Ne.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$9}' ${file}  > Nh.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$32}' ${file} > HeatingDelayMap.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,$33}' ${file} > HeatingMaxTemp.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,1E-9*($35) }' $file  > FieldX.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,1e-9*($36) }' $file  > FieldY.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,1e-9*($40) }' $file  > GradNeX.map
 awk '{ if($2<100e-6) print $2*1E6,$3*1E6,1e-9*($41) }' $file  > GradNeY.map

# filter for DepthVessel.dat
#awk '{ if($2<2e-6) print $2*1E6,$3*1E6,$4 }' $file  > Potential.map
#awk '{ if($2<2e-6) print $2*1E6,$3*1E6,1E-9*($5**2+$6**2)**0.5 }' $file  > Field.map
#awk '{ if($2<2e-6) print $2*1E6,$3*1E6,1E-9*($5) }' $file  > FieldX.map
#awk '{ if($2<2e-6) print $2*1E6,$3*1E6,1e-9*($6) }' $file  > FieldY.map

rm *.vid.eps
echo "`ls *.vid | wc -l` pictures in Video"
for i in `ls -tr -1 *.vid | sort -n`
do
	/home/thibault/Documents/LaAPT/Scripts/plotVideo.sh $i
done

/home/thibault/Documents/LaAPT/Scripts/epsToGif-TTM.sh 
