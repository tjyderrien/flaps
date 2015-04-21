#! /bin/bash

valeurX=2.40024e-09
valeurY=0
valeurT=50e-15

####### Build Video ###########

rm *.vid
# awk -f 'Scripts/video/BuildMaps.awk' DepthVessel.dat'
awk -f 'Scripts/video/BuildMaps.awk' Depth.dat
# awk -f 'Scripts/video/BuildFilesXZ.awk' DualDepth.dat
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
