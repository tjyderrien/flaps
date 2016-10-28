#! /bin/bash

rm *.gif

chain=""
for i in `ls -1 -tr 2D-Intensity*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Intensity."
convert -density 120 `cat file.tmp` 2D-Intensity-video.gif
rm file.tmp
############################

chain=""
for i in `ls -1 -tr 2D-ZoomedIntensity*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed Intensity."
convert -density 120 `cat file.tmp` 2D-ZoomedIntensity-video.gif
rm file.tmp
############################


chain=""
for i in `ls -1 -tr 2D-logIntensity*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... logIntensity."
convert -density 120 `cat file.tmp` 2D-logIntensity-video.gif
rm file.tmp
###########################

chain=""
for i in `ls -1 -tr 2D-ZoomedlogIntensity*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed logIntensity."
convert -density 120 `cat file.tmp` 2D-ZoomedlogIntensity-video.gif
rm file.tmp
###########################


chain=""
for i in `ls -1 -tr 2D-logNe*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... logNe."
convert -density 120 `cat file.tmp` 2D-logNe-video.gif
rm file.tmp
############################

chain=""
for i in `ls -1 -tr 2D-ZoomedLogNe*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... ZoomedLogNe."
convert -density 120 `cat file.tmp` 2D-ZoomedLogNe-video.gif
rm file.tmp
############################

chain=""
for i in `ls -1 -tr 2D-ZoomedNe*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed Ne."
convert -density 120 `cat file.tmp` 2D-Zoomed-Ne-video.gif
rm file.tmp
############################


chain=""
for i in `ls -1 -tr 2D-Ne*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Ne."
convert -density 120 `cat file.tmp` 2D-Ne-video.gif
rm file.tmp
############################
chain=""
for i in `ls -1 -tr 2D-Te*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Te."
convert -density 120 `cat file.tmp` 2D-Te-video.gif
rm file.tmp

############################
chain=""
for i in `ls -1 -tr 2D-ZoomedTe*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed Te."
convert -density 120 `cat file.tmp` 2D-ZoomedTe-video.gif
rm file.tmp


############################
chain=""
for i in `ls -1 -tr 2D-Ts*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Ts."
convert -density 120 `cat file.tmp` 2D-Ts-video.gif
rm file.tmp

############################
chain=""
for i in `ls -1 -tr 2D-ZoomedTs*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed Ts."
convert -density 120 `cat file.tmp` 2D-ZoomedTs-video.gif
rm file.tmp


############################
chain=""
for i in `ls -1 -tr 2D-Ex*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Ex."
convert -density 120  `cat file.tmp` 2D-Ex-video.gif
rm file.tmp

############################
chain=""
for i in `ls -1 -tr 2D-ZoomedEx*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed Ex."
convert -density 120  `cat file.tmp` 2D-ZoomedEx-video.gif
rm file.tmp


############################
chain=""
for i in `ls -1 -tr 2D-logEx*.eps`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... logEx."
convert -density 120 `cat file.tmp` 2D-logEx-video.gif
rm file.tmp

############################
chain=""
for i in `ls -1 -tr 2D-ZoomedLogEx*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... Zoomed logEx."
convert -density 120 `cat file.tmp` 2D-ZoomedlogEx-video.gif
rm file.tmp


############################
chain=""
for i in `ls -1 -tr 2D-logEnthalpy*.eps`
do
        chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] ... logEnthalpy."
convert -density 120 `cat file.tmp` 2D-logEnthalpy-video.gif
rm file.tmp

