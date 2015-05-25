#! /bin/bash

################# Video for any value ###########
Physical="Ts-needle"

################ SCRIPT to not modify #############
rm file.tmp
rm "2D-$Physical-video.gif"
chain=""
for i in `ls -1 -t *.vid.eps | sort -n`
do
	chain+=" $i"
done
echo ${chain} | tr "\n" " " > file.tmp

echo "[Building video] $Physical..."
convert -verbose -delay 50 -alpha Off -density 240 `cat file.tmp` 2D-$Physical-video.gif
rm file.tmp
