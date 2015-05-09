#!/bin/bash

#CONVERT EPS TO PNG FILES with good quality

for i in `ls -1 $1`
do
	gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=$i.png -dBATCH -dNOPAUSE $i
done
