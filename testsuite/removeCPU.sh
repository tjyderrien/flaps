#!/bin/bash

if [ -z $1 ]; then
	echo "Description: Reshape the old output data from the input file <Time*.dat>."
	echo "Usage: ./removeCPU.sh <TimeMax.dat | TimeApex.dat | TimeUp.dat | TimeDown.dat> > NewFile.dat"
else 
	if [ "$1" == "TimeMax.dat" ] || [ "$1" == "TimeApex.dat" ] || [ "$1" == "TimeUp.dat" ] || [ "$1" == "TimeBottom.dat" ]; then
		awk '{ $10=""; $19=""; $20=""; $36=""; $37=""; print }' $1 
	else 
		echo "Please use a file name Time<Max | Apex | Up | Down>.dat."
	fi
fi
