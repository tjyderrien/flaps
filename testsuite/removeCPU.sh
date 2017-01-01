#!/bin/bash

if [ -z $1 ]; then
	echo "Description: Removes the column CPU_usage from the input file <Time*.dat>."
	echo "Usage: ./removeCPU.sh <TimeMax.dat | TimeApex.dat | TimeUp.dat | TimeDown.dat>"
else 
	if [ $1== ]; then
		awk '{ $10=""; print }' $1 
	else 
		echo "Please use a file name Time<Max | Apex | Up | Down>.dat."
	fi
fi
