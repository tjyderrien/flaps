#!/bin/bash

for i in `du | awk '{ print $2 }'`
do
	echo $i
	tail -n1 $i/TimeMax.dat | awk '{ print $10 }' 
done
