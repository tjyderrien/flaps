#!/bin/bash

Simplify() {
	awk '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9 }' $1
}


if [ -z $2 ]; then
	echo "Description: Compare the 9 first columns of the <file1> and <file2>."
	echo "Required: this script needs numdiff installed. See http://www.nongnu.org/numdiff/." 
	echo "Usage: ./compareOutputs.sh <File1> <File2>"
else
	echo "Info: Files were saved in BAK, in case."
	cp $1 $1.bak
	cp $2 $2.bak
	echo "Info: Simplifying..."
	Simplify $1 > $1.simplified.dat
	Simplify $2 > $2.simplified.dat
	echo "Info: Performing the comparison..."
	numdiff $1.simplified.dat $2.simplified.dat
fi
