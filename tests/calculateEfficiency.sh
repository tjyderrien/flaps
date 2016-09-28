#!/bin/bash

python_calc () {
    python -c "print ($@)"
}

calculationDuration () { #valid
	BeginningSecond=$( date -u -d "$1 $2" +"%s.%N" )
	EndingSecond=$( date -u -d "$3 $4" +"%s.%N" )
	# echo "$1 $2: $BeginningSecond seconds"
	# echo "$3 $4: $EndingSecond seconds"
	echo "$EndingSecond - $BeginningSecond" | bc 
}

PhysicalTime () { #valid
	Starting=$( head -n1 TimeMax.dat | awk '{ print $1 }' )
	Ending=$( tail -n1 TimeMax.dat | awk '{ print $1 }' )
	echo $( python_calc "$Ending - $Starting" )
}

for i in $( du --max-depth=1 | awk '{ print $2 }' | head -n-1 )
do
	cd $i
	echo $i
	beginning=$( ls parameters.dat --full-time | awk '{ print $6, $7 }' | tr '-' '/' )
	ending=$( ls -lhtr --full-time | tail -n1 | awk '{ print $6, $7 }' | tr '-' '/' )
#	echo "$beginning, $ending"
	HumanDuration=$( calculationDuration $beginning $ending )
	echo "Human duration: $HumanDuration"
	ProgramDuration=$( PhysicalTime )
	echo "Physical duration: $ProgramDuration"
	
	Efficiency=$( python_calc " $ProgramDuration * 1E12 / $HumanDuration * ( 60 * 60 * 24 ) " ) #physical s per human second
	echo "Efficiency: $Efficiency ps / day."
	
	echo ""
	
	cd - > /dev/null	

done
