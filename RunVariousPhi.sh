#!/bin/bash

# We will run several values of rotation angles and start one parallel calculation per value

# 1. Define all the desired values 
# 2. Replace the value into the Mie_param.dat file
# 3. Modify the corresponding submit.pbs
# 4. Submit the corresponding run

values=(0 90 180 270)

counter=0
for i in ${values[*]}
do
	
	sed "s/MieValue/$i/g" input_Mie.sav > input_Mie.txt #replace in input file
	sed "s/MieValue/$i/g" submit.sav > submit.pbs #prepare batch submission
	git add input_Mie.txt submit.pbs
	git commit -m "IT4I: changed the Mie orientation plane to $i"
	qsub submit.pbs
done
