#!/bin/bash

# We will run several values of rotation angles and start one parallel calculation per value

# 1. Define all the desired values 
# 2. Replace the value into the Mie_param.dat file
# 3. Modify the corresponding submit.pbs
# 4. Submit the corresponding run

# values=(0 45 90 180 270)
values=`seq 0 45 270`

rm input_Mie.txt submit.pbs
rm submit*.pbs
for i in ${values[*]}
do
	echo "Angle = $i deg"
	# rm input_Mie.txt submit.pbs
	# sed "s/MieValue/$i/g" input_Mie.sav > input_Mie.txt #replace in input file
	sed "s/MieValue/$i/g" submit.sav > submit$i.pbs #prepare batch submission
	# git add input_Mie.txt submit.pbs
	# git commit -m "IT4I: changed the Mie orientation plane to $i"
	qsub submit$i.pbs && echo "Job $i submitted"
	# sleep 1 #let time to qsub to actually submit !
done
