##################
#CONVERT TO LIGHT FILES
##################
for i in `ls -1 $1`
do
	gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=$i.png -dBATCH -dNOPAUSE $i
done
