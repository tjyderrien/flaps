#!gnuplot

FlapsRoot="/home/thibault/Documents/Codes/Flaps/20150428-flaps-gmsh/flaps2d"

reset

set output 'Field.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22' size 8cm, 5cm

# set view equal xy
set view map

set palette rgb 33,13,10 #rainbow

set xlabel 'X [{/Symbol m}m]' 
set ylabel 'Y [{/Symbol m}m]' 

set xtics 1.0
set ytics 0.1

xscale=1E6
yscale=1E6

unset logscale cb
set key out left


splot "< awk -f '".FlapsRoot."'/Scripts/plot/plotField.awk Field.dat" u ($1*xscale):($2*xscale):3 w pm3d notitle


