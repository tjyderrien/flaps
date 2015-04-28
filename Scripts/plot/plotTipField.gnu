#!gnuplot

reset

set xlabel 'X [{/Symbol m}m]' 
set ylabel 'Y [{/Symbol m}m]' 

xscale=1E6
yscale=1E6

set view map

splot "Field.dat" u ($1*xscale):($2*xscale):3 w pm3d


