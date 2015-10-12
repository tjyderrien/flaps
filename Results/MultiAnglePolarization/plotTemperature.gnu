#!gnuplot

reset

#WhatToPlot=($2)

set output 'AngleOfOrientation.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set key bottom right

set xlabel 'Time (ps)'
set ylabel 'Lattice temperature (K)'

xscale=1E12

plot "370850.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 0', \
"370851.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 45', \
"370852.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 90', \
"370853.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 135', \
"370854.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 180', \
"370855.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 225', \
"370856.isrv5/TimeApex.dat" u ($1*xscale):4 w l t 'TE, 270'
