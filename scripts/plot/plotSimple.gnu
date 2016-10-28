#!gnuplot

reset

xscale=1E12

set key left

set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'
plot "TimeApex.dat" u ($1*xscale):4 w l t 'Apex'
