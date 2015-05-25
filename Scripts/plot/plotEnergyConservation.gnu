#!gnuplot

reset

# set output
# set terminal x11

set xlabel 'Time (s)'
set ylabel 'Energy (J.m^{-1})'

set output '20150503-EnergyConservation.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

set key outside center bottom
set grid
unset log x
unset log y
unset log y2
set xtics format "10^{%L}"
set ytics nomirror
set y2tics nomirror

xscale=1E0

plot "./TimeMax.dat" u ($1*xscale):32 w l t 'e-', \
"./TimeMax.dat" u ($1*xscale):33 w l t 'h', \
"./TimeMax.dat" u ($1*xscale):($32+$33) w l lw 3 t 'h+e-', \
"./TimeMax.dat" u ($1*xscale):34 w l t 'Lattice', \
"./TimeMax.dat" u ($1*xscale):($32+$33+$34) w l lc 8 lw 5 t 'Total in the solid', \
"./TimeMax.dat" u ($1*xscale):($30) w l lc 7 lw 3 t 'Absorbed optical energy'
# "./TimeMax.dat" u ($1*xscale):($35) w l lc 7 lw 1 t 'Laser optical energy' axis x1y2

#, \

# "./TimeMax.dat" u 1:($32+$33+$34) w l lw 5 t 'Total solid'

