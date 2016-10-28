#!gnuplot

reset

set terminal dumb

# set terminal postscript eps enhanced color font 'Helvetica, 26'
# set output '20150108-EnergyConservation.eps'

set format "%g"

xscale=1E12
set log x
unset log y
set key left

set xlabel 'Time (ps)'
set ylabel 'Energy (J)'

plot "EnergyConservation.dat" u ($1*xscale):2 w l t 'Laser', "EnergyConservation.dat" u ($1*xscale):3 w l t 'Electron', "EnergyConservation.dat" u ($1*xscale):4 w l t 'Hole', "EnergyConservation.dat" u ($1*xscale):5 w l t 'Lattice', "EnergyConservation.dat" u ($1*xscale):7 w l t 'Laser intensity energy'
