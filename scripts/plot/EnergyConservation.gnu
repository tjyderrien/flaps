#!gnuplot

reset

set output '20150509-EnergyPotvsKin.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24' 

# set key out right horizontal maxrows 4 maxcols 2
# set terminal dumb'

set format "%g"
set xrange [1e-3:10]
xscale=1E12
set log x
unset log y
set key left

set xlabel 'Time (ps)'
set ylabel 'Energy ({/Symbol m}J/m)'

# xscale=1e12
# yscale=1e6

plot "EnergyConservation.old" u 1:7 w l lc 1 lw 5 t 'Incoming laser energy'


# , \
# "EnergyConservation.old" u 1:2 w l lc 1 lw 4 t 'Abs. Laser energy', \
# "EnergyBalance.old" u 1:2 w l lc 3 lw 3 t 'Kin, e', \
# "EnergyBalance.old" u 1:5 w l lc 3 lw 4 t 'Pot. e', \
# "EnergyConservation.old" u 1:3 w l lc 3 lw 2 t 'Kin.+Pot. e-', \
# "EnergyConservation.old" u 1:4 w l lc 4 lw 2 t 'Kin.+Pot. h', \
# "EnergyConservation.old" u 1:5 w l lc 7 lw 3 t 'Lattice'


#reset
## set terminal dumb
#set terminal postscript eps enhanced color font 'Helvetica, 26'
#set output '20150108-EnergyConservation.eps'
#
#set format "%g"
#
#xscale=1e12
#yscale=1e6
#
#set log x
#unset log y
#set key left
#
#set xlabel 'Time (ps)'
#set ylabel 'Energy (uJ/m)'

#plot "EnergyConservation.dat" u ($1*xscale):($7*yscale) w l lc 1 lw 5 t 'Incoming laser energy', \
#"EnergyConservation.dat" u ($1*xscale):($2*yscale) w l lc 1 lw 4 t 'Abs. Laser energy', \
#"EnergyConservation.dat" u ($1*xscale):($8*yscale) w l lc 3 lw 3 t 'Kin, e', \
#"EnergyConservation.dat" u ($1*xscale):($9*yscale) w l lc 3 lw 4 t 'Pot. e', \
#"EnergyConservation.dat" u ($1*xscale):($3*yscale) w l lc 3 lw 2 t 'Kin.+Pot. e-', \
#"EnergyConservation.dat" u ($1*xscale):($4*yscale) w l lc 4 lw 2 t 'Kin.+Pot. h', \
#"EnergyConservation.dat" u ($1*xscale):($5*yscale) w l lc 7 lw 5 t 'Lattice'

############## Reconstructed energy
# reset
# After reconstruction using RebuildEnergy.sh, we plot kinetic energy and potential energy.

# set output
# set terminal x11

set output '20150509-EnergyPotvsKin.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24' 
# set key out right horizontal maxrows 2
set key out right vertical maxcolumn 1

plot "EnergyConservation.dat" u ($1*xscale):($7*yscale) w l lc 1 lt 2 lw 1 t 'Laser energy', \
"EnergyConservation.dat" u ($1*xscale):($2*yscale) w l lc 1 lw 5 lt 1 t 'Abs. energy', \
"EnergyBalance.dat" u ($1*xscale):($2*yscale) w l lc 3 lw 3 lt 2 t 'Kin, e', \
"EnergyBalance.dat" u ($1*xscale):($5*yscale) w l lc 3 lw 3 lt 3 t 'Pot. e', \
"EnergyConservation.dat" u ($1*xscale):($3*yscale) w l lc 3 lw 5 lt 1 t 'Kin.+Pot. e-', \
"EnergyConservation.dat" u ($1*xscale):($5*yscale) w l lc 7 lw 5 lt 1 t 'Lattice'

# "EnergyConservation.dat" u ($1*xscale):($4*yscale) w l lc 4 lw 2  t 'Kin.+Pot. h', \
