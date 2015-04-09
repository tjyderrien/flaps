#!gnuplot

reset

set output
set terminal x11

set title 'Flux evaporated'
set xlabel 'Temperature increase'
set ylabel 'Flux'


kb=1.38e-23
ec=1.6e-19
Q=0.15*ec
N=1
nu=1e3
phi(T)=N*nu*exp(-Q/(kb*T))

set xrange [1:300]
# set yrange [1e-3:1]
set log y
plot phi(x) w l t '{/Symbol p}(T)'

set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
set output '20130707-Arrhenius.eps'
replot
