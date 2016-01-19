#!gnuplot
reset

set output '20140317-CheckCoupling.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

set xtics format "10^{%L}"

tau(ne)=240e-15*(1+(ne/nth)**2e0)
Te(t)=Tl+exp(-t/tau(ne))*(Te0-Tl)

Tl =80e0
ne = 1e+27
Te0 = 3e5
nth = 6.02e+26

shift = -2e-13

set log x

set xlabel 'Time (s)'
set ylabel 'Relative error on temperature (%)'
plot "CheckCoupling/TimeMax-Te3e5K-Ne1e27const.dat" u ($1-shift):(100e0*(($2-Te($1-shift)))/(Te($1-shift))) w l t 'Error (%)' 

# , \
# Te(x) w l t 'Analytic'
