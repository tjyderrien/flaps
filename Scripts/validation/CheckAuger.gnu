#!gnuplot

###### FREE CARRIER DENSITY

# Validation of the Auger recombination

reset

set output '20140314-ValidationOfAuger.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 26'

t0=-1.9e-13-1e-14
Ne0=1E27
Nh0=1E27
CarE=2.3e-43
CarH=7.8e-44

ec=1.6e-19
Egap=1.16*ec
kb=1.38e-23

set log x

set xtics format "10^{%L}"

set xlabel 'Time (s)'
set ylabel 'Error on the free-carrier density (%)'

set key left
set grid

xscale=1E0
yscale=1E0

NAuger(t)=1e0/(2e0*t*Ne0**2*CarE+2e0*t*Ne0**2*CarH+1e0)**(1e0/2e0)*Ne0
# CarH*Ne0*Nh0/((Ne0*CarE+Nh0*CarH) * exp(CarH*Nh0**2*t) - CarE*Ne0)

plot "CheckAuger/TimeMax-SimpleNe.dat" u (($1-t0)*xscale):(100e0*($5-NAuger($1-t0))/NAuger($1-t0)*yscale) w p t 'Error on N_{e,h} (%)'

# plot NAuger(x) w l t 'Analytic', \
# "../TimeMax.dat" u ($1-t0):($5) w p t 'Numeric' 

####### Coupled density and temperature
set output '20140314-ValidationOfAugerCoupledTeNe.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

Te0=80e0
Th0=80e0
Ne0=1E27
Nh0=1E27

# CarE=2.3e-43 #in the case that Ne=Nh

NAuger(t)=(1e0/(Ne0**2) + 2e0*CarE*t + 2e0*CarH*t)**(-0.5e0)

# Non-splitted solutions
TeAuger(t)=(-(2e0/3e0)*Egap*CarE/(sqrt(2e0*t*CarE+2e0*t*CarH+1e0/Ne0**2)*(CarE+CarH)*kb)+(1e0/3e0)*(2e0*(1e0/Ne0**2)**((1e0/2e0)*CarE/(CarE+CarH))*(1e0/Ne0**2)**((1e0/2e0)*CarH/(CarE+CarH))*Egap*CarE*Ne0+3e0*Te0*kb*CarE+3e0*Te0*kb*CarH)/((1e0/Ne0**2)**((1e0/2e0)*CarE/(CarE+CarH))*(1e0/Ne0**2)**((1e0/2e0)*CarH/(CarE+CarH))*kb*(CarE+CarH)))*exp(log((2e0*CarE+2*CarH)*t+1e0/Ne0**2)*CarE/(2e0*CarE+2e0*CarH)+log((2e0*CarE+2e0*CarH)*t+1e0/Ne0**2)*CarH/(2e0*CarE+2e0*CarH))

ThAuger(t)=(-(2e0/3e0)*Egap*CarH/(sqrt(2e0*t*CarE+2e0*t*CarH+1e0/Ne0**2)*(CarE+CarH)*kb)+(1e0/3e0)*(2e0*(1e0/Ne0**2)**((1e0/2e0)*CarE/(CarE+CarH))*(1e0/Ne0**2)**((1e0/2e0)*CarH/(CarE+CarH))*Egap*CarH*Ne0+3e0*Th0*kb*CarE+3e0*Th0*kb*CarH)/((1e0/Ne0**2)**((1e0/2e0)*CarE/(CarE+CarH))*(1e0/Ne0**2)**((1e0/2e0)*CarH/(CarE+CarH))*kb*(CarE+CarH)))*exp(log((2e0*CarE+2e0*CarH)*t+1e0/Ne0**2)*CarE/(2e0*CarE+2e0*CarH)+log((2e0*CarE+2e0*CarH)*t+1e0/Ne0**2)*CarH/(2e0*CarE+2e0*CarH))

# splitted solutions (neglecting the couplnig)
# TeAuger(t)=(2e0/3e0)*Egap*CarE*Ne0*Nh0*t/kb+Te0
# ThAuger(t)=(2e0/3e0)*Egap*CarH*Nh0*Ne0*t/kb+Th0

set ylabel 'Error on Density (%)'
set y2label 'Error on Temperature (%)'
set ytics nomirror
set y2tics nomirror

set xtics format "10^{%L}"
set key out bottom center horizontal

y2scale=1e0

plot "CheckAuger/TimeMax-CoupledCeNew.dat" u (($1-t0)*xscale):(100e0*($2-TeAuger($1-t0))/TeAuger($1-t0)*y2scale) w p lw 1 lc 1 t 'Error on T_e (%)' axis x1y2, \
"CheckAuger/TimeMax-CoupledCeNew.dat" u (($1-t0)*xscale):(100e0*($3-ThAuger($1-t0))/ThAuger($1-t0)*y2scale) w p lw 1 lc 3 t 'Error on T_h (%)' axis x1y2, \
"CheckAuger/TimeMax-CoupledCeNew.dat" u (($1-t0)*xscale):(100e0*($5-NAuger($1-t0))/NAuger($1-t0)*yscale) w p lw 1 lc 7 t 'Error on N_e (%)' axis x1y1, \
"CheckAuger/TimeMax-CoupledCeNew.dat" u (($1-t0)*xscale):(100E0*($6-NAuger($1-t0))/NAuger($1-t0)*yscale) w p lw 1 lc 7 t 'Error on N_h (%)' axis x1y1





# TeAuger(x) w l lc 1 lw 5 t 'T_e, analytic' axis x1y2, \
# ThAuger(x) w l lc 3 lw 5 t 'T_h, analytic' axis x1y2, \
# NAuger(x) w l lc 7 lw 5 t 'N_{e,h}, analytic' axis x1y1, \
