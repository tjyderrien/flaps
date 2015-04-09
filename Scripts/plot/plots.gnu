#! gnuplot
# Le fichier Plasmon.gnu contient les scripts generant des figures par rapport
# au temps et espace. 
# Ce script contient de quoi faire des analyses par rapport a la fluence

#### DIRECT FIGURE OF INTENSITY PERIOD ####

# system "head -n

reset
set output '20101104-SEWperiodIntensity.eps'
set terminal postscript eps enhanced color
set autoscale xy
set key right bottom outside
set title 'Period of the interfered intensity produced by SEW as a function of time'
set xlabel 'Time [ps]'
set ylabel 'Smallest SEW period [nm]'
plot 	"< awk '{ if($2 > -100e-15 && $2<100e-15) print }' Ks6e9/Plasmon0D.dat" using ($2*1e12):($5*1e9) w d title 'Ks=6e9', \
	"< awk '{ if($2 > -100e-15 && $2<100e-15) print }' Ks8e9/Plasmon0D.dat" using ($2*1e12):($5*1e9) w d title 'Ks=8e9', \
	800 w l title 'Laser wavelength'
set output
set terminal x11 1
replot 

### FIGURE OF LATTICE TEMPERATURE PERIOD with FLUENCE ###

# system "tail -n200000 Main*.dat | awk '{ if($3<10e-6 && $3>-10e-6 && $4==2.50017e-09) print }' > Main2DX.tmp"

set terminal x11 2
set output

set autoscale xy
f(x)=A*cos(2*pi/lambda*x)+B
A=30
lambda=800e-9
B=1280
n=1

plot "< awk '{ if($1==11625 && $3<3e-6 && $3>-3e-6) print }' Main2DX.tmp" using 3:6 w lp, f(x) title "Fit"
fit f(x) "< awk '{ if($1==11625 && $3<3e-6 && $3>-3e-6)  print }' Main2DX.tmp" using 3:6 via A,B,lambda
set label 1 "{/Symbol L}=%g", lambda at 0,B-1.5*A
set label 2 "n=%g", n at 0, B+1.5*A
replot
