#!gnuplot

reset

set output "ComplexityTest.eps"
set terminal postscript eps enhanced color font "Helvetica, 26"

set format "%g"

LocalPath="/home/thibault/Documents/Codes/Flaps/flaps2d/Results/ComplexityTest"
File="TimeMax.dat"

plot LocalPath."/M-N-11/".File u 1:10 w l t 'M=N=11', \
LocalPath."/M-N-21/".File u 1:10 w l t 'M=N=21', \
LocalPath."/M-N-51/".File u 1:10 w l t 'M=N=51', \
LocalPath."/M-N-101/".File u 1:10 w l t 'M=N=101', \
LocalPath."/M-N-151/".File u 1:10 w l t 'M=N=151'
