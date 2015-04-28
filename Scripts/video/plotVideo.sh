#!/bin/sh
# module load gnuplot
gnuplot << EOF
set terminal postscript eps enhanced color font 'Helvetica, 22'
set output "$1.eps"
set view equal xy
set view map 
set xlabel 'X [{/Symbol m}m]'
set ylabel 'Y [{/Symbol m}m]'
# set cbtics format "%1.0l x 10^{%L}"

# set palette defined (0 "blue", 1 "white", 2 "red")
set palette rgb 33,13,10 #rainbow
# set palette rgb 30,31,32 #many colors

set size ratio -1.0
set xrange [0:5]
set yrange [:]
# set cbrange [0.9999e5:1.0001e5]
# set cbrange [0.99e25:1.01e25]
# set cbrange [0:100]
# set cbrange [80:1E5]
# set cbrange [1E0:1E27]
# set cbrange [1E0:1E16]
unset logscale cb
# set cbrange [1e-1:]
set key out left
# set format "%g"
# set pm3d map interpolate 4,4

set xtics 5
set ytics 2
# set cbtics 
# set key outside left 
# splot "< awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk' $1" u (\$2*1E6):(\$3*1E6):((\$5**2+\$6**2)**0.5) w pm3d t columnheader 1
# splot "< awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk' $1" u (\$2*1E6):(\$3*1E6):((\$16**2+\$17**2)**.5) w pm3d t columnheader 1
# splot "< awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk' $1" u (\$2*1E6):(\$3*1E6):(\$7) w pm3d t columnheader 1
splot "< awk -f 'Scripts/video/plotXZ.awk' $1" u (\$1):(\$2):((\$3)) w pm3d t columnheader 1

EOF
