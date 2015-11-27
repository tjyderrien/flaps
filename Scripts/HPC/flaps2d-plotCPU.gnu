#!gnuplot

reset

set format "%g"

set xlabel 'Time (ps)'
set ylabel 'CPU'

xscale=1E12

set key top right
set terminal dumb

unset log y

plot "337824.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '1 PBS, 1 OMP', \
"337853.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '1 PBS, 24 OMP', \
"337908.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '1 PBS, 24 OMP, new'

plot "337857.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '24 PBS, 1 OMP', \
"337866.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '24 PBS, 24 OMP'

# "337857.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '24 PBS, 1 OMP', \

# "337866.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '24 PBS, 24 OMP'

#, \
# "337391.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '1 PBS, 24 OMP', \
# "337621.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '24 PBS, 1 OMP', \
# "337490.isrv5/flaps2d/TimeMax.dat" u ($1*1E12):10 w l t '24 PBS, 24 OMP'
