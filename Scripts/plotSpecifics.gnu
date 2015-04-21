#!gnuplot

##### Lattice heating with laser wavelength
reset
unset log y
set key below vertical
set grid
set format "%g"
set xlabel 'Time (ns)'
set ylabel 'Temperature [K]'
set xrange [:10]
set yrange [:]
set output '20120806-LatticeHeatingWithWavelength.eps'
set terminal postscript eps enhanced monochrome
scale=1E9
plot "04-343nm-Source/Cines/CrossDiffusion/500fs/10J/TimeApex.dat" u ($1*scale):4 w l lw 2 t '500 fs, 343 nm', \
"04-515nm-Source/Cines/CrossDiffusion/500fs/dt10as/CooledBnds/10J/TimeApex.dat" u ($1*scale):4 w l lw 2 t '500 fs, 515 nm'
# "04-1030nm-Source/Direct/500fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '500 fs, 1030 nm'
# "04-343nm-Source/Cines/CrossDiffusion/40fs/10J/TimeApex.dat" u ($1*scale):4 w l lw 2 t '500 fs, 343 nm', \
# "04-1030nm-Source/Direct/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 1030 nm', \

# "04-1030nm-Source/Direct/DrudeModel2/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 1030 nm, Drude 2'

# "04-1030nm-Source/Direct/100fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '100 fs, 1030 nm', \

# "04-515nm-Source/ConstantDistribution/100fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '100 fs, 515 nm', \
# "04-515nm-Source/ConstantDistribution/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 515 nm', \



