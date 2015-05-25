
######## Lattice heating with laser wavelength
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

######## Lattice heating with pulse duration
reset
set size 0.7
unset log y
set key out below vertical
set grid
set format "%g"
set xlabel 'Time [ps]'
set ylabel 'T_{L} [K]'
set xrange [:2]
set yrange [:]
set output '20120806-LatticeHeatingWithPulseDuration-515nm.eps'
set terminal postscript eps enhanced monochrome

plot "04-515nm-Source/ConstantDistribution/500fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '500 fs, 515 nm', \
"04-515nm-Source/ConstantDistribution/100fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '100 fs, 515 nm', \
"04-515nm-Source/ConstantDistribution/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 515 nm'

# "04-1030nm-Source/Direct/500fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '500 fs, 1030 nm'
# "04-1030nm-Source/Direct/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 1030 nm', \
# "04-343nm-Source/Direct/500fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '500 fs, 343 nm', \
# "04-343nm-Source/Direct/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 343 nm', \
# "04-1030nm-Source/Direct/DrudeModel2/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 1030 nm, Drude 2'

# "04-1030nm-Source/Direct/100fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '100 fs, 1030 nm', \

######## Lattice heating with pulse intensity

reset
set size 0.7
unset log y
set key out below vertical
set grid
set format "%g"
set xlabel 'Time [ps]'
set ylabel 'T_{L} [K]'
set xrange [:2]
set yrange [:]
set output '20120806-LatticeHeatingWithPulseConstantIntensity-515nm.eps'
set terminal postscript eps enhanced monochrome

plot "04-515nm-Source/ConstantDistribution/500fs/Fluences/0.2GW/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '0.2 GW.cm^{-2}', \
"04-515nm-Source/ConstantDistribution/500fs/Fluences/0.6GW/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '0.6 GW.cm^{-2}', \
"04-515nm-Source/ConstantDistribution/500fs/Fluences/1.0GW/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '1.0 GW.cm^{-2}', \
"04-515nm-Source/ConstantDistribution/500fs/Fluences/1.2GW/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '1.2 GW.cm^{-2}'


# "04-515nm-Source/ConstantDistribution/500fs/Fluences/1.2GW/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '1.2 GW.cm^{-2}'

# "04-1030nm-Source/Direct/500fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '500 fs, 1030 nm'
# "04-1030nm-Source/Direct/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 1030 nm', \
# "04-343nm-Source/Direct/500fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '500 fs, 343 nm', \
# "04-343nm-Source/Direct/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 343 nm', \
# "04-1030nm-Source/Direct/DrudeModel2/40fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '40 fs, 1030 nm, Drude 2'

# "04-1030nm-Source/Direct/100fs/TimeMax.dat" u ($1*1E12):4 w l lw 2 t '100 fs, 1030 nm', \

######## Lattice heating with pulse intensity
unset multiplot
reset
set output 
set terminal x11 enhanced
# set size 0.7
unset log y
set key out below vertical
set grid
set format "%g"
set xlabel 'Time [ps]'
set ylabel 'T_{L} [K]'
set xrange [:2]
set yrange [:]
set output '201208013-LatticeHeatingWithPulseFDTDIntensity.eps'
set terminal postscript eps enhanced monochrome
set multiplot layout 1,3
plot "04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/1Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2 lt 1 lc 1 t '1 J.m^{-2}, 1030 nm, max', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/1Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2 lt 2 lc 1 t '1 J.m^{-2}, 1030 nm, apex', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/5Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2 lt 1 lc 2 t '5 J.m^{-2}, 1030 nm, max', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/5Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2 lt 2 lc 2 t '5 J.m^{-2}, 1030 nm, apex', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/10Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2 lt 1 lc 3 t '10 J.m^{-2}, 1030 nm, max', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/10Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2 lt 2 lc 3 t '10 J.m^{-2}, 1030 nm, apex', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/50Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2 lt 1 lc 4 t '50 J.m^{-2}, 1030 nm, max', \
"04-1030nm-Source/NonLinearAbsorption/ICPEPA-40fs-1030nm-InfraredDistrib/ChenConductivityNoDrift/50Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2 lt 2 lc 4 t '50 J.m^{-2}, 1030 nm, apex'

plot "04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/1Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2 lt 1 lc 1 t '1 J.m^{-2}, 515 nm, max', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/1Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 1 t '1 J.m^{-2}, 515 nm, apex', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/5Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 1 lc 2 t '5 J.m^{-2}, 515 nm, max', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/5Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 2 t '5 J.m^{-2}, 515 nm, apex', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/10Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 1 lc 3 t '10 J.m^{-2}, 515 nm, max', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/10Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 3 t '10 J.m^{-2}, 515 nm, apex', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/50Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 1 lc 4 t '50 J.m^{-2}, 515 nm, max', \
"04-515nm-Source/GreenDistribution/ChenConductivity/40fs/3um/50Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 4 t '50 J.m^{-2}, 515 nm, apex'

plot "04-343nm-Source/ChenConductivityNoDrift/40fs/1Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2 lt 1 lc 1 t '1 J.m^{-2}, 343 nm, max', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/1Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 1 t '1 J.m^{-2}, 343 nm, apex', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/5Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 1 lc 2 t '5 J.m^{-2}, 343 nm, max', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/5Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 2 t '5 J.m^{-2}, 343 nm, apex', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/10Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 1 lc 3 t '10 J.m^{-2}, 343 nm, max', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/10Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 3 t '10 J.m^{-2}, 343 nm, apex', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/50Jm2/TimeMax.dat" u ($1*1E12):4 w l lw 2  lt 1 lc 4 t '50 J.m^{-2}, 343 nm, max', \
"04-343nm-Source/ChenConductivityNoDrift/40fs/50Jm2/TimeApex.dat" u ($1*1E12):4 w l lw 2  lt 2 lc 4 t '50 J.m^{-2}, 343 nm, apex'

unset multiplot 

####### Depth profile

# reset
# set size 0.7
# set output '20120725-1030nm-Xprofile-35fs.eps'
# set terminal postscript eps enhanced color
# 
# set format "%g"
# set xlabel 'Depth [um]'
# set ylabel 'T [K]'
# set y2label 'N [m^{-3}]'
# set ytics nomirror format "%g"
# set y2tics nomirror format "%3.1e"
# set xtics format "%g"
# set log y2
# set grid
# set key right
# set log x
# 
# plot "< awk '{ if((($3)^2)^.5 < 10e-9) print }' 1240.vid" u 2:5 w l t 'Te [K]', \
# 	"< awk '{ if((($3)^2)^.5 < 10e-9) print }' 1240.vid" u 2:6 w l t 'Th [K]', \
# 	"< awk '{ if((($3)^2)^.5 < 10e-9) print }' 1240.vid" u 2:7 w l t 'Ts [K]', \
# 	"< awk '{ if((($3)^2)^.5 < 10e-9) print }' 1240.vid" u 2:8 axis x1y2 w l t 'Ne [m^{-3}]', \
# 	"< awk '{ if((($3)^2)^.5 < 10e-9) print }' 1240.vid" u 2:9 axis x1y2 w l t 'Nh [m^{-3}]'
# 	
# 
# set terminal x11 enhanced 1
# set size 1.0
# replot

###################################"
# 
# set terminal x11 2
# 
# reset
# 
# set grid
# set log y
# set xtics format "%g"
# set ytics format "%g"
# set xlabel 'Time (ps)'
# set ylabel 'J.m^{-3}'
# plot "TimeMax.dat" u 1:8 w l t 'Laser energy', \
# "TimeMax.dat" u 1:9 w l t 'Thermal energy'
 
 
# set terminal x11 3
# reset
# set pm3d 
# set view map
# set surface
# set xtics format "%g"
# set ytics format "%g"
# set ztics format "%3.1e"
# set xlabel 'X (um)'
# set ylabel 'Y (um)'
# set zlabel 'Te'
# splot "< awk '{ if(fabs($1) < 30e-15) print }' Depth.dat" u ($2*1E6):($3*1E6):5 w pm3d

######### contribution of ionization mechanisms ######
reset

# set output 
# set terminal x11 enhanced

set output '20120828-ProcessContribution.eps'
set terminal postscript eps enhanced color

set ytics nomirror
set y2tics nomirror

set grid

set key out below horizontal

set format x "%g"
# set format y "%2.0t{/Symbol \327}10^{%L}"
set format y2 "%2.0t{/Symbol \327}10^{%L}"

set xlabel 'Time [ps]'
set ylabel 'Temperature [K]' textcolor rgbcolor "blue"
set y2label 'Density [m^{-3}]' textcolor rgbcolor "red"

set xrange [-0.2:0.2]

unset log y
unset log y2

plot "Normal/TimeMax.dat" u ($1*1E12):2 w l lc 3 lt 1 t '1ph + 2ph + Impact + IB' axis x1y1, \
"NoLinear/TimeMax.dat" u ($1*1E12):2 w l lc 3 lt 2 t '1ph Off' axis x1y1, \
"NoMPI/TimeMax.dat" u ($1*1E12):2 w l lc 3 lt 3 t '2ph Off' axis x1y1, \
"NoImpact/TimeMax.dat" u ($1*1E12):2 w l lc 3 lt 4 t 'Impact Off' axis x1y1, \
"NoIB/TimeMax.dat" u ($1*1E12):2 w l lc 3 lt 5 lw 3 t 'IB Off' axis x1y1, \
"NoDiffusion/TimeMax.dat" u ($1*1E12):2 w l lc 3 lt 6 t 'Diffusion Off' axis x1y1, \
"Normal/TimeMax.dat" u ($1*1E12):5 w l lc 1 lt 1 t '1ph + 2ph + Impact + IB' axis x1y2, \
"NoLinear/TimeMax.dat" u ($1*1E12):5 w l lc 1 lt 2 t '1ph Off' axis x1y2, \
"NoMPI/TimeMax.dat" u ($1*1E12):5 w l lc 1 lt 3 t '2ph Off' axis x1y2, \
"NoImpact/TimeMax.dat" u ($1*1E12):5 w l lc 1 lt 4 t 'Impact Off' axis x1y2, \
"NoIB/TimeMax.dat" u ($1*1E12):5 w l lc 1 lt 5 lw 3 t 'IB Off' axis x1y2, \
"NoDiffusion/TimeMax.dat" u ($1*1E12):5 w l lc 1 lt 6 t 'Diffusion Off' axis x1y2

reset

# set output 
# set terminal x11 enhanced

set output '20120828-ProcessContribution-Apex.eps'
set terminal postscript eps enhanced color

set ytics nomirror
set y2tics nomirror

set grid

set key out below horizontal

set format x "%g"
# set format y "%2.0t{/Symbol \327}10^{%L}"
set format y2 "%2.0t{/Symbol \327}10^{%L}"

set xlabel 'Time [ps]'
set ylabel 'Temperature [K]' textcolor rgbcolor "blue"
set y2label 'Density [m^{-3}]' textcolor rgbcolor "red"

set xrange [-0.2:0.2]

unset log y
unset log y2

plot "Normal/TimeApex.dat" u ($1*1E12):2 w l lc 3 lt 1 t '1ph + 2ph + Impact + IB' axis x1y1, \
"NoLinear/TimeApex.dat" u ($1*1E12):2 w l lc 3 lt 2 t '1ph Off' axis x1y1, \
"NoMPI/TimeApex.dat" u ($1*1E12):2 w l lc 3 lt 3 t '2ph Off' axis x1y1, \
"NoImpact/TimeApex.dat" u ($1*1E12):2 w l lc 3 lt 4 t 'Impact Off' axis x1y1, \
"NoIB/TimeApex.dat" u ($1*1E12):2 w l lc 3 lt 5 lw 3 t 'IB Off' axis x1y1, \
"NoDiffusion/TimeApex.dat" u ($1*1E12):2 w l lc 3 lt 6 t 'Diffusion Off' axis x1y1, \
"Normal/TimeApex.dat" u ($1*1E12):5 w l lc 1 lt 1 t '1ph + 2ph + Impact + IB' axis x1y2, \
"NoLinear/TimeApex.dat" u ($1*1E12):5 w l lc 1 lt 2 t '1ph Off' axis x1y2, \
"NoMPI/TimeApex.dat" u ($1*1E12):5 w l lc 1 lt 3 t '2ph Off' axis x1y2, \
"NoImpact/TimeApex.dat" u ($1*1E12):5 w l lc 1 lt 4 t 'Impact Off' axis x1y2, \
"NoIB/TimeApex.dat" u ($1*1E12):5 w l lc 1 lt 5 lw 3 t 'IB Off' axis x1y2, \
"NoDiffusion/TimeApex.dat" u ($1*1E12):5 w l lc 1 lt 6 t 'Diffusion Off' axis x1y2
