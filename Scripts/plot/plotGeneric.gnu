#! gnuplot

unset multiplot

reset
set size 1
set output '20120828-Evolution.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"
set xrange [:]
set xlabel 'Time (ps)'
set ylabel 'T [K]'
set y2label 'N [m^{-3}]'
set ytics nomirror format "%g"
set y2tics nomirror #format "%3.1e"
# set format y2 "10^{%L}"
set format y "%2.0t{/Symbol \327}10^{%L}"
set xtics format "%g"
unset log x
set log y
set log y2
set grid
set key outside below vertical
scale=1E12
plot "TimeMax.dat" u ($1*scale):2 w l t 'Te', \
"TimeMax.dat" u ($1*scale):3 w l t 'Th', \
"TimeMax.dat" u ($1*scale):4 w l t 'Ts', \
"TimeMax.dat" u ($1*scale):5 w l t 'Ne' axis x1y2, \
"TimeMax.dat" u ($1*scale):6 w l t 'Nh' axis x1y2, \
"TimeMax.dat" u ($1*scale):7 w l t 'Laser Intensity' axis x1y2 

# set size 1.0
# set terminal x11 enhanced 1
# replot

########## Tip apex datas #######
reset
# 
# set size 0.7	
# set output '20130122-515nm-Apex.eps'
# set terminal postscript eps enhanced color

set format "%g"

set xrange [:]
set xlabel 'Time (ps)'
set ylabel 'T [K]'
set y2label 'N [m^{-3}]'
set ytics nomirror format "%g"
set y2tics nomirror #format "%3.1e"
# set format y2 "10^{%L}"
set format y "%2.0lx10^{%L}"
set xtics format "%g"
unset log y
unset log y2
set grid
set key below left vertical

plot "TimeApex.dat" u ($1*1E12):2 w l t 'Te', \
"TimeApex.dat" u ($1*1E12):3 w l t 'Th', \
"TimeApex.dat" u ($1*1E12):4 w l t 'Ts', \
"TimeApex.dat" u ($1*1E12):5 w l t 'Ne' axis x1y2, \
"TimeApex.dat" u ($1*1E12):6 w l t 'Nh' axis x1y2, \
"TimeApex.dat" u ($1*1E12):7 w l t 'Laser Intensity' axis x1y2 

######## Heating of the lattice

reset

scale=1e9

set terminal postscript eps enhanced color font 'Helvetica, 26' 
set output '20121128-TipHeating-Log.eps'

set multiplot 

set xlabel 'Time (ns)'
set ylabel 'Temperature (K)'
# set xtics format "%0.0t{/Symbol \327}10^{%L}"
set xtics format "%g" 5
# set xtics format "%3.1e"
unset log x
unset log y
set key out top center horizontal spacing 1.2
set xrange [-0.5:]

plot "TimeApex.dat" u ($1*scale):4 w l lw 3 t 'T_{apex}, 515 nm', \
"TimeUp.dat" u ($1*scale):4 w l lw 3 t 'T_{top}, 515 nm', \
"TimeBottom.dat" u ($1*scale):4 w l lw 3 t 'T_{bottom}, 515 nm'
# "TimeMax.dat" u ($1*scale):4 w l lw 3 t 'T_{max}, 515 nm', \

set xrange [1e-4:10]
set origin 0.45, 0.25
set size 0.5, 0.5
unset xlabel 
unset ylabel
set log x
set xtics format "10^{%L}" 100
set ytics 100
unset key
replot

unset multiplot


########## plot mesh with points top bottom and apex
set output '20130722-PositionOnMesh.dat'
set terminal postscript eps enhanced monochrome
set size 1.0,1.0
scale=1E9
set view map
splot "mesh.dat" u ($1*scale):($2*scale):(1) w l lw 0.4 notitle, \
"< awk '{ if($3==1 && $4==25) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 6 t 'Apex', \
"< awk '{ if($3==1 && $4==51) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 6 t 'Top', \
"< awk '{ if($3==1 && $4==1) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 6 t 'Bottom'

########## Free carrier density at the tip apex
reset
set output '20130711-ApexNeNh.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

set multiplot
set origin 0.0, 0.0
set size 1.0, 1.0
set xlabel 'Time (ps)'
set ylabel 'Carrier density (m^{-3})'
set y2label 'Carrier temperature (m^{-3})'
set xrange [:20]
set xtics format '%g'
set ytics nomirror format '10^{%L}'
set y2tics nomirror 
scale=1E12
shift=2e-12
unset log x
set log y
set log y2
set key bottom right
plot "TimeApex.dat" u (scale*($1+shift)):5 w l lc 1 t 'N_e', \
"TimeApex.dat" u (scale*($1+shift)):6 w l lc 1 t 'N_h', \
"TimeApex.dat" u (scale*($1+shift)):2 w l lc 3 t 'T_e' axis x1y2, \
"TimeApex.dat" u (scale*($1+shift)):3 w l lc 3 t 'T_h' axis x1y2

scale=1E9

set key bottom left
set log x
set xlabel 'Time (ns)'
set xtics format '10^{%L}'
set output '20130711-ApexNeNhLog.eps'
set xrange [1e-15*scale:10]
replot
