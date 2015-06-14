#! gnuplot

unset multiplot

reset
set size 1
set output '20150509-Evolution.eps'
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
set output '20130122-515nm-Apex.eps'
set terminal postscript eps enhanced color

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

scale=1e12

set terminal postscript eps enhanced color font 'Helvetica, 26' 
set output '20150509-TipHeating-Log.eps'

set multiplot 

set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'

unset log y
set key out top center horizontal spacing 1.2

# ### For Non-log scale
# # set xtics format "%0.0t{/Symbol \327}10^{%L}"
# set xtics format "%g" 1
# # set xtics format "%3.1e"
# set xrange [-0.5:5]
# unset log x


### For LOG scale
set log x
set xtics format "10^{%L}" 100
set ytics 100
set xrange [1e-4:10]

plot "TimeApex.dat" u ($1*scale):4 w l lc 1 lw 5 t 'T_{apex}, 515 nm', \
"TimeUp.dat" u ($1*scale):4 w l lc 3 lw 3 t 'T_{top}, 515 nm', \
"TimeBottom.dat" u ($1*scale):4 w l lc 4 lw 3 t 'T_{bottom}, 515 nm', \
"TimeMax.dat" u ($1*scale):4 w l lc 7 lw 3 t 'T_{max}, 515 nm'

## enable if you want an inset
# set xrange [1e-4:10]
# set origin 0.45, 0.25
# set size 0.5, 0.5
# unset xlabel 
# unset ylabel
# set log x
# set xtics format "10^{%L}" 100
# set ytics 100
# unset key
# replot

unset multiplot


# ########## plot mesh with points top bottom and apex
# set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
# set output "20150509-PositionOnMesh.dat"
# 
# M=2001
# N=151
# 
# set size 1.0,1.0
# scale=1E9
# 
# set key left top
# set view map
# 
# set xrange [-20:30]
# 
# splot "mesh.dat" u ($1*scale):($2*scale):(1) w l lw 0.4 notitle, \
# "< awk '{ if($3==1 && $4==76) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 6 t 'Apex', \
# "< awk '{ if($3==1 && $4==1) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 6 t 'Top', \
# "< awk '{ if($3==1 && $4==151) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 6 t 'Bottom'

########## Free carrier density at the tip apex
reset
unset multiplot
set output '20150509-ApexNeNh.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

set multiplot
set origin 0.0, 0.0
set size 1.0, 1.0

set xlabel 'Time (ps)'
set ylabel 'Carrier density (cm^{-3})'
set y2label 'Carrier temperature (K)'

set xrange [:5]

set xtics format '%g'
set ytics nomirror format '%3.1l x 10^{%L}' 2.5e20
set y2tics nomirror 

scale=1E12
yscale=1E6

shift=0e0; #2e-12

unset log x
# set log y
unset log y2
set key top right

plot "TimeApex.dat" u (scale*($1+shift)):($5/yscale) w l lw 5 lc 1 t 'N_e', \
"TimeApex.dat" u (scale*($1+shift)):($6/yscale) w l lw 5 lc 1 t 'N_h', \
"TimeApex.dat" u (scale*($1+shift)):2 w l lw 5 lc 3 t 'T_e' axis x1y2, \
"TimeApex.dat" u (scale*($1+shift)):3 w l lw 5 lc 3 t 'T_h' axis x1y2

# scale=1E12
# 
# set key bottom left
# set log x
# set xlabel 'Time (ps)'
# set xtics format '10^{%L}'
# # set output '20130711-ApexNeNhLog.eps'
# set xrange [1e-15*scale:10]
# replot

unset multiplot

########## Maximum free carrier density inside tip
reset
unset multiplot
set output '20150509-MaxNeNh.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

set multiplot
set origin 0.0, 0.0
set size 1.0, 1.0

set xlabel 'Time (ps)'
set ylabel 'Carrier density (cm^{-3})'
set y2label 'Carrier temperature (K)'

set xrange [:5]

set xtics format '%g'
set ytics nomirror format '%3.1l x 10^{%L}' 2.5e20
set y2tics nomirror 

scale=1E12
yscale=1E6

shift=0e0; #2e-12

unset log x
# set log y
unset log y2
set key top right

plot "TimeMax.dat" u (scale*($1+shift)):($5/yscale) w l lw 5 lc 1 t 'N_e', \
"TimeMax.dat" u (scale*($1+shift)):($6/yscale) w l lw 5 lc 1 t 'N_h', \
"TimeMax.dat" u (scale*($1+shift)):2 w l lw 5 lc 3 t 'T_e' axis x1y2, \
"TimeMax.dat" u (scale*($1+shift)):3 w l lw 5 lc 3 t 'T_h' axis x1y2