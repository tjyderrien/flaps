#!gnuplot

##### IN THIS FILE: all 2D meshes for Ultrasonde. 

#! gnuplot

##### plot mesh 
#scale=1E9
#set origin 0.45, 0.3
#set size 0.5,0.5
#set view map
#set xlabel 'X (nm)'
#set ylabel 'Y (nm)'
#set xtics 250e-9*scale
#set ytics 200e-9*scale
#set xrange [-0.01e-6*scale:scale*0.5e-6]
#set yrange [-0.2e-6*scale:0.2e-6*scale]
#unset key
#splot "mesh.dat" u ($1*scale):($2*scale):(1) w l lw 0.1 notitle, \
#"< awk '{ if($3==1 && $4==26) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 10 t 'Apex', \
#"< awk '{ if($3==1 && $4==51) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 10 t 'Top', \
#"< awk '{ if($3==1 && $4==1) print }' mesh.dat" u ($1*scale):($2*scale):(1) w p lw 10 t 'Bottom'

# unset multiplot


####### map of heating time
system "tail -n1 TimeApex.dat | awk '{ print $1 }'"
reset
set output '20120905-HeatingTimeMap.eps'
set terminal postscript eps enhanced color

set pm3d map interpolate 4,4
set format "%g"
set xlabel 'X [{/Symbol m}m]'
set ylabel 'Y [{/Symbol m}m]'
# set palette defined (0 "blue", 1 "white", 2 "red")
# set palette rgb 33,13,10 #rainbow
set palette rgb 30,31,32 #many colors
set size square
set xrange [:]
set yrange [:]
set log cb
splot "< awk '{if($1==0.59800E-11) print}' Depth.dat | awk -f '/home/thibault/Documents/Codes/Scripts/plotXZ.awk'" u ($2*1E6):($3*1E6):14 w pm3d t columnheader 1


set output '20120905-HeatingMaxMap.eps'
unset log cb
splot "< awk '{if($1==0.59800E-11) print}' Depth.dat | awk -f '/home/thibault/Documents/Codes/Scripts/plotXZ.awk'" u ($2*1E6):($3*1E6):15 w pm3d t columnheader 1

system "gv '20120905-HeatingMaxMap.eps' &"
system "gv '20120905-HeatingTimeMap.eps' &"


#### simple map of potential
system "tail -n1 TimeApex.dat | awk '{ print $1 }'"
reset

set terminal x11 1 enhanced
set pm3d map interpolate 4,4
set format "%g"
set xlabel 'X [{/Symbol m}m]'
set ylabel 'Y [{/Symbol m}m]'
set palette rgb 30,31,32 #many colors
set size square
set xrange [:]
set yrange [:]
set log cb
splot "< awk '{if($1==-0.19600E-12) print}' DepthVessel.dat | awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk'" u ($2*1E6):($3*1E6):4 w pm3d t columnheader 1

set terminal x11 2 enhanced
unset log cb
splot "< awk '{if($1==-0.19600E-12) print}' DepthVessel.dat | awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk'" u ($2*1E6):($3*1E6):5 w pm3d t columnheader 1

set terminal x11 3 enhanced
splot "< awk '{if($1==-0.19600E-12) print}' DepthVessel.dat | awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk'" u ($2*1E6):($3*1E6):6 w pm3d t columnheader 1


########################
### NEEDLE MESH ####
########################
# set xrange [:3e-6]

reset
set format "%g"

set output '20121018-MeshNeedle.eps'
set terminal postscript eps enhanced monochrome

set size ratio -1 0.8
set title 'Mesh'
set xlabel 'X (um)'
set ylabel 'Y (um)'
set grid
set view map
splot "mesh.dat" u ($1*1E6):($2*1E6):3 w l lw 2 t 'Needle'

set output
set terminal x11 enhanced
set size ratio -1 1.0
replot

########################
### Vacuum vessel MESH ####
########################
# set xrange [:3e-6]

reset
set format "%g"
set output '20121018-MeshVessel.eps'
set terminal postscript eps enhanced monochrome

set size ratio -1 0.8
set title 'Mesh'
set xlabel 'X (um)'
set ylabel 'Y (um)'
set grid
set view map
splot "meshVessel.dat" u ($1*1E6):($2*1E6):3 w l lw 1 notitle

set output
set terminal x11 enhanced
set size ratio -1 1.0
replot

#######################
### MESHES TOGETHER ###
#######################
reset
set format "%g"
set output '20121018-MeshBoth.eps'
set terminal postscript eps enhanced monochrome
set key opaque box outside right
set size 0.8
set title 'Mesh'
set xlabel 'X (um)'
set ylabel 'Y (um)'
set grid
set view map
set hidden
set xrange [-5:5]
splot "meshVessel.dat" u ($1*1E6):($2*1E6):3 w l lt 1 lw 1 t 'Vessel', \
"mesh.dat" u ($1*1E6):($2*1E6):3 w l lt 1 lw 1 t 'Needle'
set output
set terminal x11 enhanced
set size 1.0
replot

############## plot effect of the induction ############"
reset
set size 0.8
set output '20121214-ComparisonOfInduction.eps'
set terminal postscript eps enhanced color

set ytics nomirror
set y2tics nomirror

set xlabel 'Time (fs)'
set y2label 'E (V/m)'
set ylabel 'N_{e,h} (m^{-3})'

plot "TimeApex.dat" u ($1*1E15+200):5 w l lw 3 t 'Ne', \
"TimeApex.dat" u ($1*1E15+200):6 w l lw 3 t 'Nh', \
"TimeApexField.dat" u ($1*1E15+200):2 w l lw 3 axis x1y2 t 'Ex', \
"TimeApexInducedField.dat" u ($1*1E15+200):2 w l lw 3 axis x1y2 t 'Ex with induction'

###################

set title 'Meshes'
splot "mesh.dat" u 1:2:3 w l t 'Needle', \
"meshVessel.dat" u 1:2:3 w l t 'Vessel'

set title 'Density'
splot "88318.vid" u 2:3:9 w p lc rgbcolor "blue" t 'Nh', \
"88318.vid" u 2:3:8 w p lc rgbcolor "red" t 'Ne'

set title 'Dielectric constant'
splot "DepthVessel.map" u 2:3:7 w p lc 3 t 'epsilon_r (Vessel)', \
"mesh.dat" u 1:2:(0) w l lc 4 t 'Mesh needle'

set title 'Potential'
splot "DepthVessel.dat" u 2:3:4 w p lw 1 lc 1, \
"Depth.dat" u 2:3:34 w p lw 1 lc 3

set title 'Ex Field'
unset log xy
unset log z
splot "DepthVessel.dat" u 2:3:(($5**2)**.5) w p, \
"Depth.dat" u 2:3:(($35**2)**0.5) w p lc 3

set title 'Ey Field'
set log z
	splot "DepthVessel.dat" u 2:3:(($6**2)**.5) w p, \
"Depth.dat" u 2:3:(($35**2)**0.5) w p lc 3

set title 'Field amplitude'
set log z
unset log xy
splot "DepthVessel.dat" u 2:3:(($5**2+$6**2)**.5) w p, \
"Depth.dat" u 2:3:(($35**2+$36**2)**0.5) w p lc 3

set log cb
set format "%g"
set log z
set pm3d map interpolate 4,4
splot "< awk '{if($1==-1.617e-13 && $2<2e-6) print}' DepthVessel.dat | awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk'" u (-$2*1E6):($3*1E6):(($5**2+$6**2)**.5) w pm3d t columnheader 1

unset log cb
splot "< awk '{if($1==2.19e-13 && $2<2e-6) print}' DepthVessel.dat | awk -f '/home/thibault/Documents/LaAPT/Scripts/plotXZ.awk'" u (-$2*1E6):($3*1E6):(($4**2)**.5) w pm3d t columnheader 1
