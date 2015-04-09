#!gnuplot

set log x
set ytics nomirror
set y2tics nomirror


##### Comparison of the effect of the mesh size on electron temperature evolution

plot "./251x251-10-NoDiffusion-dt1e-16s/TimeApex.dat" u 1:2 w l t 'Te, 251-10, dt=1e-16', \
"./251x251-30-NoDiffusion-dt1e-16s/TimeApex.dat" u 1:2 w l t 'Te, 251-30, dt=1e-16', \
"351x351-40-NoDiffusion-dt1e-16s-CV1e-8/TimeApex.dat" u 1:2 w l t 'Te, 351-40, dt=1e-16, CV=1E-8', \
"451x451-30-NoDiffusion-dt1e-17/TimeApex.dat" u 1:2 w l t 'Te, 451-30, dt=1e-17', \
"651x51-20-dt1e-16-CV1e-8/TimeApex.dat" u 1:2 w l t 'Te, 651x51-20, dt=1e-16, CV=1E-8', \
"1201x101-20-dt1e-16-CV1e-8/TimeApex.dat" u 1:2 w l t 'Te, 1201x101-20, dt=1e-16, CV=1E-8', \
"2001x101-20-dt1e-16-CV1e-10/TimeApex.dat" u 1:2 w l t 'Te, 2001x101-20, dt=1e-16, CV=1E-10', \
"2001x151-20-dt1e-17-CV1e-10/TimeApex.dat" u 1:2 w l t 'Te, 2001x151-20, dt=1e-17, CV=1E-10', \
"/home/thibault/Documents/Codes/Flaps2D/MieScatteringVersion/TimeApex.dat" u 1:2 w l t 'TeDiff, 2001x151-20, dt=1e-17, CV=1e-10'


##### Comparison of the effect of the mesh size on lattice temperature evolution

plot "./251x251-10-NoDiffusion-dt1e-16s/TimeApex.dat" u 1:4 w l t 'Ts, 251-10, dt=1e-16', \
"./251x251-30-NoDiffusion-dt1e-16s/TimeApex.dat" u 1:4 w l t 'Ts, 251-30, dt=1e-16', \
"351x351-40-NoDiffusion-dt1e-16s-CV1e-8/TimeApex.dat" u 1:4 w l t 'Ts, 351-40, dt=1e-16, CV=1E-8', \
"451x451-30-NoDiffusion-dt1e-17/TimeApex.dat" u 1:4 w l t 'Ts, 451-30, dt=1e-17', \
"651x51-20-dt1e-16-CV1e-8/TimeApex.dat" u 1:4 w l t 'Ts, 651x51-20, dt=1e-16, CV=1E-8', \
"1201x101-20-dt1e-16-CV1e-8/TimeApex.dat" u 1:4 w l t 'Ts, 1201x101-20, dt=1e-16, CV=1E-8', \
"2001x101-20-dt1e-16-CV1e-10/TimeApex.dat" u 1:4 w l t 'Ts, 2001x101-20, dt=1e-16, CV=1E-10', \
"2001x151-20-dt1e-17-CV1e-10/TimeApex.dat" u 1:4 w l t 'Ts, 2001x151-20, dt=1e-17, CV=1e-10', \
"/home/thibault/Documents/Codes/Flaps2D/MieScatteringVersion/TimeApex.dat" u 1:4 w l t 'Ts, 2001x151-20, dt=1e-18, CV=1e-10'

##### Comparison of the effect of the mesh size on density carrier evolution

plot "./251x251-10-NoDiffusion-dt1e-16s/TimeApex.dat" u 1:5 w l t 'Ne, 251-10, dt=1e-16' axis x1y2, \
"./251x251-30-NoDiffusion-dt1e-16s/TimeApex.dat" u 1:5 w l t 'Ne, 251-30, dt=1e-16' axis x1y2, \
"351x351-40-NoDiffusion-dt1e-16s-CV1e-8/TimeApex.dat" u 1:5 w l t 'Ne, 351-30, dt=1e-16, CV=1E-8' axis x1y2, \
"451x451-30-NoDiffusion-dt1e-17/TimeApex.dat" u 1:5 w l axis x1y2 t 'Ne, 451-30, dt=1e-17', \
"651x51-20-dt1e-16-CV1e-8/TimeApex.dat" u 1:5 w l t 'Ne, 651x51-20, dt=1e-16, CV=1E-8' axis x1y2, \
"1201x101-20-dt1e-16-CV1e-8/TimeApex.dat" u 1:5 w l t 'Ne, 1201x101-20, dt=1e-16, CV=1E-8' axis x1y2, \
"2001x101-20-dt1e-16-CV1e-10/TimeApex.dat" u 1:5 w l t 'Ne, 2001x101-20, dt=1e-16, CV=1E-10' axis x1y2, \
"2001x151-20-dt1e-17-CV1e-10/TimeApex.dat" u 1:5 w l t 'Ne, 2001x151-20, dt=1e-17, CV=1e-10' axis x1y2, \
"/home/thibault/Documents/Codes/Flaps2D/MieScatteringVersion/TimeApex.dat" u 1:5 w l t 'NeDiff, 2001x151-20, dt=1e-17, CV=1e-10' axis x1y2

##### Comparison of the meshes

unset log x
unset y2tics
set view map
splot "/home/thibault/Documents/Codes/Flaps2D/MieScatteringVersion/mesh.dat" u 1:2:3 w d t '1201-101, eps=1e-8'

# "351x351-40-NoDiffusion-dt1e-16s-CV1e-8/mesh.dat" u 1:2:3 w d t '351-40, eps=1e-8', \

###### Visualization of the field distribution

unset log x
unset y2tics
set view 30,30
splot "/home/thibault/Documents/Codes/Flaps2D/MieScatteringVersion/Field.dat" u ($1*1E6):($2*1E6):(($3**2e0)**0.5e0) w d t 'Field'
