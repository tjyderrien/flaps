#! gnuplot

kb=1.38e-23

reset
set terminal x11 1 enhanced
set title "eta(Te,Ne)"
set xlabel 'Te [K]'
set ylabel 'Ne [m^{-3}]'
set log xy 
set yrange [1e22:]
set zrange [-50:100]
splot "Ce.dat" u 2:1:3 w d

# set ztics format "%3.0e"
reset
set terminal x11 2 enhanced 
set title 'Ce(Te,Ne)'
set xlabel 'Te [K]'
set ylabel 'Ne [m^{-3}]'
set zlabel 'C_{e}/(k_B n_e)'
# set log y
splot "Ce.dat" u 2:1:6 w d

reset
set size 1.0
unset log xy
set terminal x11 3 enhanced
set title 'eta(Te;Ne)'
set xlabel 'Te [K]'
set ylabel '{/Symbol h}'
set yrange [-10:100]
plot "Ce.dat" u 2:3 notitle
set output '20120703-EtaChenModel-Corrected.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

########### ELECTRONS 

##### plot Ce(Te,ne)

reset
set key inside right bottom
set terminal x11 4 enhanced
set grid
set log x
set title 'Ce(Te;Ne)/k_B n_e'
set xlabel 'Te [K]'
set ylabel 'Ce/(k_B n_e)'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E21) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{26} m^{-3}'

# "< awk '{ if($1==1E27) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{28} m^{-3}' 

# "< awk '{ if($1==1E19) print }' Ce.dat" u 2:6 w l t 'n_e=10^{19} m^{-3}', \
# "< awk '{ if($1==1E20) print }' Ce.dat" u 2:6 w l t 'n_e=10^{20} m^{-3}', \

set output '20120703-CeChenModel-Corrected.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

#### plot Ke(Te,ne)
reset
set terminal x11 5 enhanced
set size 1.0
set grid
set key outside right
set log xy
set title 'Ke(Te;Ne)'
set xlabel 'Te [K]'
set ylabel 'Ke [W.m^{-1}.K^{-1}]'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E19) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{19} m^{-3}', \
"< awk '{ if($1==1E20) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{20} m^{-3}', \
"< awk '{ if($1==1E21) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{26} m^{-3}'
# "< awk '{ if($1==1E27) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{28} m^{-3}' 

set output '20120703-KeChenModel.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

### plot De(Te,ne)
reset
set terminal x11 6 enhanced
set grid
set key outside right
set log xy
set title 'De(Te;Ne)'
set xlabel 'Te [K]'
set ylabel 'De [m^{2}.s^{-1}]'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E19) print }' Ce.dat" u 2:8 w l t 'n_e=10^{19} m^{-3}', \
"< awk '{ if($1==1E20) print }' Ce.dat" u 2:8 w l t 'n_e=10^{20} m^{-3}', \
"< awk '{ if($1==1E21) print }' Ce.dat" u 2:8 w l t 'n_e=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ce.dat" u 2:8 w l t 'n_e=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ce.dat" u 2:8 w l t 'n_e=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ce.dat" u 2:8 w l t 'n_e=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ce.dat" u 2:8 w l t 'n_e=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ce.dat" u 2:8 w l t 'n_e=10^{26} m^{-3}'
# "< awk '{ if($1==1E27) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{28} m^{-3}' 

set output '20120703-DeChenModel.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

reset
set terminal x11 9 enhanced
set grid
set key outside right
set log xy
set title 'De(Te;Ne) Chen expression'
set xlabel 'Te [K]'
set ylabel 'De [m^{2}.s^{-1}]'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E19) print }' Ce.dat" u 2:9 w l t 'n_e=10^{19} m^{-3}', \
"< awk '{ if($1==1E20) print }' Ce.dat" u 2:9 w l t 'n_e=10^{20} m^{-3}', \
"< awk '{ if($1==1E21) print }' Ce.dat" u 2:9 w l t 'n_e=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ce.dat" u 2:9 w l t 'n_e=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ce.dat" u 2:9 w l t 'n_e=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ce.dat" u 2:9 w l t 'n_e=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ce.dat" u 2:9 w l t 'n_e=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ce.dat" u 2:9 w l t 'n_e=10^{26} m^{-3}'
# "< awk '{ if($1==1E27) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ce.dat" u 2:7 w lp t 'n_e=10^{28} m^{-3}' 

set output '20120703-DeChenExpression.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot


########### HOLES 

##### plot Ch(Th,nh)

reset
set key inside right bottom
set terminal x11 4 enhanced
set grid
set log x
set title 'Ch(Th;Nh)/k_B n_h'
set xlabel 'Th [K]'
set ylabel 'Ch/(k_B n_h)'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E21) print }' Ch.dat" u 2:6 w l lw 2 t 'n_h=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ch.dat" u 2:6 w l lw 2 t 'n_h=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ch.dat" u 2:6 w l lw 2 t 'n_h=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ch.dat" u 2:6 w l lw 2 t 'n_h=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ch.dat" u 2:6 w l lw 2 t 'n_h=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ch.dat" u 2:6 w l lw 2 t 'n_h=10^{26} m^{-3}'

# "< awk '{ if($1==1E27) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ce.dat" u 2:6 w l lw 2 t 'n_e=10^{28} m^{-3}' 

# "< awk '{ if($1==1E19) print }' Ce.dat" u 2:6 w l t 'n_e=10^{19} m^{-3}', \
# "< awk '{ if($1==1E20) print }' Ce.dat" u 2:6 w l t 'n_e=10^{20} m^{-3}', \

set output '20120703-ChChenModel-Corrected.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

#### plot Kh(Th,nh)
reset
set terminal x11 7 enhanced
set grid
set key outside right
set log xy
set title 'Kh(Th;Nh)'
set xlabel 'Th [K]'
set ylabel 'Kh [W.m^{-1}.K^{-1}]'
# set yrange [0.01:100]
# plot "< awk '{ if($1==1E19) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{19} m^{-3}', \
# "< awk '{ if($1==1E20) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{20} m^{-3}', \
plot "< awk '{ if($1==1E21) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ch.dat" u 2:7 w lp t 'n_h=10^{26} m^{-3}'

#�"< awk '{ if($1==1E27) print }' Ch.dat" u 2:7 w lp t 'n_e=10^{27} m^{-3}' 
#�"< awk '{ if($1==1E28) print }' Ch.dat" u 2:7 w lp t 'n_e=10^{28} m^{-3}' 

set output '20120703-KhChenModel.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

### plot Dh(Th,nh)
reset
set terminal x11 8 enhanced
set grid
set key outside right
set log xy
set title 'Dh(Th;Nh)'
set xlabel 'Th [K]'
set ylabel 'Dh [m^{2}.s^{-1}]'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E19) print }' Ch.dat" u 2:8 w l t 'n_h=10^{19} m^{-3}', \
"< awk '{ if($1==1E20) print }' Ch.dat" u 2:8 w l t 'n_h=10^{20} m^{-3}', \
"< awk '{ if($1==1E21) print }' Ch.dat" u 2:8 w l t 'n_h=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ch.dat" u 2:8 w l t 'n_h=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ch.dat" u 2:8 w l t 'n_h=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ch.dat" u 2:8 w l t 'n_h=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ch.dat" u 2:8 w l t 'n_h=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ch.dat" u 2:8 w l t 'n_h=10^{26} m^{-3}'
# "< awk '{ if($1==1E27) print }' Ch.dat" u 2:7 w lp t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ch.dat" u 2:7 w lp t 'n_e=10^{28} m^{-3}' 

set output '20120703-DhChenModel.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot


### plot Dh(Th,nh) Chen expression
reset
set terminal x11 10 enhanced
set grid
set key outside right
set log xy
set title 'Dh(Th;Nh) Chen expression'
set xlabel 'Th [K]'
set ylabel 'Dh [m^{2}.s^{-1}]'
# set yrange [0.01:100]
plot "< awk '{ if($1==1E19) print }' Ch.dat" u 2:9 w l t 'n_h=10^{19} m^{-3}', \
"< awk '{ if($1==1E20) print }' Ch.dat" u 2:9 w l t 'n_h=10^{20} m^{-3}', \
"< awk '{ if($1==1E21) print }' Ch.dat" u 2:9 w l t 'n_h=10^{21} m^{-3}', \
"< awk '{ if($1==1E22) print }' Ch.dat" u 2:9 w l t 'n_h=10^{22} m^{-3}', \
"< awk '{ if($1==1E23) print }' Ch.dat" u 2:9 w l t 'n_h=10^{23} m^{-3}', \
"< awk '{ if($1==1E24) print }' Ch.dat" u 2:9 w l t 'n_h=10^{24} m^{-3}', \
"< awk '{ if($1==1E25) print }' Ch.dat" u 2:9 w l t 'n_h=10^{25} m^{-3}', \
"< awk '{ if($1==1E26) print }' Ch.dat" u 2:9 w l t 'n_h=10^{26} m^{-3}'
# "< awk '{ if($1==1E27) print }' Ch.dat" u 2:7 w lp t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ch.dat" u 2:7 w lp t 'n_e=10^{28} m^{-3}' 

set output '20120703-DhChenExpression.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot

###### Compare coupling electron - lattice with and without Fermi-Dirac statistics

reset
set terminal x11 enhanced
set size 1
set grid
set key bottom right
set log xy
set title 'Coupling rate (W.m^{-3}.K^{-1})'
set xlabel 'Ne (K)'
set ylabel 'Te (K)'

time0=240e-15
ne0=6.02e26
time(ne)=time0*(1e0+(ne/ne0)**2e0) #Sjodin PRL 

kb=1.38e-23

# set yrange [0.01:100]
# plot "< awk '{ if($1==1E19) print }' Ce.dat" u ($2/:6 w lp t 'n_e=10^{19} m^{-3}', \
# "< awk '{ if($1==1E20) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{20} m^{-3}', \
# "< awk '{ if($1==1E21) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{21} m^{-3}', \
# "< awk '{ if($1==1E22) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{22} m^{-3}', \
# "< awk '{ if($1==1E23) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{23} m^{-3}', \
# "< awk '{ if($1==1E24) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{24} m^{-3}', \
# "< awk '{ if($1==1E25) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{25} m^{-3}', \
# "< awk '{ if($1==1E26) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{26} m^{-3}'
# "< awk '{ if($1==1E27) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{27} m^{-3}' 
# "< awk '{ if($1==1E28) print }' Ce.dat" u 2:6 w lp t 'n_e=10^{28} m^{-3}' 

# splot "Ce.dat" u 1:2:($6/time($1)) w p t 'Coupling rate (W.m^{-3}.K^{-1})', \
set xrange [:1e28]
set xlabel 'Ne (m^{-3})'
set ylabel 'Coupling rate (W.m^{-3}.K^{-1})'

plot "< awk '{ if($2==81) print }' Ce.dat" u 1:($6/time($1)) w l lw 2 smooth bezier t 'Sjodin (1998), Fermi Ce', \
1.5e0*kb*x/time(x) w l lw 2 t 'Sjodin (1998), classical Ce', \
1.5e0*kb*x/time0 w l lw 2 t 'Constant coupling time, t0=240 fs'

set output '20120703-CeCoupling-Corrected.eps'
set terminal postscript eps enhanced monochrome font 'Helvetica, 24'
replot



###### Affiche les valeurs du fichier pour export des fonctions de Fermi vers le code

reset
set terminal x11 1 enhanced
set format "%g"
set grid
set key outside right
unset log x
set log y
unset log y2
set y2tics nomirror
set ytics nomirror
set xlabel 'eta'
set ylabel 'Ne/Nc'
set y2label 'Alog10(Ne/Nc)+B'
# set multiplot layout 1,2
plot "Ce.dat" u 3:10 w p t 'Ne/Nc real' axis x1y1 ,\
     "FermiDatas.dat" u 4:3 w p lw 1 t 'Ne/Nc tabulated' axis x1y1, \
     "FermiDatas.dat" u 4:1 w p t 'A log(Ne/Nc) + B' axis x1y2
     
     
