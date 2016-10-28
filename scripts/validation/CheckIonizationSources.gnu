#!gnuplot

reset

############## 
c=2.99792458e8; kb=1.3806504e-23; 
h=6.62e-34; hbar=h/(2e0*pi); 
ec=1.602176487e-19; 

lambda=1030E-9; fluence=100e0; tau=40e-15
ne=1E27; Te0=80e0
t0=5e0*tau
Ne0=1E27
sigmaT=tau/(2e0*sqrt(2e0*log(2e0)))
gap=1.16*ec
nu=1E15; #omega/3e0


# IR 1030 nm
epsilon={12.8e0,0.001414418e0}
sigma2=1.933288399e-11

# IR 800 nm
# sigma1=102024.792221655e0; sigma2=10e-11; 


# UV 343 nm
# epsilon={18.81766303e0,31.5464e0}
# sigma2=0e0

###
R=0e0 #abs((sqrt(1e0)-sqrt(epsilon))/(sqrt(1e0)+sqrt(epsilon)))**2; 

# formulas
omega=2e0*pi*c/lambda
I0=fluence/tau*sqrt(4*log(2e0)/pi)


sigma1=2e0*omega/c*imag(sqrt(epsilon))


Ce=3e0/2e0 * kb * ne
TeSig1(t) = 1e0/2e0*(-hbar*omega+gap)*sigma1*(-1e0+R)/hbar/omega*I0/Ce*pi**(1e0/2e0)*2**(1e0/2e0)*sigmaT*erf(1e0/2e0*2**(1e0/2e0)/sigmaT*t-1e0/2e0*t0/sigmaT*2**(1e0/2e0))+Te0+1e0/2e0*(-hbar*omega+gap)*sigma1*(-1e0+R)/hbar/omega*I0/Ce*pi**(1e0/2e0)*2**(1e0/2e0)*sigmaT*erf(1e0/2e0*t0/sigmaT*2**(1e0/2e0))

TeSig2(t)=-0.25e0*(-1e0+R)**2e0*sigma2*I0**2e0*(-2e0*hbar*omega+gap)*pi**0.5e0*sigmaT*erf(t/sigmaT-t0/sigmaT)/Ce/omega/hbar+Te0-0.25e0*(-1e0+R)**2e0*sigma2*I0**2e0*(-2e0*hbar*omega+gap)*pi**0.5e0*sigmaT*erf(t0/sigmaT)/Ce/omega/hbar
# TeSBDSig1(t)=4e0*sigma1*I0*(-hbar*omega+gap)/omega*(-1e0+R)/hbar*nu*pi**(0.5e0)*sigmaT**2e0*(erf(0.5e0*2e0**(0.5e0)/sigmaT*t-0.5e0*t0/sigmaT*2e0**(0.5e0))*(0.5e0*2e0**(0.5e0)/sigmaT*t-0.5e0*t0/sigmaT*2e0**(0.5e0))+1e0/pi**(0.5e0)*exp(-(0.5e0*2e0**(0.5e0)/sigmaT*t-0.5e0*t0/sigmaT*2e0**(0.5e0))**2e0))+2e0*sigma1*I0*nu*pi**(0.5e0)*2e0**(0.5e0)*sigmaT*erf(0.5e0*t0/sigmaT*2e0**(0.5e0))*(hbar*omega-hbar*omega*R-gap+gap*R)/hbar/omega*t-(-Te0*hbar*omega+4e0*sigma1*I0*nu*sigmaT**2e0*hbar*omega-4e0*sigma1*I0*nu*sigmaT**2e0*hbar*omega*R-4e0*sigma1*I0*nu*sigmaT**2e0*gap+4e0*sigma1*I0*nu*sigmaT**2e0*gap*R+2e0*sigma1*I0*nu*pi**(0.5e0)*2e0**(0.5e0)*sigmaT*erf(0.5e0*t0/sigmaT*2e0**(0.5e0))*t0*hbar*omega-2e0*sigma1*I0*nu*pi**(0.5e0)*2e0**(0.5e0)*sigmaT*erf(0.5e0*t0/sigmaT*2e0**(0.5e0))*t0*hbar*omega*R-2e0*sigma1*I0*nu*pi**(0.5e0)*2e0**(0.5e0)*sigmaT*erf(0.5e0*t0/sigmaT*2e0**(0.5e0))*t0*gap+2e0*sigma1*I0*nu*pi**(0.5e0)*2e0**(0.5e0)*sigmaT*erf(0.5e0*t0/sigmaT*2e0**(0.5e0))*t0*gap*R)/hbar/omega

############## TEMPERATURE ##############
#
### One photon transition check
#shift=500e-15-3.1e-13+1E-14
##
#set xlabel 'Time (s)'
#set ylabel 'Relative error (%)'
#set xrange [1e-14:1E-12]
#set grid
#set log x
#unset log y
#set key out center bottom
## plot  "TimeMax.dat" u ($1+shift):(100e0*($2-TeSig1($1+shift))/TeSig1($1+shift)) w p lw 1 t 'Relative error (%)'
#
### Two photon transition check
#plot  "TimeMax.dat" u ($1+shift):(100e0*($2-TeSig2($1+shift))/TeSig2($1+shift)) w p lw 1 t 'Relative error (%)'
#
#
## TeSig1(x) t 'Analytic, 1 photon transition' w l lw 2, \
## TeSig2(x) t 'Analytic, 2 photon transition', \
## 6.70481583500663e5 w l t 'Asymptotic limit for 2 photon transition', \
## "TimeSig2-2e15.dat" u ($2+shift):3 w p lw 2 t 'Numerical, 2 photon transition'
#
## 3.35e17 w l t 'Asymptotic limit for 2 photon transition', \
#
## plot TeSBDSig1(x) t 'Analytic, 1 photon transition' w l lw 2, \
#
###
#set output '20130905-HeatingIonization2-1E27.eps'
#set terminal postscript eps enhanced monochrome font 'Helvetica, 26'
#replot


########### DENSITY ##############

Ne0=1E5

 set output '20130907-NePhotonic.eps'
 set terminal postscript eps enhanced color font 'Helvetica, 26'
 
 Ne1(t)=-sigma1 / hbar * I0 / omega * (-1e0 + R) * sqrt(0.3141592654e1) * sqrt(0.2e1) * sigmaT * erf(t / sigmaT * sqrt(0.2e1) / 0.2e1) / 0.2e1 #* sqrt(4e0*log(2e0)/pi)
 Ne2(t)=0.1e1 / hbar * I0 ** 2 / omega * sigma2 * ((-1 + R) ** 2) * sqrt(0.3141592654e1) * sigmaT * erf(t / sigmaT) / 0.4e1
 
 set key outside bottom center
 set xlabel 'Time (s)'
 set ylabel 'Density (m3)'
 set xrange [0e0:10*tau]
 
# plot Ne1(x-t0)+Ne1(t0)+Ne0 w l t '1 photon, Analytic', \
# "TimeNeSig1.dat" u ($2+shift):5 w p t '1 photon, Numeric', \
# Ne2(x-t0)+Ne2(t0)+Ne0 w l t '2 photons, Analytic', \
# "TimeNeSig2.dat" u ($2+shift):5 w p t '2 photons, Numeric'
shift=-2e-13
t0=2e-13

plot "TimeMax.dat" u ($1-shift):(100e0*($5-(Ne2(($1-shift)-t0)+Ne2(t0)+Ne0))/((Ne2(($1-shift)-t0)+Ne2(t0)+Ne0))) w l t 'Relative error (%)'

plot "TimeMax.dat" u ($1-shift):(Ne2(($1-shift)-t0)+Ne2(t0)+Ne0) w l t 'Analytic', \
"TimeMax.dat" u ($1-shift):5 w p t 'Numeric'
