#!/usr/bin/env python
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker

rc('text', usetex=True)
rc('font', family='serif', size='18')

columnNumX=1
columnNumY=2
columnNumZ=3

Ny=101 #nombre de points en profondeur
# Ny=101 #nombre de points en profondeur

NxReg=2000; NyReg=2000; #mesh for interpolation
isolines=6

# isolevels=[1,2] #temperature
# isolevels=[1e5, 1e6, 2e6, 3e6, 4e6, 5e6]
# isolevels=[0,0.2,0.4,0.6,0.8,1] #enthalpy
# isolevels=[1,10,100,1E3,1E4,1E5,1E6,1E7,1E8,1E9,2E9,3E9,4E9,5E9,1E10] #log
# isolevels=[1E3,5E3,10e3,15e3]
isolevels=[3]
# isolevels=[1E10,1E25,1E27]

#### read file
f = open('Field-515nm-TE-0.5pi.map','r')
data = f.readlines()
N =len(data)
print N
x = np.zeros(N)
y = np.zeros(N)
z  = np.zeros(N)
for i in range(N):
    values = data[i].split()
    x[i] = float(values[columnNumX-1])
    y[i] = float(values[columnNumY-1])
    z[i] = float(values[columnNumZ-1])


# Ny=N/Nx
Nx=N/Ny

k=0

#### define boundaries
xmin=np.amin(x); xmax=np.amax(x)
ymin=np.amin(y); ymax=np.amax(y)

if xmax > ymax :
	ymax=xmax/2.
	ymin=-xmax/2.

# xmin=x[0]; xmax=x[Nx*Ny-1]
# ymin=y[0]; ymax=y[Ny*Nx-1]

#xmin=-1.; xmax=100.;
#ymin=0.; ymax=1; 

print "Points=", Nx, Ny
print xmin, xmax, ymin, ymax

### define range of the ticks
xminV=xmin; xmaxV=xmax
yminV=ymin; ymaxV=ymax

#### define mesh ticks
xn=np.arange(xminV, xmaxV, (xmaxV-xminV)/5.)
yn=np.arange(yminV, ymaxV, (ymaxV-yminV)/5.)

xS=np.zeros([Nx,Ny])
yS=np.zeros([Nx,Ny])
zn=np.zeros([Nx,Ny])
zn2=np.zeros([Nx,Ny])
zreg=np.zeros([NxReg,NyReg])
zreg2=np.zeros([NxReg,NyReg])

#### from series to matrix
for j in range(Ny):
	for i in range(Nx):
		xS[i,j]=x[k]
		yS[i,j]=y[k]
	        zn[i,j]=z[k]
	        k=k+1

#### PLOT ON REGULAR MESH 
#imgplot = plt.imshow(zn.T, extent=[xmin,xmax,ymin,ymax], aspect='auto')
## imgplot.set_interpolation('nearest')
# imgplot.set_interpolation('bicubic')
# imgplot.set_cmap('RdBu')
#
#print xn
#print yn
#
#plt.xticks(xn)
#plt.yticks(yn)
#
#plt.xlim(xmin,xmax)
#plt.ylim(ymin,ymax)
#
#plt.xlabel("Time (ps)")
#plt.ylabel("Depth (nm)")
#plt.colorbar()
#for j in range(Ny):
#	for i in range(Nx):
#	        zn2[i,j]=zn[i,Ny-j-1]
#CS=plt.contour(zn2,extent=[xmin,xmax,ymin,ymax],colors='black')
#plt.clabel(CS, inline=1, fontsize=12)
#plt.show()
#plt.savefig('plotted')

#### PLOT SCATTERED CONTOURS

### lets change the irregular mesh into a regular one
xi=np.linspace(xmin,xmax,NxReg)
yi=np.linspace(ymin,ymax,NyReg)
Xreg,Yreg=np.meshgrid(xi,yi)
zreg=griddata((x,y),z,(Xreg,Yreg),method='cubic')

for j in range(NyReg):
	for i in range(NxReg):
	        zreg2[i,j]=zreg[NxReg-i-1,j]
		
imgplot=plt.imshow(zreg2,extent=[xmin,xmax,ymin,ymax], aspect='auto')
imgplot.set_interpolation('bicubic')
# imgplot.set_cmap('RdBu')
imgplot.set_cmap('YlOrRd')
# imgplot.set_cmap('Blues')
# imgplot.set_cmap('Paired')
plt.colorbar()

## contourplot creation
CS=plt.contour(Xreg,Yreg,zreg,isolevels,linewdith=0.5,colors='k')

CS.levels=isolevels

## modification of contour ticks format
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
     def __repr__(self):
         str = '%.1f' % (self.__float__(),)
         if str[-1]=='0':
             return '%.0f' % self.__float__()
         else:
             return '%.1f' % self.__float__()

# Recast levels to new class
CS.levels = [nf(val) for val in CS.levels ]

#### Format control
fmt='%r'
# fmt = ticker.LogFormatterMathtext()
# fmt.create_dummy_axis()

# plt.clabel(CS, CS.levels, inline=1, fmt=fmt, fontsize=14, colors='black')

# plt.title(r'$Re( \varepsilon ) $')
# plt.title(r'$(H-H_{ref})/\Delta H $')
# plt.title('Potential (V)')
# plt.title('Temperature (K)')
plt.title('Normalized electric field $E_{int}/E_0$')
# plt.title('Liquid ratio')
# plt.title('$N_{e,h}$ ($m^{-3}$)')
# plt.title('Intensity ($W/m^{2}$)')

plt.xlabel("X ($\mu m$)")
plt.ylabel("Y ($\mu m$)")
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.show()

# plt.savefig('H-50fs-490mJ.eps',format='eps')


