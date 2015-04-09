#!/usr/bin/env python
import pylab
import matplotlib.pyplot as plt
import numpy as np
# from numpy import * 

filename='LaplaceMatrix.dat'
data=np.loadtxt('LaplaceMatrix.dat')
Z=data[:,:]

mesh=np.loadtxt('meshVessel.dat')
x,y=mesh[:,1], mesh[:,2]


# potential
p=plt.imshow(Z)
fig=plt.gcf()
plt.clim()
plt.title("Potential map")
#plt.pause(3)

Fields=np.loadtxt('SliceDepthVessel.dat', unpack=True)

scale=10000.

Ex,Ey=Fields[:,5]/scale, Fields[:,6]/scale
plt.figure()

Q=plt.quiver(x,y,Ex,Ey)
#qk = plt.quiverkey(Q, 0.9, 0.95, 2, r'$2 \frac{m}{s}$',
               #labelpos='E',
               #coordinates='figure',
               #fontproperties={'weight': 'bold'})
#plt.title('Field in the needle')

plt.pause(0)
