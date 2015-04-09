# Author: Thibault J.-Y. Derrien <thibault.derrien@gmail.com>

# Financial supports
1. [15/02/2012 - 14/12/2012] French National Research Agency (ANR, Agence Nationale de la Recherche), project "Ultrasonde".
1. [01/02/2013 - 31/01/2014] Adolf Martens Fellowship, German Ministry of Economy and Technology, Germany
1. [01/02/2014 - 30/04/2014] BAM Federal Institute for Material Research and Testing, Berlin, Germany
1. [01/09/2015 - 31/08/2017] Marie Sklodowska Curie Actions, European Commision, project "QuantumLaP"

#Description
This code calculates the excitation, heating and transport of free-carriers (electrons and holes, separately) of silicon irradiated by a laser pulse. The coupling with phonons and their diffusive transport is also calculated in time and space. Geometry is expressed in cartesian plane (x,y) and assumes invariance per translation. 

#Features
* The shape of the sample can be changed using a parametric function for the contour. 
* Equations are discretized using a *finite volume formulation* with mesh-free formulation.
* Calculations are performed in the (x,y) plane, thus assuming the invariance by translation on Z axis. This hypothesis was required since laser does not originates from the tip symmetry axis, but from aside.
* The Mie scattering theory for cylinders is locally applied for irradiations by the apex side. 
* The Fermi-Dirac statistics of electrons and holes are taken into account as a function of carrier density and temperature. 