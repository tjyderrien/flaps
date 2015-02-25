# Author: Thibault J.-Y. Derrien <thibault.derrien@gmail.com>
# Sponsored by french "Agence Nationale de la Recherche" under the project name "Ultrasonde" from 15-02-2012 - 14-12-2012.
# Now under development supported by Marie Curie Actions, project "QuantumLaP", from 1-09-2015 - 31-08-2017

This code calculates the excitation, heating and transport of free-carriers (electrons and holes, separately) in silicon, coupling with phonons and their diffusive transport. 

Features:
- The shape of Silicon can be changed.
- Finite volume formulation has been used with a free-mesh formulation.
- Calculations are performed in the (x,y) plane, thus assuming the invariance by translation on Z axis. This hypothesis was required since laser does not originates from the tip symmetry axis, but from aside.
- The Fermi-Dirac statistics have been introduced in the calculations.