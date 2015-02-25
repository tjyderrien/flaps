# Author: Thibault J.-Y. Derrien <thibault.derrien@gmail.com>

# Financial supports
[15/02/2012 - 14/12/2012] French National Research Agency (ANR, Agence Nationale de la Recherche), project "Ultrasonde".
[01/09/2015 - 31/08/2017] Marie Sklodowska Curie Actions, European Commision, project "QuantumLaP"

#Description
This code calculates the excitation, heating and transport of free-carriers (electrons and holes, separately) in silicon, coupling with phonons and their diffusive transport. 

#Features
* The shape of Silicon can be changed.
* Finite volume formulation has been used with a free-mesh formulation.
* Calculations are performed in the (x,y) plane, thus assuming the invariance by translation on Z axis. This hypothesis was required since laser does not originates from the tip symmetry axis, but from aside.
* The Fermi-Dirac statistics are taken into account in the calculations.