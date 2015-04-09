#! gawk

BEGIN {
lambda=1030e-9
pi=3.141592654
Xvalue=-4.24661e-05
Xrange=-10e-6 
Trange=1e-12
# Xrange=2*lambda
# Yrange=Xrange #BonseSipe fourier space
Yrange=20e-6 #real space
}

{
####### CONSTRUCTION DE CARTES X Z
# On parcourt le fichier : quel que soit X, on parcourt Z
# on voudrait sauter une ligne lorsque X change pour faire du PM3D
# if($2!=Xvalue && $2<Xrange && $2>-Xrange && $3<Yrange) {
if($1 != Yvalue) 
{ print ""; }
# }
# else {
print $1, $2, $3;
# }

#definir la valeur de X
Yvalue=$1; 

}
