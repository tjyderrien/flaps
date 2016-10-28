#! gawk

# pour chaque t correspondant à ce qu'on veut sortir, on peut essayer de sortir les données dans fichier différent, ça permettra de ne lire le fichier qu'une seule fois au lieu de N

BEGIN {
	tau=40e-15
	dt0=1E-18; #dt of the code
	dt=1e-15; #dt between export
	tmax=1e-9
	tmin=-0.1900e-12
	initial=dt/dt0
	Each=3
}

#a chaque nouvelle valeur de temps, un nouveau fichier
{
	Condition1=(($1-tmin)/dt0)
	Condition2=(Each*dt/dt0)
#	print Condition1, Condition2
	Condition=int(Condition1 % Condition2)
	if(Condition==0 && $1 < tmax && $1>tmin) {
	print $0 >> int(($1-tmin)/dt)".vid"; 

	}
}
