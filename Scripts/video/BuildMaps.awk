#! gawk

# pour chaque t correspondant à ce qu'on veut sortir, on peut essayer de sortir les données dans fichier différent, ça permettra de ne lire le fichier qu'une seule fois au lieu de N

BEGIN {
	tau=40e-15
	dt0=100E-18; #dt of the code
	dt=1e-12; #dt between export
	tmax=0.47200E-11 #100*tau #0.2e-12
	tmin=-0.19000E-12 #-0.24990E-11 #-4.9978000000e-13 #0.159E-12 #min of the one exported!
	M=2001; N=151
	NbPoints=100
	NbPointsDec=NbPoints*1e0
#	print "NbPointsDec", 1e0/NbPointsDec
	
#	print ((tmax-tmin))/(dt)^(1/)
	Reason=((tmax-tmin)/(dt))**(1e0/NbPointsDec)
	print "Reason", Reason
	
	Outputs[0]=dt
	for(i=1; i<=NbPoints; i++) {
		Outputs[i]=Outputs[i-1]*Reason
	}
	for(i=0; i<=NbPoints; i++) {
		Outputs[i]=Outputs[i]-tmin-dt
	}
#	Outputs[0]=-4e0*tau; Outputs[1]=-3e0*tau; Outputs[2]=-2e0*tau; Outputs[3]=-tau; Outputs[4]=0e0; Outputs[5]=tau; Outputs[6]=2e0*tau; Outputs[7]=3e0*tau; Outputs[8]=4e0*tau; Outputs[9]=10e0*tau; Outputs[10]=20e0*tau; Outputs[11]=40e0*tau; Outputs[12]=50e0*tau; Outputs[13]=100e0*tau; Outputs[14]=200e0*tau; Outputs[15]=400e0*tau; Outputs[16]=500e0*tau; Outputs[17]=1000e0*tau; 	Outputs[18]=2000e0*tau; Outputs[19]=5000e0*tau; Outputs[20]=1e4*tau; Outputs[21]=2e4*tau; 	Outputs[22]=5e4*tau; Outputs[23]=1e5*tau; Outputs[24]=2e5*tau; Outputs[25]=5e5*tau; Outputs[26]=1e6*tau;


	indice=0
	PreviousTime=tmin
	Exportedlines=0
	MapSize=M*N
}

#a chaque nouvelle valeur de temps, un nouveau fichier
{
#	Condition1=int(($1-tmin)/dt0)
#	Condition2=int(Each*dt/dt0)
#	Test=Condition1/Condition2
#	Condition=int(Condition1 % Condition2)
#	print Condition1, Condition2, Condition, Test
	if($1 < tmax && $1>tmin) {
#enable me	if($1 >= Outputs[indice]) {
			print $0 >> ($1-tmin)/dt".vid";
			Exportedlines=Exportedlines+1
#			print $1 " exported, id=", indice
			if(Exportedlines>=MapSize) {
				Exportedlines=0
				indice=indice+1
			}
#enable me		}
	}
	
}
