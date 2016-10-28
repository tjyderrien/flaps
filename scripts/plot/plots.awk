#! /bin/gawk

BEGIN {
	CurrentFluence=0; 
}


{
	if(CurrentFluence==$1) {
		if($3>-3e-6 && $3<3e-6) print ;
	}
	else {
		print "\n" #saut de ligne si on change de fluence
	}
}
