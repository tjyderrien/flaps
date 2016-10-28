#! gawk

BEGIN {
previousX=0
previousY=0
}

{
	dx=$1-previousX
	dy=$2-previousY
#	derivative=dy/dx 
	previousX=$1
	previousY=$2
	
	if(dx!=0) { printf("%le\t%lg\n", $1, dy/dx); }
}
