#! gawk

# script pour extraire les infos interessantes des fichiers passes en argument
BEGIN {
nom=FILENAME
}
{
# print FILENAME
	if($3==0) { print $0 >> FILENAME".tmp" }
}
