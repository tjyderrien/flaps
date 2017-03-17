#!/usr/bin/perl
# ! @PERL@
# 
#

use Term::ANSIColor;

$ntests = 1;
for ($test=1;$test<=$ntests;$test++)
{
        $tol[$test] = 0.0001; #relative value
}

#$args=scalar(@ARGV);
#if($args==0)
#{
#   $starttest = 1;
#   $ntests = $number_of_tests;
#}


# Variable d'environnement
$omp="export OMP_NUM_THREADS=2; ";
# $EXE ="@PACKAGE_NAME@-@PACKAGE_VERSION@";
$EXE="flaps";
#Initialisation
for ($i=1 ; $i <= $ntests ; $i++) 
{
        $failed[$i]=1.0;
        $message[$i]="Test $i:";
}

# set the message array (default test1 test2 .. testN)
$message[1] ="Test  1:  [Analytic] Zero energy integration ..................";
$message[2] ="Test  2:  [Analytic] Excitation rate with homogeneous source ..";
$message[3] ="Test  3:  [Analytic] Heating with homogeneous source ..........";
$message[4] ="Test  4:  [Analytic] Heating via electron-lattice heat transfer";
$message[5] ="Test  5:  [Analytic] Auger recombination rate .................";
$message[6] ="Test  6:  [Analytic] Heating induced by Auger recombination ...";
$message[7] ="Test  7:  [Regression] Mie scattering TM total energy .........";
$message[8] ="Test  8:  [Regression] Mie scattering TE total energy .........";
$message[9] ="Test  9:  [Regression] Mie induced electron heating ...........";
$message[10]="Test 10:  [Regression] Mie induced lattice heating ............";
$message[11]="Test 11:  [Analytic] Fermi-Dirac integral convergence .........";
$message[12]="Test 12:  [Regression] Fermi-Dirac regressive test ............";

#Commandes pour les tests
for ($i=1 ; $i <= $ntests ; $i++) 
{
	#TODO: copy test.$i.in flaps.in
	$testscript[$i]="$omp cp Test".$i.".in flaps.in; ../src/".$EXE." >log";
}

# clean the directory from previous failed tests
system("rm -f *faile*");
system("clear");

print color "cyan";
print "*****************************\n";
print "Performing Unit-Tests        \n";
print "*****************************\n\n";
print color "reset";

for ($test = 1 ; $test <= $ntests ; $test++) {
        print $message[$test];
        system($testscript[$test]);
        $ODT="Test".$test."/TimeMax.dat"; #ODT = reference !
	$file="output/TimeMax.dat";
        @old = `cat $ODT | grep -v \"#\"`;
        @new = `cat $file | grep -v \"#\"`;
	#on utilise wc pour avoir accès au nombre de lignes et de colonnes
        @size = split(/ +/,`cat $file | grep -v \"#\" | wc`);
        $lignes = @size[1];
        $colonnes = @size[2]/$lignes;
        for ($i=0 ; $i < $lignes ; $i++)
	{
            @val_old = split(/ +/,$old[$i]);
            @val_new = split(/ +/,$new[$i]);
            chomp($val_old[$colonnes]);
            chomp($val_new[$colonnes]);
            for ($j=1 ; $j <= $colonnes ; $j++)
	    {
                if ($val_new[$j]=~ m/NaN/ or $val_new[$j]=~ m/NAN/)
		{
			goto failed;
		}
                if (abs($val_old[$j]) > 0.1)
		{
                    $diff = abs(($val_old[$j] - $val_new[$j])/$val_old[$j]) ;
                }
		else
		{
                    $diff = abs($val_old[$j] - $val_new[$j]) ;
                }
                if (not ($diff < $tol[$test]) )
		{
			goto failed;
                }
            }
        }

	#Le test est passé
        system("rm -rf output");
	system("rm -f log");
        print color "bold green";
        print "OK";
        print color "reset";
        print "\n";
	$failed[$test]=0.0;
        goto nexttest;
        failed: #Le test n'est pas passé
        print color "red";
        print "Failed";
        print color "reset";
        print "\n";
        $mes="mv output output".$test.".failed";
        system($mes);
        $mes="mv log Test".$test.".failed.log";
        system($mes);
        nexttest:
}

$sum=0.;
for ($i=$starttest;$i<=$ntests;$i++)
{
        $sum=$sum+$failed[$i];
}

if ($sum == 0)
{
        print "\n";
        print "*************** Unit Tests OK  ***************\n";
        print "***********  You can now use Flaps  **********\n\n";
}
else
{
        print "\n";
        print "Beware: some tests didn't work. Maybe it is only matter of\n";
        print "machine's precision. We suggest to graphically compare the\n";
        print "following files: \n\n";
        for ($i=$starttest;$i<=$ntests;$i++)
	{
        	if ($failed[$i] > 0.)
		{
                      print "\t 'Result".$i.".failed' and 'ODT/Test".$i."/Result'\n";
                }
        }
        print "\n which are in ./tests directory\n\n";
}

end:



