CC=gfortran
CFLAGS=-I

flaps: main_explicit_Mie.o Bivariate.o zeroin.o amos/*.o 
	$(CC) -o Flaps_explicit.out
