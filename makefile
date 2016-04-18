CC=gfortran
CFLAGS=-I/usr/include -I/opt/intel/mkl/include
LDFLAGS=-L/opt/intel/mkl/lib/intel64 -R/opt/intel/mkl/lib/intel64 -shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -openmp -lpthread -lm -o

flaps: src/main_explicit_Mie.o libs/Bivariate.o libs/zeroin.o libs/amos/*.o libs/gmsh/libmsh2vf.o libs/control_file.o 
	$(CC) -o Flaps_explicit.out
