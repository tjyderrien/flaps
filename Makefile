DebugLine=-g -pg -fbounds-check -w -fbacktrace -Wno-unused-variable -Wall
Optimization=$(DebugLine) -O2 -ffree-line-length-none -ffixed-line-length-none -finteger-4-integer-8 -mcmodel=small -I/usr/include -llapack -fopenmp -lm

CC=gfortran
CFLAGS=$(Optimization)
EXEC=Flaps_explicit.out

#LDFLAGS=-L/opt/intel/mkl/lib/intel64 -R/opt/intel/mkl/lib/intel64 -shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -openmp -lpthread -lm

OBJ_LIBS=libs/Bivariate.o libs/zeroin.o $( wildcards libs/amos/*.o) libs/gmsh/libmsh2vf.o libs/control_file.o

all: linking main_explicit_Mie.o

linking: main_explicit_Mie.o $(OBJ_LIBS)
	$(CC) $(CFLAGS) $(OBJ_LIBS) -o $(EXEC)

main_explicit_Mie.o: $(OBJ_LIBS)
	$(CC) -c $(CFLAGS) -I./libs/gmsh src/main_explicit_Mie.f90

$(OBJ_LIBS): 
	cd libs && $(MAKE)

clean: 
	rm -fr *.o
	cd src && rm -fr *.o 
	cd libs && $(MAKE) clean
