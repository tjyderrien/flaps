#DebugLine=-g -check all -warn all -backtrace 
DebugLine= -pg -g -fbounds-check -w -fbacktrace -Wno-unused-variable -Wall
Optimization=$(DebugLine) -O3 -I/usr/include -llapack -lm  -fopenmp #-mcmodel=large

#-finteger-4-integer-8

#CC=ifort
CC=gfortran
CFLAGS=$(Optimization)
EXEC=Flaps_explicit.out

#LDFLAGS=-L/opt/intel/mkl/lib/intel64 -R/opt/intel/mkl/lib/intel64 -shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -openmp -lpthread -lm
LDFLAGS=-L/usr/lib/lapack -llapack

OBJ_LIBS=libs/Bivariate.o libs/zeroin.o $(wildcard libs/amos/*.o) libs/gmsh/libmsh2vf.o 

OBJS = $(wildcard src/*.o)

 
all: $(EXEC)

$(EXEC): external_libs src_files
	@echo 'Building target: $@'
	$(CC) $(CFLAGS) -o $@  $(OBJ_LIBS) $(OBJS) $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo ' '
	@echo '************  Compilation OK  ************';
	@echo '******  You can now do >make test< *******';

src_files: external_libs 
	cd src && $(MAKE)


external_libs:
	@echo 'Building libraries' 
	cd libs && $(MAKE)
	@echo 'Finished building libraries'
	@echo ' '
 

clean: 
	rm -fr *.o *.mod
	cd src && $(MAKE) clean
	cd libs && $(MAKE) clean

