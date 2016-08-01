FC90=gfortran
FCFLAGS=-O2 -g -pg -Wall -ffree-line-length-none -I../libs/gmsh
#LDFLAGS=-L/opt/intel/mkl/lib/intel64 -R/opt/intel/mkl/lib/intel64 -shared-intel -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -openmp -lpthread -lm
#LDFLAGS=-L/usr/lib/lapack -llapack

F90_SRCS= grid.f90     \
          material.f90 \
          maths.f90    \
          mesh.f90     \
          mie.f90      \
          poisson.f90  \
          source.f90   \
          solver.f90   \
          types.f90   
          
         
OBJ90= $(F90_SRCS:.f90=.o)


all: main_explicit_Mie.o $(OBJ90)

%.o: %.f90
	@echo 'Building source file $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC90) $(FCFLAGS) -c -o $@ $<
	@echo 'Finished building: $<'
	@echo ' '

clean:
	rm -fr *.o *.mod

.PHONY: clean


main_explicit_Mie.o: main_explicit_Mie.f90 $(OBJ90)
mesh.o: mesh.f90 types.o
mie.o: mie.f90 maths.o
source.o: source.f90 maths.o types.o
solver.o: solver.f90 types.o
grid.o: grid.f90 maths.o
material.o: material.f90 types.o maths.o 
