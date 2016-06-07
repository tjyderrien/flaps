FC90=gfortran
FCFLAGS=-O2 -g -pg -Wall -ffree-line-length-none -I../libs/gmsh

F90_SRCS=mesh.f90 maths.f90 types.f90 mie.f90
OBJ90= $(F90_SRCS:.f90=.o)


all: main_explicit_Mie.o $(OBJ90)

%.o: %.f90
	@echo 'Building source file $<'
	@echo 'Invoking: Fortran Compiler'
	$(FC90) $(FCFLAGS) -c -o $@ $<
	@echo 'Finished building: $<'
	@echo ' '

clean:
	rm -fr *.o

.PHONY: clean


main_explicit_Mie.o: main_explicit_Mie.f90 $(OBJ90)
mesh.o: mesh.f90 types.o
mie.o: mie.f90 maths.o
