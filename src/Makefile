CC=gfortran
LDFLAGS=-O2 -g -pg

F90_SRCS=mesh.f90 maths.f90 types.f90
OBJ90= $(F90_SRCS:.f90=.o)


all: $(OBJ90)

%.o: %.f90
	@echo 'Building source file $<'
	@echo 'Invoking: Fortran Compiler'
	$(CC) $(LDFLAGS) -c -o $@ $<
	@echo 'Finished building: $<'
	@echo ' '

clean:
	rm -fr *.o

.PHONY: clean

mesh.o: mesh.f90 types.o
