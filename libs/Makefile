CC=gfortran
LDFLAGS=-O2 -g -pg

SRC=Bivariate.f zeroin.f
SRC90=control_file.f90
OBJ=Bivariate.o zeroin.o
OBJ90=control_file.o

all: libsf libsf90 libamos libgmsh

libsf: $(OBJ)

$(OBJ): %.o: %.f
	$(CC) -o $@ $^ $(LDFLAGS) -c

libsf90: $(OBJ90)

$(OBJ90): %.o: %.f90
	$(CC) -o $@ $^ $(LDFLAGS) -c

libamos:
	cd amos && $(MAKE)

libgmsh:
	cd gmsh && $(MAKE)

clean:
	rm -fr *.o
	cd amos && $(MAKE) clean
	cd gmsh && $(MAKE) clean

.PHONY: clean
