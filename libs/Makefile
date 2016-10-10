CC=gfortran #ifort
LDFLAGS=-O2 -g -pg

SRC=Bivariate.f zeroin.f
OBJ=Bivariate.o zeroin.o

all: libsf libamos libgmsh

libsf: $(OBJ)

$(OBJ): %.o: %.f
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
