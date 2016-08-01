CC=ifort #gfortran
LDFLAGS=-O3 #-g -pg
SRC=$(wildcard *.f)
OBJ=$(SRC:.f=.o)

all: amos

amos: $(OBJ)

$(OBJ): %.o: %.f
	$(CC) -o $@ $^ $(LDFLAGS) -c

clean: 
	rm -fr *.o
