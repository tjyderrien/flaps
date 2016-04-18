CC=gfortran
LDFLAGS=-O2 -g -pg
SRC=$(wildcard *.f)
OBJ=$(SRC:.f=.o)

all: amos

amos: $(OBJ)

$(OBJ): %.o: %.f
	$(CC) -o $@ $^ $(LDFLAGS) -c

clean: 
	rm -fr *.o
