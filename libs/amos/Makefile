CC=gfortran
CFLAGS=-O2 -g -pg

all: amos

amos: $(OBJ)

$(OBJ): %.o: %.f
	$(CC) $(CFLAGS) -c $< -o $@
