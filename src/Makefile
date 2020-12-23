CC=gcc
MPICC=mpicc
CFLAGS=-O3
LIBRARIES = -lm -lopenblas

V0:
	$(CC) $(CFLAGS) V0.c utilities.c -o V0 $(LIBRARIES)

V1:
	$(MPICC) V1.c utilities.c -o V1 $(LIBRARIES)

all: V0 V1

.PHONY: clean

clean:
	rm -f V0 V1
