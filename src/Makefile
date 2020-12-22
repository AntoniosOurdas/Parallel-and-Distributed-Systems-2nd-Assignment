CC=gcc
CFLAGS=-O3
LIBRARIES = -lm -lopenblas

V0:
	$(CC) $(CFLAGS) V0.c utilities.c -o V0 $(LIBRARIES)

all: V0

.PHONY: clean

clean:
	rm -f V0
