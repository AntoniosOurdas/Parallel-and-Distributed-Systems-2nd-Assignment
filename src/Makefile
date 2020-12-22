CC=gcc
CFLAGS=-I Docker-main/pds-codebase/2nd_Assignment -O3
LIBRARIES = -lm -lopenblas

V0:
	$(CC) $(CFLAGS) Docker-main/pds-codebase/2nd_Assignment/V0.c Docker-main/pds-codebase/2nd_Assignment/utilities.c -o Docker-main/pds-codebase/2nd_Assignment/V0 $(LIBRARIES)

all: V0

.PHONY: clean

clean:
	rm -f V0
