CC=gcc
MPICC=mpicc
CFLAGS=-O3
LIBRARIES = -lm -lopenblas

V0:
	$(CC) $(CFLAGS) V0.c utilities.c -o V0 $(LIBRARIES)

V0_BLAS:
	$(CC) $(CFLAGS) V0_BLAS.c utilities.c -o V0_BLAS $(LIBRARIES)

V1:
	$(MPICC) $(CFLAGS) V1.c utilities.c -o V1 $(LIBRARIES)

V1_BLAS:
	$(MPICC) $(CFLAGS) V1_BLAS.c utilities.c -o V1_BLAS $(LIBRARIES)

V2_BLAS:
	$(MPICC) $(CFLAGS) V2_BLAS.c utilities.c -o V2_BLAS $(LIBRARIES)

all: V0 V0_BLAS V1 V1_BLAS V2_BLAS

.PHONY: clean

clean:
	rm -f V0 V0_BLAS V1 V1_BLAS V2_BLAS
