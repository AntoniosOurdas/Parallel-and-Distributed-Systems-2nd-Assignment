#!/bin/bash
make clean
make all
mpirun -np "$4" V1 "$1" "$2" "$3"
