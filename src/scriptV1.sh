#!/bin/bash

if [ "$#" -lt 4 ] || [ "$#" -gt 4 ]; then
  echo "Usage: ./scriptV1 n d k np"
  exit 1
fi

# check if n is negative
if [ "$1" -lt 1 ]; then
  echo "n must be a positive integer"
  exit 1
fi

# check if d is negative
if [ "$2" -lt 1 ]; then
  echo "d must be a positive integer"
  exit 1
fi

# check if k is negative
if [ "$3" -lt 1 ]; then
  echo "k must be a positive integer"
  exit 1
fi

# check if k is greater than n
# (you can't search for more neighbors than there are)
if [ "$3" -gt "$1" ]; then
  echo "k can't be greater than n"
  exit 1
fi

# check if k is negative
if [ "$4" -lt 1 ]; then
  echo "np must be a positive integer"
  exit 1
fi

make clean
make all
mpirun -np "$4" V1 "$1" "$2" "$3"
