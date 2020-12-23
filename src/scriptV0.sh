#!/bin/bash
make clean
make all
./V0 "$1" "$2" "$3" "$4"
