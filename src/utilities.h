#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cblas.h>

#ifndef UTILITIES_H_
#define UTILITIES_H_

double calculateExecutionTime(struct timespec start_time);

void printMatrixInt(int* A, int n, int m);

void printMatrixDouble(double* A, int n, int m);

void swapInt(int* a, int* b);

void swapDouble(double* a, double* b);

int min(int a, int b);

int partition(double* D, int* sortedIndexes, int start, int end);

void quicksort(double* D, int* sortedIndexes, int start, int end, int k);

void findDistanceMatrix(double* X, double* Y, int n, int m, int d, double* D);

void findDMatrix(double* X, double* Y,int n, int m, int d, double* D);

#endif
