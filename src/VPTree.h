#ifndef VPTREE_H_
#define VPTREE_H_
#include "utilities.h"

// Used for printing tree
#define COUNT 10

typedef struct VPT {
  double* VP;       // Vantage point
  double median;    // Median distance
  struct VPT* in;   // Inner tree
  struct VPT* out;  // Outer tree
} VPT;

//Auxiliary functions to print in tree like format
void print2DUtil(VPT* root, int space);
void print2D(VPT* root);

// Delete VPT tree recursively
void delete(VPT* T);

// Partition function used by quick select
int partitionQS(double* arr, int l, int r);

// k-select to find k smallest element of arr
double kthSmallest(double* arr, int l, int r, int k);

// Recursively make VPT tree
void makeVPT(VPT* T, double* X, int n, int d);

void searchVPT(VPT *T, double *X, int d, double* ndist, int* nidx);

#endif
