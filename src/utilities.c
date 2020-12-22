#include "utilities.h"

#define BILLION 1000000000L;

struct timespec start_time;
struct timespec stop_time;

double calculateExecutionTime()
{

    clock_gettime(CLOCK_MONOTONIC, &stop_time);

    double dSeconds = (stop_time.tv_sec - start_time.tv_sec);

    double dNanoSeconds = (double)(stop_time.tv_nsec - start_time.tv_nsec) / BILLION;

    return dSeconds + dNanoSeconds;
}

/*
    simple functions used for printing matrix contents
*/
void printMatrixInt(int* A, int n, int m) {
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			printf("%d ", A[i*m+j]);
		}
		printf("\n");
	}
}

void printMatrixDouble(double* A, int n, int m) {
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			printf("%lf ", A[i*m+j]);
		}
		printf("\n");
	}
}

void printMatrixMatlabFormat(double* A, int n, int m) {
  printf("A = [");
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      printf("%lf ", A[i*m+j]);
    }
    if(i < n -1)
      printf(";");
  }
  printf("]\n");
}

/* 
   =================================================================
   =================================================================
   This quicksort is used for k-select implementation
*/
void swapInt(int* a, int* b)  
{  
    int t = *a;  
    *a = *b;  
    *b = t;
}

void swapDouble(double* a, double* b) {
	double t = *a;
	*a = *b;
	*b = t;
}

int partition(double* D, int* indices, int start, int end)
{  
    double pivot = D[end]; // pivot  
    int i = (start - 1); // Index of smaller element  
    for (int j = start; j <= end - 1; j++)  
    {  
        // If current element is smaller than the pivot 
        if (D[j] < pivot)  
        {  
            i++; // increment index of smaller element  
            swapDouble(D+i, D+j);
            swapInt(indices+i, indices+j);
        }
    }  
    swapDouble(D+i+1, D+end);
    swapInt(indices+i+1, indices+end);
    return (i + 1);  
}


// Quick sort implementation to sort only k first elements
void quicksort(double* D, int* indices, int start, int end, int k)  
{  
    if (start < end)  
    {  
        // Find partition. After this command arr[p] is
        // in the right position
        int pi = partition(D, indices, start, end);
        // Separately sort elements before  
        // partition and after partition
        // check which to sort according to k
        
        // If k is smaller than partition then we don't need 
        // to sort elements biggest than p
        // else sort both parts (before and after p)
        if((k - 1) <= pi)  
        {
        	quicksort(D, indices, start, pi-1, k);  
        } else if((k - 1) > pi)
        {
        	quicksort(D, indices, start, pi - 1, k);
        	quicksort(D, indices, pi + 1, end, k);
    	}
    }
}

/* 
   =================================================================
   =================================================================
*/


/* 
   =================================================================
   =================================================================
   This is a simple distance calculator function
   used for calculating euclidean distance between 
   every X and Y point the old fashioned way
   X [m-by-d]
   Y [m-by-d]
   D [m-by-n]
*/
void findDistanceMatrix(double* X, double* Y, int n, int m, int d, double* D, int low, int high) {
  
  // Block size of query Y
  int size = high - low + 1;

  // e*e^T matrix [d-by-n]
  double* e1 = (double*)malloc(d*n*sizeof(double));
  // e^T*e matrix [size-by-d]
  double* e2 = (double*)malloc(size*d*sizeof(double));

  // matrices have all elements 1.0
  for(int i = 0; i < d*n; ++i)
    e1[i] = 1.0;

  for (int i = 0; i < size*d; ++i)
    e2[i] = 1.0;

  // Y.*Y product
  double* dotP_Y = (double*)malloc(size*d*sizeof(double));
  // X.*X product
  double* dotP_X = (double*)malloc(d*n*sizeof(double));

  // Calculate
  for(int i = low; i < low + size; ++i) {
    for(int j = 0; j < d; ++j) {
      dotP_Y[i*d+j] = Y[i*d+j] * Y[i*d+j];
    }
  }

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < d; ++j) {
      dotP_X[i*d+j] = X[i*d+j] * X[i*d+j];
    }
  }

  double* A = (double*)malloc(size*n*sizeof(double));
  for(int i = 0; i < n*size; ++i)
    A[i] = 0.0;

  printf("size=%d\tn=%d\td=%d\n", size, n, d);
  // A <- e2*dotP_X^T
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, size, n, d, 1.0, e2, size, dotP_X, d, 0.0, A, size);
  // A <- -2*X*Y^T + A
  // cblas_dgemm();
  // A <- dotP_Y*e2 + A
  // cblas_dgemm();

  printf("X = \n");
  printMatrixDouble(X, n, d);
  printf("X.*X = \n");
  printMatrixDouble(dotP_X, n, d);
  printf("e2 = \n");
  printMatrixDouble(e2, size, d);  
  printf("A = \n");
  printMatrixDouble(A, size, n);
  printf("Y = \n");
  printMatrixDouble(Y, size, d);
  printf("Y.*Y = \n");
  printMatrixDouble(dotP_Y, size, d);

  free(e1);
  free(e2);
  free(dotP_X);
  free(dotP_Y);
}

void findDMatrix(double* X, double* Y,int n, int m, int d, double* D) {
  // Compute D matrix using d^2 = (x1-x2)^2 + (y1-y2)^2 + ...
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
			D[i*n+j] = 0.0;
			for(int k = 0; k < d; ++k) {
				D[i*n+j] += (X[j*d+k] - Y[i*d+k])*(X[j*d+k] - Y[i*d+k]);
			}
			D[i*n+j] = sqrt(D[i*n+j]);
		}
	}
}