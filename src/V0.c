#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// #include "mpi.h"

// Definition of the kNN result struct
typedef struct knnresult {
	int * nidx;		 //!< Indices (0-based) of nearest neighbors   [m-by-k]
	double * ndist;  //!< Distance of nearest neighbors 		   [m-by-k]
	int m;			 //!< Number of query points 				   [scalar]
	int k;			 //!< Number of nearest neighbors 			   [scalar]
} knnresult;


// ================================================================= //
// ================================================================= //
// Functions used for printing matrices for testing 			    

void printMatrixInt(int* A, int m, int n) {
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			printf("%d ", A[i*m+j]);
		}
		printf("\n");
	}
}

void printMatrixDouble(double* A, int m, int n) {
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			printf("%lf ", A[i*m+j]);
		}
		printf("\n");
	}
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

int partition(double* D, int* sortedIndexes, int start, int end)
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
            swapInt(sortedIndexes+i, sortedIndexes+j);
        }
    }  
    swapDouble(D+i+1, D+end);
    swapInt(sortedIndexes+i+1, sortedIndexes+end);
    return (i + 1);  
}

void quicksort(double* D, int* sortedIndexes, int start, int end)  
{  
    if (start < end)  
    {  
        /* pi is partitioning index, arr[p] is now  
        at right place */
        int pi = partition(D, sortedIndexes, start, end);  
  
        // Separately sort elements before  
        // partition and after partition  
        quicksort(D, sortedIndexes, start, pi - 1);  
        quicksort(D, sortedIndexes, pi + 1, end);
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
   used for calculating euclidean distance the old fashioned way
*/
void findDMatrix(double* X, double* Y,int n, int m, int d, double* D) {
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			D[i*m+j] = 0;
			for(int k = 0; k < d; ++k) {
				D[i*n+j] += (X[i*d+k] - Y[j*d+k])*(X[i*d+k] - Y[j*d+k]);
			}
			D[i*m+j] = sqrt(D[i*m+j]);
		}
	}	
}

/* 
   =================================================================
   =================================================================
*/

//! Compute k nearest neighbors of each point in X [n-by-d]
/*!


	\param X Corpus data points 			[n-by-d]
	\param Y Query data points 				[m-by-d]
	\param n Number of corpus points		[scalar]
	\param m Number of query points 		[scalar]
	\param d Number of dimensions 			[scalar]
	\param k Number of neighbors 			[scalar]
	

	\return The kNN result
*/

knnresult kNN(double * X, double * Y, int n, int m, int d, int k, double* D) 
{
	findDMatrix(X, Y, n, m, d, D);
	knnresult result;
	int* sortedIndexes = NULL;
	// For each element j in corpus set Y find k nearest neighbors 
	// in query set X by finding k least distances of D(:,j)
	// This can be achieved by sorting D(:,j) 
	// and selecting first k elements
	for(int j = 0; j < n; ++j) {
		quicksort(D, sortedIndexes, j*n, (j+1)*n);
		for(int i = 0; i < k; ++i) {
			*(result.nidx+i) = sortedIndexes[i]; 
		}
	}

	return result;
}


// Main function
int main(int argc, char* argv[]) {
	
	srand(time(NULL));
	if(argc < 3) {
		printf("Usage: ./V0 m n k\n");
		return -1;
	}

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	int d = 2;
	int k = atoi(argv[3]);

	// int X2[4][2] = {{1, 2}, {-1, 3}, {0, 1}, {4, 1}};
	// int Y2[3][2] = {{3, 4}, {2, -1}, {4, 5}};

	double* X = (double*)malloc(n*d*sizeof(double));
	double* Y = (double*)malloc(m*d*sizeof(double));

	
	for(int i = 0; i < n; ++i) {
		for(int k = 0; k < d; ++k) {
			X[i*d+k] = -5.0 + (double)rand() / RAND_MAX * 10.0;
		}
	}
	
	for(int i = 0; i < m; ++i) {
		for(int k = 0; k < d; ++k) {
			Y[i*d+k] = -5.0 + (double)rand() / RAND_MAX * 10.0;
		}
	}
	
	printf("\nX = \n");
	printMatrixDouble(X, n, d);

	printf("\nY = \n");
	printMatrixDouble(Y, m, d);

	double* D = (double*)malloc(n*m*sizeof(double));
	if(D == NULL) {
		printf("Couldn't allocate memory for D\n");
		return -1;
	}

	int* sortedIndexes = (int*)malloc(m*sizeof(int));
	if(sortedIndexes == NULL) {
		printf("Couldn't allocate memory for sortedIndexes\n");
		return -1;
	}

	findDMatrix(X, Y, n, m, d, D);
	printf("\nD = \n");
	printMatrixDouble(D, n, m);
	printf("\nPrinting line by line:\n");
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
			sortedIndexes[j] = j;
		}
		printf("\n=========================================\n");
		printMatrixInt(sortedIndexes, 1, m);
		printMatrixDouble(D+i*m, 1, m);
		printf("\nSorting D(%d,:):\n", i);
		quicksort(D+i*m, sortedIndexes, 0, m-1);
		printf("Done. Printing %d nearest neighbors:\n", k);
		printMatrixInt(sortedIndexes, 1, k);
		printMatrixDouble(D+i*m, 1, k);
		printf("\n=========================================\n");
	}

	free(X);
	free(Y);
	free(D);
	free(sortedIndexes);
}
