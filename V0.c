#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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
    int pivot = D[end]; // pivot  
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
void findDMatrix(double* X, double* Y,int n,int m,int d, double* D) {
	printf("\nD = \n");
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			D[i*n+j] = 0;
			for(int k = 0; k < d; ++k) {
				D[i*n+j] += (X[i*d+k] - Y[j*d+k])*(X[i*d+k] - Y[j*d+k]);
			}
			printf("%03d ", D[i*n+j]);
		}
		printf("\n");
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

knnresult kNN(double * X, double * Y, int n, int m, int d, int k) 
{
	int* D = (int*)malloc(n*m*sizeof(int));
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
	
	/*srand(time(NULL));
	if(argc < 3) {
		printf("Usage: ./V0 m n\n");
		return -1;
	}

	int m = atoi(argv[1]);
	int n = atoi(argv[2]);
	int d = 2;

	// int X2[4][2] = {{1, 2}, {-1, 3}, {0, 1}, {4, 1}};
	// int Y2[3][2] = {{3, 4}, {2, -1}, {4, 5}};

	int* X = (int*)malloc(m*d*sizeof(int));
	int* Y = (int*)malloc(n*d*sizeof(int));

	
	for(int i = 0; i < m; ++i) {
		for(int k = 0; k < d; ++k) {
			X[i*d+k] = rand() % 100;
		}
	}
	
	for(int i = 0; i < n; ++i) {
		for(int k = 0; k < d; ++k) {
			Y[i*d+k] = rand() % 10;
		}
	}
	
	printf("\nX = \n");
	for(int i = 0; i < m; ++i) {
		for(int k = 0; k < d; ++k) {
			printf("%d ", X[i*d+k]);
		}
		printf("\n");
	}

	printf("\nY = \n");
	for(int i = 0; i < n; ++i) {
		for(int k = 0; k < d; ++k) {
			printf("%d ", Y[i*d+k]);
		}
		printf("\n");
	}

	knnresult result = kNN(X, Y, n, m, d, 1);*/

	srand(time(NULL));


	int n = atoi(argv[1]);
	double* D = (double*)malloc(n*sizeof(double));
	if(D == NULL) {
		printf("Failed to allocate memory\n");
		return -1;
	}

	for(int i = 0; i < n; ++i) {
		D[i] = (double)rand() / RAND_MAX * atoi(argv[2]);
	}
	printMatrixDouble(D, 1, n);

	int* sortedIndexes = (int*)malloc(10*sizeof(int));
	if(sortedIndexes == NULL) {
		printf("Failed to allocate memory\n");
		return -1;
	}

	for(int i = 0; i < n; ++i)
		sortedIndexes[i] = i;

	printMatrixInt(sortedIndexes, 1, n);
	printf("\n");
	quicksort(D, sortedIndexes, 0, n-1);

	printMatrixDouble(D, 1, n);
	printMatrixInt(sortedIndexes, 1, n);
	printf("\n");
	for(int i = 0; i < n-1; ++i) {
		if(D[i] > D[i+1]) {
			printf("False\n");
			break;
		}
	}

	free(D);
	free(sortedIndexes);
}