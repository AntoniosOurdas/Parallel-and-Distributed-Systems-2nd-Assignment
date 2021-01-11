#include "utilities.h"

// Definition of the kNN result struct
typedef struct knnresult {
	int * nidx;		 //!< Indices (0-based) of nearest neighbors   [m-by-k]
	double * ndist;  //!< Distance of nearest neighbors 		   [m-by-k]
	int m;			 //!< Number of query points 				   [scalar]
	int k;			 //!< Number of nearest neighbors 			   [scalar]
} knnresult;


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
	knnresult result;

	// Allocate memory for distance matrix D
	double* D = (double*)malloc(m*n*sizeof(double));
	if(D == NULL) {
		printf("Couldn't allocate memory for D\n");
		return result;
	}

	// Allocate memory for indices matrix
	int* indices = (int*)malloc(n*sizeof(int));
	if(indices == NULL) {
		printf("Couldn't allocate memory for indices\n");
		return result;
	}

	// Calculate distance matrix D
	findDistanceMatrix(X, Y, n, m, d, D);

	result.nidx = (int*)malloc(m*k*sizeof(int));
	if(result.nidx == NULL) {
		printf("Couldn't allocate memory for nidx\n");
		return result;
	}

	result.ndist = (double*)malloc(m*k*sizeof(double));
	if(result.ndist == NULL) {
		printf("Couldn't allocate memory for ndist\n");
		return result;
	}

	result.m = m;
	result.k = k;

	// For every point in query set Y
	// find k nearest points from corpus set X
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			indices[j] = j;
		}
		// printf("\n=========================================\n");
		// printMatrixInt(indices, 1, n);
		// printMatrixDouble(D+i*n, 1, n);
		// printf("\nSorting D(%d,:):\n", i);

		quicksort(D+i*n, indices, 0, n-1, k);

		// printf("Done. Printing %d nearest neighbors:\n", k);
		// printMatrixInt(indices, 1, n);
		// printMatrixDouble(D+i*n, 1, n);
		// printf("\n=========================================\n");

		for(int l = 0; l < k; ++l) {
			result.ndist[i*k+l] = D[i*n+l];
			result.nidx[i*k+l] = indices[l];
		}
	}

	free(D);
	free(indices);
	return result;
}


// Main testing function
// Takes n, m, k as command line arguments
int main(int argc, char* argv[]) {

	srand(time(NULL));
	if(argc < 3) {
		printf("Usage: ./V0 m n d k\n");
		return -1;
	}

	// Set variables
	int n, m, d, k;

	if((n = atoi(argv[1])) <= 0) {
		printf("n must be positive integer\n");
		return -1;
	}

	if((m = atoi(argv[2])) <= 0) {
		printf("m must be positive integer\n");
		return -1;
	}

	if((d = atoi(argv[3])) <= 0) {
		printf("d must be positive integer\n");
		return -1;
	}

	if((k = atoi(argv[4])) <= 0) {
		printf("k must be positive integer\n");
		return -1;
	}

	if(k > n) {
		printf("You are asking for more neighbors than there are\n");
		return -1;
	}

	// Initialize corpus and query sets (X and Y respectively)
	double* X = (double*)malloc(n*d*sizeof(double));
	if(X == NULL) {
		printf("Couldn't allocate memory for X\n");
		return -1;
	}

	double* Y = (double*)malloc(m*d*sizeof(double));
	if(Y == NULL) {
		printf("Couldn't allocate memory for Y\n");
		return -1;
	}


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


	// printf("Hi %d %d %d %d\n", m, n, d, k);
	// printMatrixMatlabFormatDouble(X, n, d, "X");
	// printMatrixMatlabFormatDouble(Y, m, d, "Y");

	// printMatrixDouble(X, n, d);
	// printMatrixDouble(Y, m, d);
	// Find and print knnresult
	struct timespec start_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	knnresult kNNresult = kNN(X, Y, n, m, d, k);
	printf("\nV0 run time: %f\n", calculateExecutionTime(start_time));

	// printMatrixMatlabFormatDouble(X, n, d, "X");
	// printMatrixMatlabFormatDouble(Y, m, d, "Y");
	// printf("%d nearest neighbors distances:\n", k);
	// printMatrixDouble(kNNresult.ndist, m, k);
	// printf("\n%d nearest neighbors indexes:\n", k);
	// printMatrixInt(kNNresult.nidx, m, k);

	free(X);
	free(Y);
	free(kNNresult.nidx);
	free(kNNresult.ndist);

	return 0;
}
