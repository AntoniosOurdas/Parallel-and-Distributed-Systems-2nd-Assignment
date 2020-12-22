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

//! Compute distributed all-kNN of points in X
/*!


	\param X Corpus data points 			[n-by-d]
	\param n Number of corpus points		[scalar]
	\param d Number of dimensions 			[scalar]
	\param k Number of neighbors 			[scalar]
	

	\return The kNN result
*/

knnresult distrAllkNN(double* X, int n, int d, int k) 
{
	knnresult result;
	
	int  numtasks, rank, len, rc; 
   	char hostname[MPI_MAX_PROCESSOR_NAME];
	// initialize MPI  
   	MPI_Init(&argc,&argv);

   	// get number of tasks 
   	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

   	// get my rank  
   	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   	// this one is obvious  
   	MPI_Get_processor_name(hostname, &len);
   	printf ("Number of tasks= %d My rank= %d Running on %s Length %d\n", numtasks,rank,hostname,len);

        // do some work with message passing 

	// Allocate memory for distance matrix D
	double* D = (double*)malloc(n*n*sizeof(double));
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
	findDMatrix(X, X, n, n, d, D);

	result.nidx = (int*)malloc(n*k*sizeof(int));
	if(result.nidx == NULL) {
		printf("Couldn't allocate memory for nidx\n");
		return result;
	}

	result.ndist = (double*)malloc(n*k*sizeof(double));
	if(result.ndist == NULL) {
		printf("Couldn't allocate memory for ndist\n");
		return result;
	}

	result.m = n;
	result.k = k;

	// For every point in query set Y
	// find k nearest points from corpus set X
	for(int i = 0; i < n; ++i) {
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
	// done with MPI  
   	MPI_Finalize();
}


// Main testing function
// Takes n, m, k as command line arguments
int main(int argc, char* argv[]) {
	
	srand(time(NULL));
	if(argc < 3) {
		printf("Usage: ./V1 n d k\n");
		return -1;
	}

	// Set variables
	int n, d, k;

	if((n = atoi(argv[1])) <= 0) {
		printf("n must be positive integer\n");
		return -1;
	}
	
	if((d = atoi(argv[2])) <= 0) {
		printf("d must be positive integer\n");
		return -1;
	}

	if((k = atoi(argv[3])) <= 0) {
		printf("k must be positive integer\n");
		return -1;
	}
	
	if(k > n) {
		printf("You are asking for more neighbors than there are\n");
		return -1;
	}

	// Initialize X	
	double* X = (double*)malloc(n*d*sizeof(double));
	if(X == NULL) {
		printf("Couldn't allocate memory for X\n");
		return -1;
	}

	for(int i = 0; i < n; ++i) {
		for(int k = 0; k < d; ++k) {
			X[i*d+k] = -5.0 + (double)rand() / RAND_MAX * 10.0;
		}
	}

	printMatrixMatlabFormat(X, n, d);
	// Find and print knnresult
	knnresult kNNresult = distrAllkNN(X, n, d, k);

	printf("%d nearest neighbors distances:\n", k);
	printMatrixDouble(kNNresult.ndist, n, k);
	printf("\n%d nearest neighbors indexes:\n", k);
	printMatrixInt(kNNresult.nidx, n, k);

	free(X);
	free(kNNresult.nidx);
	free(kNNresult.ndist);

	return 0;
}