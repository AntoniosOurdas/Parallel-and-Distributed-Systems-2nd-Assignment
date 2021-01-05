#include "utilities.h"
#include "mpi.h"
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

void kNN(double * X, double * Y, int n, int m, int d, int k, knnresult* result, int step)
{

	knnresult tempResult;

	tempResult.nidx = (int*)malloc(m*k*sizeof(int));
	if(tempResult.nidx == NULL) {
		printf("Couldn't allocate memory for nidx\n");
		return;
	}

	tempResult.ndist = (double*)malloc(m*k*sizeof(double));
	if(tempResult.ndist == NULL) {
		printf("Couldn't allocate memory for ndist\n");
		return;
	}

	tempResult.m = m;
	tempResult.k = k;


	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int numtasks = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

	// Allocate memory for distance matrix D
	double* D = (double*)malloc(m*n*sizeof(double));
	if(D == NULL) {
		printf("Couldn't allocate memory for D\n");
		return;
	}

	// Allocate memory for indices matrix
	int* indices = (int*)malloc(n*sizeof(int));
	if(indices == NULL) {
		printf("Couldn't allocate memory for indices\n");
		return;
	}

	// Calculate distance matrix D
	findDMatrix(X, Y, n, m, d, D);

	// For every point in query set Y
	// find k nearest points from corpus set X
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			if((rank - step) >= 0) {
				indices[j] = j + (rank - step)*n;
			} else {
				indices[j] = j + ((numtasks - 2) - (rank - step))*n;
			}
		}

		quicksort(D+i*n, indices, 0, n-1, k);

		for(int l = 0; l < k; ++l) {
			tempResult.ndist[i*k+l] = D[i*n+l];
			tempResult.nidx[i*k+l] = indices[l];
		}

	}

	// Update k nearest neighbors
	// by merging two structs
	// until k shortest distances
	// (i.e. nearest points)
	// are found

	// Temp variables to store data to be merged
	int* tempIndex = (int*)malloc(k*sizeof(int));
	double* tempDist = (double*)malloc(k*sizeof(double));

	// Indices used for merging
	int index1 = 0;
	int index2 = 0;
	int count = 0;

	for(int i = 0; i < m; ++i) {

		index1 = 0;
		index2 = 0;
		count = 0;

		// Copy results data in temp variables
		for(int j = 0; j < k; ++j) {
			tempIndex[j] = result->nidx[i*k+j];
			tempDist[j] = result->ndist[i*k+j];
		}

		// Merge
		while((index1 < k) && (index2 < k) && (count < k)) {

			if(tempResult.ndist[i*k+index1] < tempDist[index2]) {

				result->nidx[i*k+count] = tempResult.nidx[i*k+index1];
				result->ndist[i*k+count] = tempResult.ndist[i*k+index1];
				++index1;

			} else {

				result->nidx[i*k+count] = tempIndex[index2];
				result->ndist[i*k+count] = tempDist[index2];
				++index2;

			}
			++count;

		}

		// if(index1 == k) {
		// 	for(int l = index2; l < k; ++l) {
		// 		result->nidx[i*k+count] = tempIndex[l];
		// 		result->ndist[i*k+count] = tempDist[l];
		// 	}
		// } else if(index2 == k) {
		// 	for(int l = index1; l < k; ++l) {
		// 		result->nidx[i*k+count] = tempResult.nidx[i*k+index1];
		// 		result->ndist[i*k+count] = tempResult.ndist[i*k+index1];
		// 	}
		// }

	}


	free(tempDist);
	free(tempIndex);

	free(tempResult.nidx);
	free(tempResult.ndist);

	free(D);
	free(indices);

	return;

}

knnresult distrAllkNN(double* X, int n, int d, int k)
{

	// define knn struct
	// each proccess has its own struct
	// afterawds proccess 0
	// will merge all structs together
	knnresult* result = (knnresult*)malloc(1*sizeof(knnresult));

	// Variables used for MPI communication
	int  numtasks, rank, len, rc;
  char hostname[MPI_MAX_PROCESSOR_NAME];

	MPI_Status status;
	MPI_Comm intercomm;
	// initialize MPI
  // MPI_Init(NULL, NULL);

 	// MPI_Comm_spawn("worker_program", MPI_ARGV_NULL, n_spawns, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
  // get number of
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	// get my rank
 	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 	// this one is obvious
 	MPI_Get_processor_name(hostname, &len);
 	// printf ("Number of tasks = %d My rank = %d Running on %s Length %d\n", numtasks, rank, hostname, len);

	// We use two requests for each task
	// one for sending data
	// and one for receiving data
	// both asychronous
	MPI_Request request[2*numtasks];

	// define previous and next task
	int prev = 0, next = 0;
	if((prev = rank - 1) == -1)
		prev = numtasks - 1;
	if((next = rank + 1) == numtasks)
		next = 0;

	int chunk_size_est = n / numtasks;
	// printf("Chunck size: %d\n", chunk_size);

	// Local data pointers
	double* Xlocal = NULL;
	double* Y = NULL;
	double* Z = NULL;

	// Here we store chunk size of each proccess so that
	// we know how much data to send to other proccesses
	int* chunk_size = (int*)malloc(numtasks*sizeof(int));
	for(int i = 0; i < numtasks; ++i) {
		// If n can't be evenly divided in all proccesses
		// then assign each of remaining (n % numtasks) to first (n % numtasks) tasks
		chunk_size[i] = ( ( (n % numtasks) == 0 ) ? chunk_size_est : ( ( i <= (n % numtasks - 1) ) ? (chunk_size_est + 1) : (chunk_size_est)));
		// printf("Proccess %d chunck_size: %d\n", i, chunk_size[i]);
	}

	int max_chunk_size = chunk_size[0];

	result->nidx = (int*)malloc(chunk_size[rank]*k*sizeof(int));
	if(result->nidx == NULL) {
		printf("Couldn't allocate memory for nidx\n");
		return;
	}

	result->ndist = (double*)malloc(chunk_size[rank]*k*sizeof(double));
	if(result->ndist == NULL) {
		printf("Couldn't allocate memory for ndist\n");
		return;
	}

	result->m = chunk_size[rank];
	result->k = k;

	// ====================================================================== //
	// Part 1
	// Split X evenly in multiple proccesses
	// ====================================================================== //

	// Work to be done for mother proccess
	// Note: proccess 0 is called mother proccesses
	// because it is the one who creates X
	// and then distributes it among other proccesses
	// while keeping a part for herself
	if(rank == 0) {

		// Allocate space for local data
		Xlocal = (double*)malloc(max_chunk_size*d*sizeof(double*));
		Y = (double*)malloc(max_chunk_size*d*sizeof(double*));
		Z = (double*)malloc(max_chunk_size*d*sizeof(double*));

		// Copy first chunk of X in local storage
		for(int i = 0; i < chunk_size[rank]; ++i) {
			for(int j = 0; j < d; ++j) {
				Xlocal[i*d+j] = X[i*d+j];
			}
		}

		// Then send each other chunk to corresponding proccess
		for(int i = 0; i < numtasks - 1; ++i) {
			// Send each i-th chunk of X to
			// ith proccess
			// we keep the first chunk i=0
			// printf("Sending X from %d to %d\n", chunk_size[rank]+chunk_size[i+1]*i, i+1);
			MPI_Isend((X+chunk_size[0]*d+chunk_size[i+1]*d*i), chunk_size[i+1]*d, MPI_DOUBLE, i+1, 99, MPI_COMM_WORLD, &request[2*(i+1)]);
		}

	// Work to be done for other proccesses

	}
	else if (rank != 0) {

		// Allocate space for local data
		Xlocal = (double*)malloc(max_chunk_size*d*sizeof(double*));
		Y = (double*)malloc(max_chunk_size*d*sizeof(double*));
		Z = (double*)malloc(max_chunk_size*d*sizeof(double*));

		// printf("I am daughter and my chunk size is %d\n", chunk_size[rank]);

		MPI_Recv(Xlocal, chunk_size[rank]*d, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);

	}
	// ====================================================================== //
	// End of Part 1
	// ====================================================================== //

	for(int i = 0; i < chunk_size[rank]; ++i) {
		for(int j = 0; j < k; ++j) {
			result->nidx[i*k+j] = 0;
			result->ndist[i*k+j] = 100.0;
		}
	}

	// ====================================================================== //
	// Part 2
	// Start finding k nearest neighbors for local corpus set
	// and then send exchange data in a ring topology
	// i.e. send to next and receive from previous
	// ====================================================================== //
	for(int i = 0; i < chunk_size[rank]; ++i) {
		for(int j = 0; j < d; ++j) {
			Y[i*d+j] = Xlocal[i*d+j];
		}
	}

	// Start calculating while receiving from others
	for(int i = 0; i < numtasks; ++i) {

		// While proccessing receive from others if they have finished before
		MPI_Irecv(Z, max_chunk_size*d, MPI_DOUBLE, prev, 99, MPI_COMM_WORLD, &request[2*rank+1]);

		// Calculate
		kNN(Y, Xlocal, chunk_size[rank], chunk_size[rank], d, k, result, i);

		// Send data when finished
		MPI_Isend(Xlocal, max_chunk_size*d, MPI_DOUBLE, next, 99, MPI_COMM_WORLD, &request[2*rank]);
		// Wait to receive data
		MPI_Wait(&request[2*rank+1], &status);

		// Move Z to Y to be proccessed
		for(int i = 0; i < max_chunk_size; ++i) {
			for(int j = 0; j < d; ++j) {
				Y[i*d+j] = Z[i*d+j];
			}
		}

		// Wait to send data before proccessing again
		MPI_Wait(&request[2*rank], &status);

	}
	// ====================================================================== //
	// End of Part 2
	// ====================================================================== //



	// ====================================================================== //
	// Part 3
	// At this point mother proccess 0
	// will merge all knnresult structs in one
	// and return it to main
	// ====================================================================== //
	if(rank == 0) {

		knnresult final;
		final.nidx = (int*)malloc(n*k*sizeof(int));
		final.ndist = (double*)malloc(n*k*sizeof(double));
		final.m = n;
		final.k = k;

		for(int i = 0; i < chunk_size[rank]; ++i) {
			for(int j = 0; j < k; ++j) {
					final.nidx[i*k+j] = result->nidx[i*k+j];
					final.ndist[i*k+j] = result->ndist[i*k+j];
			}
		}
		// Receive structs from all other procceses
		// and save them into final struct to be returned
		for(int i = 1; i < numtasks; ++i) {

			MPI_Recv(final.nidx + chunk_size[0]*k + (i-1)*chunk_size[i]*k, chunk_size[i]*k,
			MPI_INT, i, 99, MPI_COMM_WORLD, &status);

			MPI_Recv(final.ndist + chunk_size[0]*k + (i-1)*chunk_size[i]*k, chunk_size[i]*k,
			MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);

		}

		return final;

	} else {

		// Send knnresult in mother
		MPI_Send(result->nidx, max_chunk_size*k, MPI_INT, 0, 99, MPI_COMM_WORLD);
		MPI_Send(result->ndist, max_chunk_size*k, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
		free(result->nidx);
		free(result->ndist);
		free(result);
		return *result;

	}

	// ====================================================================== //
	// End of Part 3
	// ====================================================================== //
}


// Main driver function
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

	int rank = 0;
	// int rank = 0;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double *X = NULL;

	if(rank == 0)
	{

		// printf("I am mother. Creating X array...\n");
		// Initialize X
		X = (double*)malloc(n*d*sizeof(double));
		if(X == NULL) {
			printf("Couldn't allocate memory for X\n");
			return -1;
		}

		for(int i = 0; i < n; ++i) {
			for(int k = 0; k < d; ++k) {
				X[i*d+k] = -5.0 + (double)rand() / RAND_MAX * 10.0;
			}
		}

		// printf("Matrix created\n");
		// printf("X=\n");
		// printMatrixDouble(X, n, d);

	} else if(rank != 0)
	{

		// printf("I am daughter no %d. Just waiting to calculate\n", rank);

	}

	// Find and print knnresult
	struct timespec start_time;
	if(rank == 0) {
		clock_gettime(CLOCK_MONOTONIC, &start_time);
	}

	knnresult kNNresult = distrAllkNN(X, n, d, k);

	if(rank == 0) {
		printf("\nV1 run time: %f\n", calculateExecutionTime(start_time));
	}
	// printf("%d nearest neighbors distances:\n", k);
	// printMatrixDouble(kNNresult.ndist, n, k);
	// printf("\n%d nearest neighbors indexes:\n", k);
	// printMatrixInt(kNNresult.nidx, n, k);

	if(rank == 0) {

		// printf("Done. Printing final kknresult\n");
		// printMatrixInt(kNNresult.nidx, n, k);
		// printMatrixDouble(kNNresult.ndist, n, k);

		// printf("Printing in Matlab format\n");
		// printMatrixMatlabFormatDouble(X, n, d, "X");
		// printMatrixMatlabFormatInt(kNNresult.nidx, n, k, "knn");
		// printf("k = %d;\n", k);
		// printf("n = %d;\n", n);
		// printf("WTF\n");
		// printf("knnCheck\n");

		free(kNNresult.nidx);
		free(kNNresult.ndist);
		free(X);

	} else {


	}
	// free(kNNresult.nidx);
	// free(kNNresult.ndist);

	MPI_Finalize();

	return 0;
}
