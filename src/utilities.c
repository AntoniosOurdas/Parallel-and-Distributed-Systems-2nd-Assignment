#include "utilities.h"

#define BILLION 1000000000L;

// struct timespec start_time;
struct timespec stop_time;

double calculateExecutionTime(struct timespec start_time)
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
  for(int j = 0; j < m*2+3; ++j) {
    printf("-");
  }
  printf("\n");
  for(int i = 0; i < n; ++i) {
    printf("| ");
    for(int j = 0; j < m; ++j) {
			printf("%d ", A[i*m+j]);
		}
		printf("|\n");
	}
  for(int j = 0; j < m*2+3; ++j) {
    printf("-");
  }
  printf("\n");
}

void printMatrixDouble(double* A, int n, int m) {
  for(int j = 0; j < m*11; ++j) {
    printf("-");
  }
  printf("\n");
  for(int i = 0; i < n; ++i) {
    printf("| ");
    for(int j = 0; j < m; ++j) {
      if(A[i*m+j] < 0.0)
        printf("%lf ", A[i*m+j]);
      else
        printf(" %lf ", A[i*m+j]);
		}
		printf(" |\n");
	}
  for(int j = 0; j < m*11; ++j) {
    printf("-");
  }
  printf("\n");
}

void printMatrixMatlabFormatInt(int* A, int n, int m, char* name) {
  printf("%s = [", name);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      printf("%d ", A[i*m+j]);
    }
    if(i < n -1)
      printf(";");
  }
  printf("];\n");
}

void printMatrixMatlabFormatDouble(double* A, int n, int m, char* name) {
  printf("%s = [", name);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      printf("%lf ", A[i*m+j]);
    }
    if(i < n -1)
      printf(";");
  }
  printf("];\n");
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

int min(int a, int b) {
  return (a > b) ? b : a;
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
void findDistanceMatrix(double* X, double* Y, int n, int m, int d, double* D) {

  int blocksNo = min(m, 20);
  int block_size = m / blocksNo;
  double* X2 = (double*)malloc(n*d*sizeof(double));
	double* Y2 = (double*)malloc(block_size*d*sizeof(double));

  double* e1 = (double*)malloc(d*n*sizeof(double)); // [d-by-m] matrix with 1s (e*e^T)
  double* e2 = (double*)malloc(block_size*d*sizeof(double));	// [m-by-d] matrix with 1s (e*e^T)

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < d; ++j) {
      X2[i*d+j] = X[i*d+j]*X[i*d+j];
    }
  }

  for(int b = 0; b < blocksNo; ++b) {

  	for(int i = 0; i < block_size; ++i) {
  		for(int j = 0; j < d; ++j) {
  			Y2[i*d+j] = Y[b*d*block_size+i*d+j]*Y[b*d*block_size+i*d+j];
  		}
  	}

  	for(int i = 0; i < d; ++i) {
  		for(int j = 0; j < n; ++j) {
  			e1[i*n+j] = 1.0;
  		}
  	}

  	for(int i = 0; i < block_size; ++i) {
  		for(int j = 0; j < d; ++j) {
  			e2[i*d+j] = 1.0;
  		}
  	}
    // printf("b = %d\n", b);
  	// printMatrixDouble(X,n,d);
  	// printMatrixDouble(Y+b*d*block_size,block_size,d);

  	// printMatrixDouble(X2,n,d);
  	// printMatrixDouble(Y2,block_size,d);
  	// printMatrixDouble(e1,d,n);
  	// printMatrixDouble(e2,block_size,d);

    struct timespec start_time;

  	clock_gettime(CLOCK_MONOTONIC, &start_time);

  	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, block_size, n, d, 1, Y2, d, e1, n, 0.0, D+b*block_size*n, n);
  	// printMatrix(D, m, n);

  	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, block_size, n, d, -2.0, Y+b*block_size*d, d, X, d, 1.0, D+b*block_size*n, n);
  	// printMatrix(D, m, n);

  	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, block_size, n, d, 1.0, e2, d, X2, d, 1.0, D+b*block_size*n, n);
  	// printMatrix(D, m, n);

    // printf("D for b = %d\n", b);
    // printMatrixDouble(D+b*n*block_size, block_size, n);
  	// printf("BLAS run time: %f\n", calculateExecutionTime(start_time));

  }

  free(Y2);
  free(e1);
  free(e2);
  free(X2);

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
