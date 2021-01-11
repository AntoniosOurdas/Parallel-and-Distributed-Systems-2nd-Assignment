#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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


void printMatrix(double* A, int n, int m) {
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

void printMatrixMatlab(double* A, int m, int n) {
	printf("[");
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			printf("%lf ", A[i*n+j]);
		}
		printf(";\n");
	}
	printf("]\n");
}

void multiply(double* A, double* B, double* C, int m, int n, int k) {
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			// C[i*n+j] = 0.0;
			for(int l = 0; l < k; ++l) {
				C[i*n+j] += A[i*k+l]*B[l*n+j];
			}
		}
	}
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
  int m = atoi(argv[1]);
	int k = atoi(argv[2]);
  int n = atoi(argv[3]);
	//
	// double* A = (double*)malloc(m*k*sizeof(double));
	// double* B = (double*)malloc(k*n*sizeof(double));
	// double* C = (double*)malloc(m*n*sizeof(double));
	//
	// for(int i = 0; i < m; ++i) {
	// 	for(int j = 0; j < k; ++j) {
	// 		// A[i*k+j] = (double)(i+1);
	// 		A[i*k+j] = -(double)rand()/RAND_MAX*10.0+5.0;
	// 	}
	// }
	//
	// for(int i = 0; i < k; ++i) {
	// 	for(int j = 0; j < n; ++j) {
	// 		// B[i*n+j] = (double)(i+1);
	// 		B[i*n+j] = -(double)rand()/RAND_MAX*10.0+5.0;;
	// 	}
	// }
	//
	// for(int i = 0; i < m; ++i) {
	// 	for(int j = 0; j < n; ++j) {
	// 		C[i*n+j] = -(double)rand()/RAND_MAX*10.0+5.0;
	// 	}
	// }
	//
	// // printMatrixMatlab(A, m, k);
	// // printMatrixMatlab(B, k, n);
	// // printMatrixMatlab(C, m, n);
	// struct timespec start_time;
	// clock_gettime(CLOCK_MONOTONIC, &start_time);
	// // multiply(A, B, C, m, n, k);
	// printf("Simple run time: %f\n", calculateExecutionTime(start_time));
	//
	// // struct timespec start_time;
	// clock_gettime(CLOCK_MONOTONIC, &start_time);
	// cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 1, C, n);
	// printf("BLAS run time: %f\n", calculateExecutionTime(start_time));

	int d = 2;
	double* X = (double*)malloc(n*d*sizeof(double));
	double* Y = (double*)malloc(m*d*sizeof(double));

	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < d; ++j) {
			Y[i*d+j] = (double)rand()/RAND_MAX*10.0-5.0;
		}
	}

	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < d; ++j) {
			X[i*d+j] = (double)rand()/RAND_MAX*10.0-5.0;
		}
	}

	double* X2 = (double*)malloc(n*d*sizeof(double));
	double* Y2 = (double*)malloc(m*d*sizeof(double));

	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < d; ++j) {
			X2[i*d+j] = X[i*d+j]*X[i*d+j];
		}
	}

	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < d; ++j) {
			Y2[i*d+j] = Y[i*d+j]*Y[i*d+j];
		}
	}

	double* e1 = (double*)malloc(d*n*sizeof(double)); // [d-by-m] matrix with 1s (e*e^T)
	double* e2 = (double*)malloc(m*d*sizeof(double));	// [m-by-d] matrix with 1s (e*e^T)

	for(int i = 0; i < d; ++i) {
		for(int j = 0; j < n; ++j) {
			e1[i*n+j] = 1.0;
		}
	}

	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < d; ++j) {
			e2[i*d+j] = 1.0;
		}
	}

	double* D = (double*)malloc(n*m*sizeof(double));

	// printMatrixMatlab(X	,n,d);
	// printMatrixMatlab(Y,n,d);

	// printMatrix(X2,n,d);
	// printMatrix(Y2,m,d);
	// printMatrix(e1,d,n);
	// printMatrix(e2,m,d);

	struct timespec start_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, d, 1, Y2, d, e1, n, 0.0, D, n);
	// printMatrix(D, m, n);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, d, -2.0, Y, d, X, d, 1.0, D, n);
	// printMatrix(D, m, n);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, d, 1.0, e2, d, X2, d, 1.0, D, n);
	// printMatrix(D, m, n);

	printf("BLAS run time: %f\n", calculateExecutionTime(start_time));
	free(X);
	free(Y);
	free(X2);
	free(Y2);
	free(e1);
	free(e2);
	free(D);

  return 0;
}
