#include <cblas.h>
#include <stdio.h>

void printMatrix(double* A, int m, int n) {
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			printf("%lf ", A[i+j*m]);
		}
		printf("\n");
	}
	printf("\n");
}

int main()
{
  int i=0;
  int m = 3;
  int n = 3;
  int k = 2;
  double A[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};         
  double B[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};  
  double C[9] = {.5,.5,.5,.5,.5,.5,.5,.5,.5}; 
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,m,n,k,1,A, m, B, m,2,C,m);

  printMatrix(A, m, k);
  printMatrix(B, m, k);
  printMatrix(C, m, n);

  return 0;
}