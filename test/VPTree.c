#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include "../src/utilities.h"
#define COUNT 10

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

typedef struct VPT {
  double* VP;       // Vantage point
  double median;    // Median distance
  struct VPT* in;   // Inner tree
  struct VPT* out;  // Outer tree
} VPT;

void addChild(VPT* node, int child, double* VP, double median) {

  // Insert point to left tree
  if(child == 1) {

    node->in = (VPT*)malloc(1*sizeof(VPT));
    node->in->in = NULL;
    node->in->out = NULL;

    node->in->VP = VP;
    // node->in->VP = (double*)malloc(d*sizeof(double));
    // for(int i = 0; i < d; ++i)
      // node->in->VP[i] = VP[i];

    node->in->median = median;

  } else if(child == 2) {

      node->out = (VPT*)malloc(1*sizeof(VPT));
      node->out->in = NULL;
      node->out->out = NULL;

      node->out->VP = VP;
      // node->out->VP = (double*)malloc(d*sizeof(double));
      // for(int i = 0; i < d; ++i)
      // node->out->VP[i] = VP[i];

      node->out->median = median;
  }
  return;
}

void printTree(VPT* T, int d) {
  if(T == NULL) {
    return;
  } else {
    printf("VP: ");
    for(int i = 0; i < d; ++i) {
      printf("%lf ", T->VP[i]);
    }
    printf("\n");
    printf("median: %lf\n", T->median);
  }
}

void print2DUtil(VPT* root, int space)
{
    // Base case
    if (root == NULL)
        return;

    // Increase distance between levels
    space += COUNT;

    // Process right child first
    print2DUtil(root->out, space);

    // Print current node after space
    // count
    printf("\n");
    for (int i = COUNT; i < space; i++)
        printf(" ");
    printf("%lf\n", *root->VP);

    // Process left child
    print2DUtil(root->in, space);
}

void print2D(VPT* root)
{
   // Pass initial space count as 0
   print2DUtil(root, 0);
}

void delete(VPT* T) {
  if(T == NULL) {
    return;
  } else {
    delete(T->in);
    delete(T->out);
    free(T);
    return;
  }
}

int partitionD(double* arr, int l, int r)
{
    double x = arr[r];
    int i = l;
    for (int j = l; j <= r - 1; j++) {
        if (arr[j] <= x) {
            swapDouble(&arr[i], &arr[j]);
            i++;
        }
    }
    swapDouble(&arr[i], &arr[r]);
    return i;
}

// This function returns k'th smallest
// element in arr[l..r] using QuickSort
// based method.  ASSUMPTION: ALL ELEMENTS
// IN ARR[] ARE DISTINCT
double kthSmallest(double* arr, int l, int r, int k)
{
  // If k is smaller than number of
  // elements in array
  if (k > 0 && k <= r - l + 1) {

      // Partition the array around last
      // element and get position of pivot
      // element in sorted array
      int index = partitionD(arr, l, r);

      // If position is same as k
      if (index - l == k - 1)
          return arr[index];

      // If position is more, recur
      // for left subarray
      if (index - l > k - 1)
        return kthSmallest(arr, l, index - 1, k);

      // Else recur for right subarray
      return kthSmallest(arr, index + 1, r, k - index + l - 1);
  }

  // If k is more than number of
  // elements in array
  printf("You are asking for more than available\n");
  exit(1);
}

void makeVPT(VPT* T, double* X, int n, int d) {
    T->VP = X;
    if(n == 1) {
      return;
    } else {

      // Find all distances from VP
      // except from VP of course
      double* dVP = (double*)malloc((n-1)*sizeof(double));
      for(int i = 1; i < n; ++i) {
        dVP[i] = 0;
        for(int j = 0; j < d; ++j) {
          dVP[i] += (X[i*d+j] - T->VP[j])*(X[i*d+j] - T->VP[j]);
        }
        dVP[i] = sqrt(dVP[i]);
      }

      printMatrixDouble(dVP, n, 1);

      // Find median and create subtrees
      T->median = kthSmallest(X, 0, n-1, (n-1)/ 2);
      printf("median = %lf\n", T->median);

      double* X1 = NULL;
      double* X2 = NULL;

      int sizeX1 = 0;
      int sizeX2 = 0;

      for(int i = 1; i < n; ++i) {
        if(dVP[i] < T->median) {

          ++sizeX1;
          X1 = (double*)realloc(X1, sizeX1*d*sizeof(double));
          for(int j = 0; j < d; ++j) {
            X1[(sizeX1-1)*d+j] = X[i*d+j];
          }

        } else {

          ++sizeX2;
          X2 = (double*)realloc(X2, sizeX2*d*sizeof(double));
          for(int j = 0; j < d; ++j) {
            X2[(sizeX2-1)*d+j] = X[i*d+j];
          }

        }
      }

      // Create left subtree for nearest neighbors
      T->in = (VPT*)malloc(1*sizeof(VPT));
      // makeVPT(T->in, X1, sizeX1, d);

      // Create right subtree for furthest neighbors
      T->out = (VPT*)malloc(1*sizeof(VPT));
      // makeVPT(T->out, X2, sizeX2, d);

      free(X1);
      free(X2);
      return;
  }
}

int main(int argc, char* argv[]) {

  int n = atoi(argv[1]);
  int d = atoi(argv[2]);

  VPT* T = (VPT*)malloc(1*sizeof(VPT));

  double* X = (double*)malloc(n*d*sizeof(double));
  int* indices = (int*)malloc(n*sizeof(int));
  for(int i = 0; i < n; ++i)
    indices[i] = i;

  for(int i = 0; i < n; ++i)
    for(int j = 0; j < d; ++j)
      X[i*d+j] = (double)rand()/RAND_MAX;

  T->VP = X;
  T->in = NULL;
  T->out = NULL;
  // print2D(T);
  makeVPT(T, X, n, d);
  
  delete(T);
  free(X);

  // printf("%d %d\n", (n-1)/2, (n-1)/2+(n-1)%2);
  return 0;

}
