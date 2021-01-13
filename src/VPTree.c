#include "VPTree.h"

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

int partitionQS(double* arr, int l, int r)
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
      int index = partitionQS(arr, l, r);

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
        dVP[i-1] = 0.0;
        for(int j = 0; j < d; ++j) {
          dVP[i-1] += (X[i*d+j] - T->VP[j])*(X[i*d+j] - T->VP[j]);
        }
        dVP[i-1] = sqrt(dVP[i-1]);
      }

      double* temp = (double*)malloc((n-1)*sizeof(double));
      for(int i = 0; i < n-1; ++i) {
        temp[i] = dVP[i];
      }

      // Find median and create subtrees
      T->median = kthSmallest(temp, 0, n-2, (n-1)/2);

      free(temp);

      double* X1 = NULL;
      double* X2 = NULL;

      int sizeX1 = 0;
      int sizeX2 = 0;
      double* tempX = NULL;

      for(int i = 1; i < n; ++i) {
        if(dVP[i-1] < T->median) {
          ++sizeX1;
          if(sizeX1 == 1) {
            X1 = (double*)malloc(sizeX1*d*sizeof(double));
          } else {
            tempX = X1;
            X1 = (double*)malloc(sizeX1*d*sizeof(double));
            for(int l = 0; l < (sizeX1-1); ++l) {
              for(int j = 0; j < d; ++j) {
                X1[l*d+j] = tempX[l*d+j];
              }
            }
            free(tempX);
          }

          for(int j = 0; j < d; ++j) {
            X1[(sizeX1-1)*d+j] = X[i*d+j];
          }

        } else if(dVP[i-1] > T->median){
          ++sizeX2;
          if(sizeX2 == 1) {
            X2 = (double*)malloc(sizeX2*d*sizeof(double));
          } else {
            tempX = X2;
            X2 = (double*)malloc(sizeX2*d*sizeof(double));
            for(int l = 0; l < (sizeX2-1); ++l) {
              for(int j = 0; j < d; ++j) {
                X2[l*d+j] = tempX[l*d+j];
              }
            }
            free(tempX);
          }

          for(int j = 0; j < d; ++j) {
            X2[(sizeX2-1)*d+j] = X[i*d+j];
          }
        }
      }

      // Create left subtree for nearest neighbors
      T->in = (VPT*)malloc(1*sizeof(VPT));
      makeVPT(T->in, X1, sizeX1, d);

      // Create right subtree for furthest neighbors
      T->out = (VPT*)malloc(1*sizeof(VPT));
      makeVPT(T->out, X2, sizeX2, d);

      free(X1);
      free(X2);
      return;
  }
}

void searchVPT(VPT *T, double *X, int d, double* ndist, int* nidx) {
  
}
