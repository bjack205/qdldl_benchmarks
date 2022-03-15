#include "csc.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int csc_Nonzeros(const SparseMatrixCSC *A) { return A->colptr[A->n]; }

void csc_FreeSparseMatrixCSC(SparseMatrixCSC *A) {
  if (A->colptr) {
    free(A->colptr);
  }
  if (A->rowval) {
    free(A->rowval);
  }
  if (A->nzval) {
    free(A->nzval);
  }
}

void csc_ConvertToCSR(const SparseMatrixCSC* A, int* ia, int* ja, double* a) {
  const int n = A->n;
  ia[0] = 0;
  int cnt = 0;
  for (int r = 0; r < n; ++r) {
    // Search columns for an entry with the current row
    for (int c = 0; c < n; ++c) { 
      int istart = A->colptr[c];
      int istop = A->colptr[c + 1];
      for (int i = istart; i < istop; ++i) {
        if (A->rowval[i] == r) {
          a[cnt] = A->nzval[i];
          ja[cnt] = c;
          cnt += 1;
        }
      }
    }
    ia[r + 1] = cnt;
  }
}
