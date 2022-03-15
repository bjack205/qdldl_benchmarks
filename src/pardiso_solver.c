#include "pardiso_solver.h"
#include "csc.h"

#include <stdlib.h>
#include <stdio.h>

PardisoWorkspace solvers_InitializePardisoWorkspace(KKTSystem* kkt) {   
  (void) kkt;

  SparseMatrixCSC* A = &kkt->A;
  const int n = A->n;
  const int nnz = csc_Nonzeros(A);

  // Convert data in A to upper triangular CSR 
  double* a = (double*) malloc(nnz * sizeof(double));
  int* ia = (int*) malloc((n + 1) * sizeof(int));
  int* ja = (int*) malloc(nnz * sizeof(int));

  csc_ConvertToCSR(A, ia, ja, a);

  // Convert to 1-based indexing
  for (int i = 0; i < n + 1; ++i) {
    ia[i] += 1;
  }
  for (int i = 0; i < nnz; ++i) {
    ja[i] += 1;
  }

  // Create permuation vector
  int* perm = (int*) malloc(n * sizeof(int));
  for (int i = 0; i < n; ++i) {
    perm[i] = i + 1;
  }

  // Allocate x and b arrays
  const int nrhs = 1;
  double* b = (double*) malloc(nrhs * n * sizeof(double));
  double* x = (double*) malloc(nrhs * n * sizeof(double));
  for (int i = 0; i < n; ++i) {
    b[i] = kkt->b[i];
    x[i] = 0.0;
  }

  PardisoWorkspace ws = {
    .a = a,
    .ia = ia,
    .ja = ja,
    .perm = perm,
    .b = b,
    .x = x,
  };
  return ws;
}

void solvers_FreePardisoWorkspace(PardisoWorkspace* ws) {
  free(ws->a);
  free(ws->ia);
  free(ws->ja);
  free(ws->perm);
  free(ws->b);
  free(ws->x);
}

int solvers_GetOmpThreads() {
  int num_procs;
  char* var = getenv("OMP_NUM_THREADS");
  if (var != NULL) {
    sscanf(var, "%d", &num_procs);
  } else {
    num_procs = 1;
  }
  return num_procs;
}
