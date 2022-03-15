#include "pardiso_solver.h"

#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "csc.h"

PardisoWorkspace solvers_InitializePardisoWorkspace(const KKTSystem* kkt) {   
  (void) kkt;

  const SparseMatrixCSC* A = &kkt->A;
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
  ws.iparm[0] = 0; // Use default params
  ws.iparm[2] = solvers_GetOmpThreads();
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

double solvers_SolvePardiso(const KKTSystem* kkt) {
  const SparseMatrixCSC *A = &kkt->A;
  (void)A;

  PardisoWorkspace ws = solvers_InitializePardisoWorkspace(kkt);
  int mtype = PARDISO_RINDEF;
  int solver = PARDISO_SPARSE_DIRECT_SOLVER;
  int n = kkt->A.n;
  int nrhs = 1;
  int maxfct = 1; // Maximum number of numerical factorization
  int mnum = 1;   // which factorization to use
  int msglvl = 0; // Print statistical information
  int error = 0;

  // Initialization
  pardisoinit(ws.pt, &mtype, &solver, ws.iparm, ws.dparm, &error);

  // Parameters
  ws.iparm[1] = 0;  // Reordering. I found this makes significant difference (default = 2) 
  ws.iparm[32] = 1; // Compute the determinant
  ws.iparm[7] = 1;  // Max number of iterative refinement steps (default = 0)

  // Analysis (symbolic factorization)
  int phase = PARDISO_ANALYSIS;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);

  // Numerical factorization
  phase = PARDISO_FACTOR;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);

  // Back substitution and interative refinement
  phase = PARDISO_SOLVE_REFINE;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);

  double err = SumOfSquaredError(ws.x, kkt->x, n);

  // Free all memory (still has a memory leak)
  phase = PARDISO_RELEASE_ALL;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);

  solvers_FreePardisoWorkspace(&ws);

  return err;
}
