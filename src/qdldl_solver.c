#include "qdldl_solver.h"

#include <stdlib.h>

#include "csc.h"
#include "qdldl_types.h"

QDLDLWorkspace solvers_InitializeQDLDLWorkspace(const KKTSystem* kkt) {
  const SparseMatrixCSC* A = &kkt->A;
  const int n = A->n; 
  // int nnzA = csc_Nonzeros(A);

  // Assign pointers to data in A 
  const QDLDL_int* const Ap = A->colptr; 
  const QDLDL_int* const Ai = A->rowval; 
  const QDLDL_float* const Ax = A->nzval; 

  // Allocate integer data
  int idata_size = 3 * n;
  QDLDL_int* idata = (QDLDL_int*) malloc(idata_size * sizeof(QDLDL_int));

  // Partition integer data 
  QDLDL_int* work = idata;
  QDLDL_int* Lnz = work + n; 
  QDLDL_int* etree = Lnz + n; 

  // Partially initialize the workspace
  QDLDLWorkspace ws = {
    .n = n,
    .Ap = Ap,
    .Ai = Ai,
    .Ax = Ax,
    .work = work,
    .Lnz = Lnz,
    .etree = etree,
    .Lp = NULL,
    .Li = NULL,
    .Lx = NULL,
    .D = NULL,
    .Dinv = NULL,
    .bwork = NULL,
    .iwork = NULL,
    .fwork = NULL,
  };
  return ws;
}

void solvers_FreeQDLDLWorkspace(QDLDLWorkspace* ws) {
  free(ws->work);
  if (ws->Lp) {
    free(ws->Lp);
    free(ws->Lx);
    free(ws->bwork);
  }
}

void solvers_AppendQDLDLWorkspace(QDLDLWorkspace* ws) { 
  const int n = ws->n;
  QDLDL_int nnzL = 0; 
  for (int i = 0; i < ws->n; ++i) {
    nnzL += ws->Lnz[i];
  }

  // Allocate integer data
  int idata_size = (n + 1) + nnzL + 3 * n;
  QDLDL_int* idata = (QDLDL_int*) malloc(idata_size * sizeof(QDLDL_int));

  // Parition integer data
  ws->Lp = idata;
  ws->Li = ws->Lp + (n + 1);
  ws->iwork = ws->Li  + nnzL;

  // Allocate float data
  int fdata_size = nnzL + 3 * n;
  QDLDL_float* fdata = (QDLDL_float*) malloc(fdata_size * sizeof(QDLDL_float));

  // Partition float data
  ws->Lx = fdata; 
  ws->D = ws->Lx + nnzL;
  ws->Dinv = ws->D + n;
  ws->fwork = ws->Dinv + n;

  // Create boolean data
  ws->bwork = (QDLDL_bool*) malloc(n * sizeof(QDLDL_bool));
}
