#include "qdldl_solver.h"

#include <stdlib.h>

#include "csc.h"
#include "qdldl_types.h"

QDLDLWorkspace solvers_InitializeQDLDLWorkspace(const KKTSystem* kkt) {
  const SparseMatrixCSC* A = &kkt->A;
  int n = kkt->n;
  int nnzA = csc_Nonzeros(A);

  // Allocate integer data
  int idata_size = (n + 1) + nnzA + 3 * n;
  QDLDL_int* idata = (QDLDL_int*) malloc(idata_size * sizeof(QDLDL_int));

  // Partition integer data 
  QDLDL_int* Ap = idata;
  QDLDL_int* Ai = idata + (n + 1);
  QDLDL_int* work = Ai + nnzA; 
  QDLDL_int* Lnz = work + n; 
  QDLDL_int* etree = Lnz + n; 

  // Allocate Floating data
  int fdata_size = nnzA;
  QDLDL_float* fdata = (QDLDL_float*) malloc(fdata_size * sizeof(QDLDL_float));

  // Partition float data
  QDLDL_float* Ax = fdata;

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
  free(ws->Ap);
  free(ws->Ax);
}
