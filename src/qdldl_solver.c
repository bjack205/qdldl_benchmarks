#include "qdldl_solver.h"

#include <stdlib.h>

#include "csc.h"
#include "qdldl_types.h"

QDLDLWorkspace solvers_InitializeQDLDLWorkspace(const KKTSystem* kkt) {
  const SparseMatrixCSC* A = &kkt->A;
  int n = kkt->n;
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
}
