#include "csc.h"
#include "simpletest/simpletest.h"
#include "csc.h"
#include "pardiso_solver.h"
#include "test_utils.h"

void PardisoInit() {
  KKTSystem kkt = GetTestSystem();

  PardisoWorkspace ws = solvers_InitializePardisoWorkspace(&kkt);
  int mtype = RINDEF;
  int solver = SPARSE_DIRECT_SOLVER;
  int error;
  pardisoinit(ws.pt, &mtype, &solver, ws.iparm, ws.dparm, &error);
  TEST(error == PARDISO_NO_ERROR);

  kkt_FreeKKTSystem(&kkt);
}

void PardisoChecks() {
  KKTSystem kkt = GetTestSystem();

  PardisoWorkspace ws = solvers_InitializePardisoWorkspace(&kkt);
  int mtype = RINDEF;
  int solver = SPARSE_DIRECT_SOLVER;
  int error;
  pardisoinit(ws.pt, &mtype, &solver, ws.iparm, ws.dparm, &error);
  TEST(error == PARDISO_NO_ERROR);

  // SparseMatrixCSC* A = &kkt.A;
  // pardiso_chkmatrix(&mtype, &A->n, A->nzval, A->colptr, A->rowval, &error);

  kkt_FreeKKTSystem(&kkt);
}

int main() {
  PrintTestResult();
  return TestResult();
}
