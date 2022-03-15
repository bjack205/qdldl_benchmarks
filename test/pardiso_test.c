#include <stdio.h>

#include "csc.h"
#include "pardiso_solver.h"
#include "simpletest/simpletest.h"
#include "test_utils.h"

void PardisoInit() {
  KKTSystem kkt = GetTestSystem();

  PardisoWorkspace ws = solvers_InitializePardisoWorkspace(&kkt);
  int mtype = PARDISO_RINDEF;
  int solver = PARDISO_SPARSE_DIRECT_SOLVER;
  int error;

  pardisoinit(ws.pt, &mtype, &solver, ws.iparm, ws.dparm, &error);
  TEST(error == PARDISO_NO_ERROR);
  solvers_FreePardisoWorkspace(&ws);

  kkt_FreeKKTSystem(&kkt);
}

void PardisoChecks() {
  KKTSystem kkt = GetTestSystem();
  SparseMatrixCSC *A = &kkt.A;
  (void)A;

  PardisoWorkspace ws = solvers_InitializePardisoWorkspace(&kkt);
  int mtype = PARDISO_RINDEF;
  int solver = PARDISO_SPARSE_DIRECT_SOLVER;
  int error;

  pardisoinit(ws.pt, &mtype, &solver, ws.iparm, ws.dparm, &error);
  TEST(error == PARDISO_NO_ERROR);

  int n = kkt.A.n;
  pardiso_chkmatrix(&mtype, &n, ws.a, ws.ia, ws.ja, &error);
  TEST(error == 0);

  int nrhs = 1;
  pardiso_chkvec(&n, &nrhs, ws.b, &error);
  TEST(error == 0);

  pardiso_printstats(&mtype, &n, ws.a, ws.ia, ws.ja, &nrhs, ws.b, &error);
  TEST(error == 0);

  solvers_FreePardisoWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

void PardisoSolve() {
  KKTSystem kkt = GetTestSystem();
  SparseMatrixCSC *A = &kkt.A;
  (void)A;

  PardisoWorkspace ws = solvers_InitializePardisoWorkspace(&kkt);
  int mtype = PARDISO_RINDEF;
  int solver = PARDISO_SPARSE_DIRECT_SOLVER;
  int n = kkt.A.n;
  int nrhs = 1;
  int maxfct = 1; // Maximum number of numerical factorization
  int mnum = 1;   // which factorization to use
  int msglvl = 1; // Print statistical information
  int error = 0;

  // Initialization
  pardisoinit(ws.pt, &mtype, &solver, ws.iparm, ws.dparm, &error);
  TEST(error == PARDISO_NO_ERROR);

  // Analysis (symbolic factorization)
  int phase = PARDISO_ANALYSIS;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);
  TEST(error == 0);

  // Numerical factorization
  phase = PARDISO_FACTOR;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);
  TEST(error == 0);

  // Back substitution and interative refinement
  phase = PARDISO_SOLVE_REFINE;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);
  TEST(error == 0);

  double err = SumOfSquaredError(ws.x, kkt.x, n);
  TEST(err < 1e-10);

  // Free all memory (still has a memory leak)
  phase = PARDISO_RELEASE_ALL;
  pardiso(ws.pt, &maxfct, &mnum, &mtype, &phase, &n, ws.a, ws.ia, ws.ja,
          ws.perm, &nrhs, ws.iparm, &msglvl, ws.b, ws.x, &error, ws.dparm);

  solvers_FreePardisoWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

int main() {
  PardisoInit();
  PardisoChecks();
  PardisoSolve();
  PrintTestResult();
  return TestResult();
}
