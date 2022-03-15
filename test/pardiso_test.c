#include <stdio.h>

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
  solvers_FreePardisoWorkspace(&ws);

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

  int n = kkt.A.n;
  pardiso_chkmatrix(&mtype, &n, ws.a, ws.ia, ws.ja, &error);

  kkt_FreeKKTSystem(&kkt);
}

int main() {
  PardisoInit();
  PrintTestResult();
  return TestResult();
}
