#include <stdio.h>
#include <stdlib.h>

#include "csc.h"
#include "kkt.h"
#include "qdldl.h"
#include "qdldl_solver.h"
#include "qdldl_types.h"
#include "simpletest/simpletest.h"
#include "test_utils.h"

void CreateWorkspace() {
  KKTSystem kkt = GetTestSystem();
  QDLDLWorkspace ws = solvers_InitializeQDLDLWorkspace(&kkt);
  solvers_FreeQDLDLWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

void ComputeEtree() {
  KKTSystem kkt = GetTestSystem();
  QDLDLWorkspace ws = solvers_InitializeQDLDLWorkspace(&kkt);

  int nnzL = QDLDL_etree(ws.n, ws.Ap, ws.Ai, ws.work, ws.Lnz, ws.etree);
  TEST(nnzL > 0);
  int sum = 0;
  for (int i = 0; i < kkt.A.n; ++i) {
    sum += ws.Lnz[i];
  }
  TEST(sum == nnzL);

  solvers_FreeQDLDLWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

void Factorize() {
  KKTSystem kkt = GetTestSystem();
  QDLDLWorkspace ws = solvers_InitializeQDLDLWorkspace(&kkt);

  // Compute the elimination tree
  int nnzL = QDLDL_etree(ws.n, ws.Ap, ws.Ai, ws.work, ws.Lnz, ws.etree);
  TEST(nnzL > 0);

  // Allocate the rest of the data now that we know nnzL
  solvers_AppendQDLDLWorkspace(&ws);

  // Compute the factorization
  QDLDL_int res =
      QDLDL_factor(ws.n, ws.Ap, ws.Ai, ws.Ax, ws.Lp, ws.Li, ws.Lx, ws.D,
                   ws.Dinv, ws.Lnz, ws.etree, ws.bwork, ws.iwork, ws.fwork);
  TEST(res > 0);
  TEST(res == kkt.nprimals);

  solvers_FreeQDLDLWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

void Solve() {
  KKTSystem kkt = GetTestSystem();
  QDLDLWorkspace ws = solvers_InitializeQDLDLWorkspace(&kkt);

  // Compute the elimination tree
  int nnzL = QDLDL_etree(ws.n, ws.Ap, ws.Ai, ws.work, ws.Lnz, ws.etree);
  TEST(nnzL > 0);

  // Allocate the rest of the data now that we know nnzL
  solvers_AppendQDLDLWorkspace(&ws);

  // Compute the factorization
  QDLDL_int res =
      QDLDL_factor(ws.n, ws.Ap, ws.Ai, ws.Ax, ws.Lp, ws.Li, ws.Lx, ws.D,
                   ws.Dinv, ws.Lnz, ws.etree, ws.bwork, ws.iwork, ws.fwork);
  TEST(res > 0);
  TEST(res == kkt.nprimals);

  // Solve Ax = b
  const int n = kkt.A.n;
  QDLDL_float* x = (QDLDL_float*) malloc(n * sizeof(QDLDL_float));
  for (int i = 0; i < n; ++i) {
    x[i] = kkt.b[i];
  }
  double err = SumOfSquaredError(x, kkt.x, n);
  TEST(err > 0);

  QDLDL_solve(n, ws.Lp, ws.Li, ws.Lx, ws.Dinv, x);
  err = SumOfSquaredError(x, kkt.x, n);
  TEST(err < 1e-8);

  free(x);
  solvers_FreeQDLDLWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

int main() {
  CreateWorkspace();
  ComputeEtree();
  Factorize();
  Solve();
  PrintTestResult();
  return TestResult();
}
