#include <stdio.h>

#include "csc.h"
#include "qdldl.h"
#include "qdldl_solver.h"
#include "test_utils.h"
#include "simpletest/simpletest.h"

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
  for (int i = 0; i < kkt.n; ++i) {
    sum += ws.Lnz[i];
  }
  TEST(sum == nnzL);

  solvers_FreeQDLDLWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

int main() {
  CreateWorkspace();
  ComputeEtree();
  PrintTestResult();
  return TestResult(); 
}
