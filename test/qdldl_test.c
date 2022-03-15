
#include "csc.h"
#include "qdldl_solver.h"
#include "test_utils.h"
#include "simpletest/simpletest.h"

void CreateWorkspace() {
  KKTSystem kkt = GetTestSystem();
  QDLDLWorkspace ws = solvers_InitializeQDLDLWorkspace(&kkt);
  solvers_FreeQDLDLWorkspace(&ws);
  kkt_FreeKKTSystem(&kkt);
}

int main() {
  CreateWorkspace();
  PrintTestResult();
  return TestResult(); 
}
