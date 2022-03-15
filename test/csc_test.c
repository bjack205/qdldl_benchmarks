#include <stdio.h>

#include "test_utils.h"
#include "csc.h"
#include "simpletest/simpletest.h"

void CheckTestSystem() { 
  KKTSystem kkt = GetTestSystem();
  const int nnz = 409;
  const int n = 159;
  const int nprimals = 95;
  TEST(kkt.A.n == n);
  TEST(csc_Nonzeros(&kkt.A) == nnz);
  TEST(kkt.nprimals == nprimals);
  TEST(kkt.nduals == n - nprimals);
  kkt_FreeKKTSystem(&kkt);
}

int main() {
  CheckTestSystem();
  PrintTestResult();
  return TestResult();
}
