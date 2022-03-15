#include <stdio.h>

#include "test_utils.h"
#include "csc.h"
#include "simpletest/simpletest.h"

void CheckTestSystem() { 
  KKTSystem kkt = GetTestSystem();
  const int nnz = 659;
  const int n = 159;
  TEST(kkt.n == n);
  TEST(kkt.nnz == nnz);
  TEST(csc_Nonzeros(&kkt.A) == nnz);
}

int main() {
  CheckTestSystem();
  PrintTestResult();
  return TestResult();
}