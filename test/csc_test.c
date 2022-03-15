#include <stdio.h>

#include "test_utils.h"
#include "csc.h"
#include "kkt.h"
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

void CSRConversion() { 
  const int n = 5;
  const int nnz = 11;
  csc_int colptr[6] = {1,2,3,6,8,12};
  csc_int rowval[11] = {1,2,1,2,3,1,4,1,2,3,5};
  for (int i = 0; i < n + 1; ++i) {
    colptr[i] -= 1;
  }
  for (int i = 0; i < nnz; ++i) {
    rowval[i] -= 1;
  }
  double nzval[11] = {2.6708, 1.5248, -0.4627, 1.2921, 1.0, -0.6546, -1.4721, -0.4212, -0.423, -1.2201, -0.0934};
  SparseMatrixCSC A = {n, colptr, rowval, nzval};
  int ia[6];
  int ja[11];
  double a[11];
  csc_ConvertToCSR(&A, ia, ja, a);
  double nzval_ans[11] = {2.6708, -0.4627, -0.6546, -0.4212, 1.5248, 1.2921, -0.423, 1.0, -1.2201, -1.4721, -0.0934};
  TEST(SumOfSquaredError(a, nzval_ans, nnz) < 1e-6);

}

int main() {
  // CheckTestSystem();
  CSRConversion();
  PrintTestResult();
  return TestResult();
}
