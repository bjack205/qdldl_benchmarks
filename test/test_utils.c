#include "test_utils.h"
#include "csc.h"

#include <stdio.h>
#include <math.h>

#ifdef JSONDATA
const char* kTestFile = JSONDATA;
#endif

double SumOfSquaredError(const double* x, const double* y, int len) {
  double err = 0;
  for (int i = 0; i < len; ++i) {
    double diff = x[i] - y[i];
    err += diff * diff;
  }
  return sqrt(err);
}

KKTSystem GetTestSystem() {
  return kkt_ReadFromFile(kTestFile);
}
