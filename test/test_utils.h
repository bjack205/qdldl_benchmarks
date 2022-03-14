#pragma once

#include "csc.h"

#ifdef JSONDATA
const char* kTestFile = JSONDATA;
#endif

double SumOfSquaredError(const double* x, const double* y, int len);

