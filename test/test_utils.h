#pragma once

#include "kkt.h"

double SumOfSquaredError(const double* x, const double* y, int len);

KKTSystem GetTestSystem();
