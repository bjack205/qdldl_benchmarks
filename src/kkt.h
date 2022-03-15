#pragma once

#include "csc.h"

typedef struct {
    csc_int nprimals;
    csc_int nduals;
    SparseMatrixCSC A;
    double* b;
    double* x;
} KKTSystem;

KKTSystem kkt_ReadFromFile(const char* filename);

void kkt_FreeKKTSystem(KKTSystem* kkt);
