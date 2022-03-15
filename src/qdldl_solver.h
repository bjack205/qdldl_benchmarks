#pragma once

#include "qdldl.h"
#include "csc.h"

typedef struct {
    QDLDL_int n;
    QDLDL_int* Ap;       // n+1
    QDLDL_int* Ai;       // nnzA
    QDLDL_float* Ax;     // nnzA
    QDLDL_int* work;     // n
    QDLDL_int* Lnz;      // n
    QDLDL_int* etree;    // n
    QDLDL_int* Lp;       // n+1
    QDLDL_int* Li;       // nnzL
    QDLDL_float* Lx;     // nnzL
    QDLDL_float* D;      // n
    QDLDL_float* Dinv;   // n
    QDLDL_bool* bwork;   // n
    QDLDL_int* iwork;    // 3n
    QDLDL_float* fwork;  // n
} QDLDLWorkspace;

QDLDLWorkspace solvers_InitializeQDLDLWorkspace(const KKTSystem* kkt);

void solvers_FreeQDLDLWorkspace(QDLDLWorkspace* ws);
