#pragma once

#include "qdldl.h"
#include "csc.h"
#include "kkt.h"

typedef struct {
    QDLDL_int n;
    const QDLDL_int* const Ap;       // n+1
    const QDLDL_int* const Ai;       // nnzA
    const QDLDL_float* const Ax;     // nnzA
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

/**
 * @brief Initializes the rest of the workspace needed by QDLDL after the 
 *        number of nonzeros in L is known.
 * 
 * Must be called after the elimination tree is computed.
 * 
 * @param ws 
 */
void solvers_AppendQDLDLWorkspace(QDLDLWorkspace* ws);

double solvers_SolveQDLDL(const KKTSystem* kkt);
