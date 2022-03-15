#pragma once

#include <inttypes.h>
#include <stdint.h>

typedef long long csc_int;

typedef struct {
    csc_int n;
    csc_int* colptr;
    csc_int* rowval;
    double* nzval;
} SparseMatrixCSC;

/**
 * @brief Get the number of nonzeros in the sparse matrix
 * 
 * @param A 
 * @return Number of nonzeros
 */
int csc_Nonzeros(const SparseMatrixCSC* A);

void csc_FreeSparseMatrixCSC(SparseMatrixCSC* A);

typedef struct {
    csc_int nprimals;
    csc_int nduals;
    SparseMatrixCSC A;
    double* b;
    double* x;
} KKTSystem;

KKTSystem kkt_ReadFromFile(const char* filename);

void kkt_FreeKKTSystem(KKTSystem* kkt);
