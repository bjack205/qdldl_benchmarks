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

/**
 * @brief Convert a CSC matrix to CSR
 * 
 * @param A Original matrix
 * @param ia (n+1,) vector of row pointers
 * @param ja (nnz,) vector of column values
 * @param a  (nnz,) vector of nonzero values
 */
void csc_ConvertToCSR(const SparseMatrixCSC* A, int* ia, int* ja, double* a);
