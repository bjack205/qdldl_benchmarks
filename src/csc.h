#pragma once

typedef struct {
    int n;
    int* colptr;
    int* rowval;
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
    int n;
    int nnz;
    SparseMatrixCSC A;
    double* b;
    double* x;
} KKTSystem;

KKTSystem kkt_ReadFromFile(const char* filename);

void kkt_FreeKKTSystem(KKTSystem* kkt);