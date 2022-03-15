#include "csc.h"

// Forward-declare PARDISO functions
void pardisoinit(void* pt, int* mtype, int* solver, int* iparm, double* dparm, int* error);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);


enum MTYPE {
   RSYM = 1,     // real symmetric
   RPSD = 2,     // real positive definite
   RINDEF = -2,  // real symmetric indefinite
   CSYM = 3,     // complex structurally symmetric
   CPSD = 4,     // complex Hermitian positive definite
   CINDEF = -4,  // complex Hermitian indefinite
   CSYM2 = 6,    // complex symmetric
   RGEN = 11,    // real nonsymmetric
   CGEN = 13,    // complex nonsymmetric
};

enum PARDISO_SOLVER {
   SPARSE_DIRECT_SOLVER = 0,
   RECURSIVE_ITERATIVE_SOLVER = 1,
};

enum PARDISO_ERROR {
   PARDISO_NO_ERROR = 0,
   PARDISO_NO_LICENSE = -10,
   PARDISO_LICENSE_EXPIRED = -11,
   PARDISO_WRONG_USERNAME = -12, 
};

typedef struct {
   void* pt;        // (64,)
   int* iparm;      // (64,)
   double* dparm;   // (64,)
} PardisoWorkspace;

PardisoWorkspace solvers_InitializePardisoWorkspace(KKTSystem* kkt);
