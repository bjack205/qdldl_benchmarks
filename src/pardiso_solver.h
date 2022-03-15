#include "csc.h"
#include "kkt.h"

// Forward-declare PARDISO functions
void pardisoinit(void* pt, int* mtype, int* solver, int* iparm, double* dparm, int* error);
void pardiso     (void* pt, int* maxfct, int* mnum, int* mtype,  int* phase, int* n, 
                  double* a, int* ia, int* ja, int* perm, int* nrhs, int* iparm,
                  int* msglvl, double* b, double* x, int* error, double* dparm);
void pardiso_chkmatrix  (int* mtype, int* n, double* a, int* ia, int* ja, int* error);
void pardiso_chkvec     (int* n, int* nrhs, double* b, int* error);
void pardiso_printstats (int* mtype, int* n, double* a, int* ia, int* ja, int* nrhs,
                           double* b, int* error);


enum MTYPE {
   PARDISO_RSYM = 1,     // real symmetric
   PARDISO_RPSD = 2,     // real positive definite
   PARDISO_RINDEF = -2,  // real symmetric indefinite
   PARDISO_CSYM = 3,     // complex structurally symmetric
   PARDISO_CPSD = 4,     // complex Hermitian positive definite
   PARDISO_CINDEF = -4,  // complex Hermitian indefinite
   PARDISO_CSYM2 = 6,    // complex symmetric
   PARDISO_RGEN = 11,    // real nonsymmetric
   PARDISO_CGEN = 13,    // complex nonsymmetric
};

enum PARDISO_SOLVER {
   PARDISO_SPARSE_DIRECT_SOLVER = 0,
   PARDISO_RECURSIVE_ITERATIVE_SOLVER = 1,
};

enum PARDISO_INIT_ERROR {
   PARDISO_NO_ERROR = 0,
   PARDISO_NO_LICENSE = -10,
   PARDISO_LICENSE_EXPIRED = -11,
   PARDISO_WRONG_USERNAME = -12, 
};

enum PARDISO_PHASE {
   PARDISO_ANALYSIS = 11,
   PARDISO_ANALYSIS_FACTOR = 12,
   PARDISO_ANALYSIS_FACTOR_SOLVE_REFINE = 13,
   PARDISO_FACTOR = 22,
   PARDISO_INVERSION = -22,
   PARDISO_FACTOR_SOLVE_REFINE = 23,
   PARDISO_SOLVE_REFINE = 33,
   PARDISO_RELEASE_LU = 0,
   PARDISO_RELEASE_ALL = -1,
};

typedef struct {
   void* pt[64];
   int iparm[64];
   double dparm[64];
   double* a;          // (nnz,)
   int* ia;            // (n+1,)
   int* ja;            // (nnz,)
   int* perm;          // (n,)
   double* b;          // (n, nrhs)
   double* x;          // (n, nrhs)
} PardisoWorkspace;

PardisoWorkspace solvers_InitializePardisoWorkspace(const KKTSystem* kkt);

void solvers_FreePardisoWorkspace(PardisoWorkspace* ws);

/**
 * @brief Get the value of OMP_NUM_THREADS
 * 
 * @return The value of OMP_NUM_THREADS, or 1 if not defined
 */
int solvers_GetOmpThreads();

double solvers_SolvePardiso(const KKTSystem* kkt);
