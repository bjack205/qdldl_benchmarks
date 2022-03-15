#include "pardiso_solver.h"

#include <stdlib.h>
#include <stdio.h>

PardisoWorkspace solvers_InitializePardisoWorkspace(KKTSystem* kkt) {
  (void) kkt;

  void* pt[64];
  int iparm[64];
  double dparm[64];
  PardisoWorkspace ws = {
    .pt = pt,
    .iparm = iparm,
    .dparm = dparm,
  };
  return ws;
}

int solvers_GetOmpThreads() {
  int num_procs;
  char* var = getenv("OMP_NUM_THREADS");
  if (var != NULL) {
    sscanf(var, "%d", &num_procs);
  } else {
    num_procs = 1;
  }
  return num_procs;
}
