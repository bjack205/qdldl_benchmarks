#include "pardiso_solver.h"

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
