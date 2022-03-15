#include <stdlib.h>
#include <stdio.h>

#include "pardiso_solver.h"
#include "utils.h"

int main(int argc, char* argv[]) {
  // Parse filename argument
  const char* filename;
  if (argc != 2) {
    printf("Usage Error: Expecting a single argument holding the json file name.\n");
  } else {
    filename = argv[1];
  }

  KKTSystem kkt = kkt_ReadFromFile(filename);
  double err = solvers_SolvePardiso(&kkt);
  printf("Pardiso solve error = %g\n", err);
  kkt_FreeKKTSystem(&kkt);

  return EXIT_SUCCESS;
}
