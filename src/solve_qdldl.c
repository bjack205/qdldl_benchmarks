#include <stdio.h>
#include <stdlib.h>

#include "csc.h"
#include "kkt.h"
#include "utils.h"
#include "qdldl_solver.h"


int main(int argc, char* argv[]) {
  // Parse filename argument
  char* filename = "";
  if (argc != 2) {
    printf("Usage Error: Expecting a single argument holding the json file name.\n");
  } else {
    filename = argv[1];
  }

  KKTSystem kkt = kkt_ReadFromFile(filename);
  solvers_SolveQDLDL(&kkt);
  kkt_FreeKKTSystem(&kkt);

  return EXIT_SUCCESS;
}
