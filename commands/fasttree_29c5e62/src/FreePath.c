#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

int *FreePath(int *path, NJ_t *NJ) {
  myfree(path, sizeof(int)*NJ->maxnodes);
  return(NULL);
}