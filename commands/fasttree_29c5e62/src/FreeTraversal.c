#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

traversal_t FreeTraversal(traversal_t traversal, NJ_t *NJ) {
  myfree(traversal, sizeof(bool)*NJ->maxnodes);
  return(NULL);
}