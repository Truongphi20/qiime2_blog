#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

int ActiveAncestor(/*IN*/NJ_t *NJ, int iNode) {
  if (iNode < 0)
    return(iNode);
  while(NJ->parent[iNode] >= 0)
    iNode = NJ->parent[iNode];
  return(iNode);
}