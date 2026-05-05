#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

profile_t *DeleteUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node) {
  assert(node>=0 && node < NJ->maxnodes);
  if (upProfiles[node] != NULL)
    upProfiles[node] = FreeProfile(upProfiles[node], NJ->nPos, NJ->nConstraints); /* returns NULL */
  return(NULL);
}