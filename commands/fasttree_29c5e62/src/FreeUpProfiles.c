#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

profile_t **FreeUpProfiles(profile_t **upProfiles, NJ_t *NJ) {
  int i;
  int nUsed = 0;
  for (i=0; i < NJ->maxnodes; i++) {
    if (upProfiles[i] != NULL)
      nUsed++;
    DeleteUpProfile(upProfiles, NJ, i);
  }
  myfree(upProfiles, sizeof(profile_t*)*NJ->maxnodes);
  if (verbose >= 3)
    fprintf(stderr,"FreeUpProfiles -- freed %d\n", nUsed);
  return(NULL);
}