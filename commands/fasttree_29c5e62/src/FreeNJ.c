#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


NJ_t *FreeNJ(NJ_t *NJ) {
  if (NJ==NULL)
    return(NJ);

  int i;
  for (i=0; i < NJ->maxnode; i++)
    NJ->profiles[i] = FreeProfile(NJ->profiles[i], NJ->nPos, NJ->nConstraints);
  NJ->profiles = myfree(NJ->profiles, sizeof(profile_t*) * NJ->maxnodes);
  NJ->outprofile = FreeProfile(NJ->outprofile, NJ->nPos, NJ->nConstraints);
  NJ->diameter = myfree(NJ->diameter, sizeof(numeric_t)*NJ->maxnodes);
  NJ->varDiameter = myfree(NJ->varDiameter, sizeof(numeric_t)*NJ->maxnodes);
  NJ->selfdist = myfree(NJ->selfdist, sizeof(numeric_t)*NJ->maxnodes);
  NJ->selfweight = myfree(NJ->selfweight, sizeof(numeric_t)*NJ->maxnodes);
  NJ->outDistances = myfree(NJ->outDistances, sizeof(numeric_t)*NJ->maxnodes);
  NJ->nOutDistActive = myfree(NJ->nOutDistActive, sizeof(int)*NJ->maxnodes);
  NJ->parent = myfree(NJ->parent, sizeof(int)*NJ->maxnodes);
  NJ->branchlength = myfree(NJ->branchlength, sizeof(numeric_t)*NJ->maxnodes);
  NJ->support = myfree(NJ->support, sizeof(numeric_t)*NJ->maxnodes);
  NJ->child = myfree(NJ->child, sizeof(children_t)*NJ->maxnodes);
  NJ->transmat = myfree(NJ->transmat, sizeof(transition_matrix_t));
  AllocRateCategories(&NJ->rates, 0, NJ->nPos);
  return(myfree(NJ, sizeof(NJ_t)));
}