#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void SetCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join) {
  if(join->i < 0
     || join->j < 0
     || NJ->parent[join->i] >= 0
     || NJ->parent[join->j] >= 0)
    return;
  assert(NJ->nOutDistActive[join->i] >= nActive);
  assert(NJ->nOutDistActive[join->j] >= nActive);

  int nDiffAllow = tophitsMult > 0 ? (int)(nActive*staleOutLimit) : 0;
  if (NJ->nOutDistActive[join->i] - nActive > nDiffAllow)
    SetOutDistance(NJ, join->i, nActive);
  if (NJ->nOutDistActive[join->j] - nActive > nDiffAllow)
    SetOutDistance(NJ, join->j, nActive);
  double outI = NJ->outDistances[join->i];
  if (NJ->nOutDistActive[join->i] != nActive)
    outI *= (nActive-1)/(double)(NJ->nOutDistActive[join->i]-1);
  double outJ = NJ->outDistances[join->j];
  if (NJ->nOutDistActive[join->j] != nActive)
    outJ *= (nActive-1)/(double)(NJ->nOutDistActive[join->j]-1);
  join->criterion = join->dist - (outI+outJ)/(double)(nActive-2);
  if (verbose > 2 && nActive <= 5) {
    fprintf(stderr, "Set Criterion to join %d %d with nActive=%d dist+penalty %.3f criterion %.3f\n",
	    join->i, join->j, nActive, join->dist, join->criterion);
  }
}