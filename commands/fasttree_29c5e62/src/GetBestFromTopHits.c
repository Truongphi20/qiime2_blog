#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void GetBestFromTopHits(int iNode,
			/*IN/UPDATE*/NJ_t *NJ,
			int nActive,
			/*IN*/top_hits_t *tophits,
			/*OUT*/besthit_t *bestjoin) {
  assert(iNode >= 0);
  assert(NJ->parent[iNode] < 0);
  top_hits_list_t *l = &tophits->top_hits_lists[iNode];
  assert(l->nHits > 0);
  assert(l->hits != NULL);

  if(!fastest)
    SetOutDistance(NJ, iNode, nActive); /* ensure out-distances are not stale */

  bestjoin->i = -1;
  bestjoin->j = -1;
  bestjoin->dist = 1e20;
  bestjoin->criterion = 1e20;

  int iBest;
  for(iBest=0; iBest < l->nHits; iBest++) {
    besthit_t bh = HitToBestHit(iNode, l->hits[iBest]);
    if (UpdateBestHit(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/&bh, /*update dist*/true)) {
      SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/&bh); /* make sure criterion is correct */
      if (bh.criterion < bestjoin->criterion)
	*bestjoin = bh;
    }
  }
  assert(bestjoin->j >= 0);	/* a hit was found */
  assert(bestjoin->i == iNode);
}