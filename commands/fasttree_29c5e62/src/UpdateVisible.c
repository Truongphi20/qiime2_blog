#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void UpdateVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		   /*IN*/besthit_t *tophitsNode,
		   int nTopHits,
		  /*IN/OUT*/top_hits_t *tophits) {
  int iHit;

  for(iHit = 0; iHit < nTopHits; iHit++) {
    besthit_t *hit = &tophitsNode[iHit];
    if (hit->i < 0) continue;	/* possible empty entries */
    assert(NJ->parent[hit->i] < 0);
    assert(hit->j >= 0 && NJ->parent[hit->j] < 0);
    besthit_t visible;
    bool bSuccess = GetVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, hit->j, /*OUT*/&visible);
    if (!bSuccess || hit->criterion < visible.criterion) {
      if (bSuccess)
	nVisibleUpdate++;
      hit_t *v = &tophits->visible[hit->j];
      v->j = hit->i;
      v->dist = hit->dist;
      UpdateTopVisible(NJ, nActive, hit->j, v, /*IN/OUT*/tophits);
      if(verbose>5) fprintf(stderr,"NewVisible %d %d %f\n",
			    hit->j,v->j,v->dist);
    }
  } /* end loop over hits */
}