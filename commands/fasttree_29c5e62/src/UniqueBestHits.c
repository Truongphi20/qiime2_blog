#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

besthit_t *UniqueBestHits(/*IN/UPDATE*/NJ_t *NJ, int nActive,
			  /*IN/SORT*/besthit_t *combined, int nCombined,
			  /*OUT*/int *nUniqueOut) {
  int iHit;
  for (iHit = 0; iHit < nCombined; iHit++) {
    besthit_t *hit = &combined[iHit];
    UpdateBestHit(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/hit, /*update*/false);
  }
  qsort(/*IN/OUT*/combined, nCombined, sizeof(besthit_t), CompareHitsByIJ);

  besthit_t *uniqueList = (besthit_t*)mymalloc(sizeof(besthit_t)*nCombined);
  int nUnique = 0;
  int iSavedLast = -1;

  /* First build the new list */
  for (iHit = 0; iHit < nCombined; iHit++) {
    besthit_t *hit = &combined[iHit];
    if (hit->i < 0 || hit->j < 0)
      continue;
    if (iSavedLast >= 0) {
      /* toss out duplicates */
      besthit_t *saved = &combined[iSavedLast];
      if (saved->i == hit->i && saved->j == hit->j)
	continue;
    }
    assert(nUnique < nCombined);
    assert(hit->j >= 0 && NJ->parent[hit->j] < 0);
    uniqueList[nUnique++] = *hit;
    iSavedLast = iHit;
  }
  *nUniqueOut = nUnique;

  /* Then do any updates to the criterion or the distances in parallel */
#ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 50)
#endif
  for (iHit = 0; iHit < nUnique; iHit++) {
    besthit_t *hit = &uniqueList[iHit];
    if (hit->dist < 0.0)
      SetDistCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/hit);
    else
      SetCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/hit);
  }
  return(uniqueList);
}