#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Recompute the topvisible list */
void ResetTopVisible(/*IN/UPDATE*/NJ_t *NJ,
		     int nActive,
		     /*IN/OUT*/top_hits_t *tophits) {
  besthit_t *visibleSorted = mymalloc(sizeof(besthit_t)*nActive);
  int nVisible = 0;		/* #entries in visibleSorted */
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    /* skip joins involving stale nodes */
    if (NJ->parent[iNode] >= 0)
      continue;
    besthit_t v;
    if (GetVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, iNode, /*OUT*/&v)) {
      assert(nVisible < nActive);
      visibleSorted[nVisible++] = v;
    }
  }
  assert(nVisible > 0);
    
  qsort(/*IN/OUT*/visibleSorted,nVisible,sizeof(besthit_t),CompareHitsByCriterion);
    
  /* Only keep the top m items, and try to avoid duplicating i->j with j->i
     Note that visible(i) -> j does not necessarily imply visible(j) -> i,
     so we store what the pairing was (or -1 for not used yet)
   */
  int *inTopVisible = malloc(sizeof(int) * NJ->maxnodes);
  int i;
  for (i = 0; i < NJ->maxnodes; i++)
    inTopVisible[i] = -1;

  if (verbose > 2)
    fprintf(stderr, "top-hit search: nActive %d nVisible %d considering up to %d items\n",
	    nActive, nVisible, tophits->m);

  /* save the sorted indices in topvisible */
  int iSave = 0;
  for (i = 0; i < nVisible && iSave < tophits->nTopVisible; i++) {
    besthit_t *v = &visibleSorted[i];
    if (inTopVisible[v->i] != v->j) { /* not seen already */
      tophits->topvisible[iSave++] = v->i;
      inTopVisible[v->i] = v->j;
      inTopVisible[v->j] = v->i;
    }
  }
  while(iSave < tophits->nTopVisible)
    tophits->topvisible[iSave++] = -1;
  myfree(visibleSorted, sizeof(besthit_t)*nActive);
  myfree(inTopVisible, sizeof(int) * NJ->maxnodes);
  tophits->topvisibleAge = 0;
  if (verbose > 2) {
    fprintf(stderr, "Reset TopVisible: ");
    for (i = 0; i < tophits->nTopVisible; i++) {
      int iNode = tophits->topvisible[i];
      if (iNode < 0)
	break;
      fprintf(stderr, " %d=>%d", iNode, tophits->visible[iNode].j);
    }
    fprintf(stderr,"\n");
  }
}