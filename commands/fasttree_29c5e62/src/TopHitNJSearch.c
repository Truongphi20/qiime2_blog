#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void TopHitNJSearch(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		    /*IN/OUT*/top_hits_t *tophits,
		    /*OUT*/besthit_t *join) {
  /* first, do we have at least m/2 candidates in topvisible?
     And remember the best one */
  int nCandidate = 0;
  int iNodeBestCandidate = -1;
  double dBestCriterion = 1e20;

  int i;
  for (i = 0; i < tophits->nTopVisible; i++) {
    int iNode = tophits->topvisible[i];
    besthit_t visible;
    if (GetVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, iNode, /*OUT*/&visible)) {
      nCandidate++;
      if (iNodeBestCandidate < 0 || visible.criterion < dBestCriterion) {
	iNodeBestCandidate = iNode;
	dBestCriterion = visible.criterion;
      }
    }
  }
  
  tophits->topvisibleAge++;
  /* Note we may have only nActive/2 joins b/c we try to store them once */
  if (2 * tophits->topvisibleAge > tophits->m
      || (3*nCandidate < tophits->nTopVisible && 3*nCandidate < nActive)) {
    /* recompute top visible */
    if (verbose > 2)
      fprintf(stderr, "Resetting the top-visible list at nActive=%d\n",nActive);

    /* If age is low, then our visible set is becoming too sparse, because we have
       recently recomputed the top visible subset. This is very rare but can happen
       with -fastest. A quick-and-dirty solution is to walk up
       the parents to get additional entries in top hit lists. To ensure that the
       visible set becomes full, pick an arbitrary node if walking up terminates at self.
    */
    if (tophits->topvisibleAge <= 2) {
      if (verbose > 2)
	fprintf(stderr, "Expanding visible set by walking up to active nodes at nActive=%d\n", nActive);
      int iNode;
      for (iNode = 0; iNode < NJ->maxnode; iNode++) {
	if (NJ->parent[iNode] >= 0)
	  continue;
	hit_t *v = &tophits->visible[iNode];
	int newj = ActiveAncestor(NJ, v->j);
	if (newj >= 0 && newj != v->j) {
	  if (newj == iNode) {
	    /* pick arbitrarily */
	    newj = 0;
	    while (NJ->parent[newj] >= 0 || newj == iNode)
	      newj++;
	  }
	  assert(newj >= 0 && newj < NJ->maxnodes
		 && newj != iNode
		 && NJ->parent[newj] < 0);

	  /* Set v to point to newj */
	  besthit_t bh = { iNode, newj, -1e20, -1e20, -1e20 };
	  SetDistCriterion(NJ, nActive, /*IN/OUT*/&bh);
	  v->j = newj;
	  v->dist = bh.dist;
	}
      }
    }
    ResetTopVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits);
    /* and recurse to try again */
    TopHitNJSearch(NJ, nActive, tophits, join);
    return;
  }
  if (verbose > 2)
    fprintf(stderr, "Top-visible list size %d (nActive %d m %d)\n",
	    nCandidate, nActive, tophits->m);
  assert(iNodeBestCandidate >= 0 && NJ->parent[iNodeBestCandidate] < 0);
  bool bSuccess = GetVisible(NJ, nActive, tophits, iNodeBestCandidate, /*OUT*/join);
  assert(bSuccess);
  assert(join->i >= 0 && NJ->parent[join->i] < 0);
  assert(join->j >= 0 && NJ->parent[join->j] < 0);

  if(fastest)
    return;

  int changed;
  do {
    changed = 0;

    besthit_t bestI;
    GetBestFromTopHits(join->i, NJ, nActive, tophits, /*OUT*/&bestI);
    assert(bestI.i == join->i);
    if (bestI.j != join->j && bestI.criterion < join->criterion) {
      changed = 1;
      if (verbose>2)
	fprintf(stderr,"BetterI\t%d\t%d\t%d\t%d\t%f\t%f\n",
		join->i,join->j,bestI.i,bestI.j,
		join->criterion,bestI.criterion);
      *join = bestI;
    }

    besthit_t bestJ;
    GetBestFromTopHits(join->j, NJ, nActive, tophits, /*OUT*/&bestJ);
    assert(bestJ.i == join->j);
    if (bestJ.j != join->i && bestJ.criterion < join->criterion) {
      changed = 1;
      if (verbose>2)
	fprintf(stderr,"BetterJ\t%d\t%d\t%d\t%d\t%f\t%f\n",
		join->i,join->j,bestJ.i,bestJ.j,
		join->criterion,bestJ.criterion);
      *join = bestJ;
    }
    if(changed) nHillBetter++;
  } while(changed);
}