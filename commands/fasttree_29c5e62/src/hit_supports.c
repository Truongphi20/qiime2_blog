#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

int CompareHitsByCriterion(const void *c1, const void *c2) {
  const besthit_t *hit1 = (besthit_t*)c1;
  const besthit_t *hit2 = (besthit_t*)c2;
  if (hit1->criterion < hit2->criterion) return(-1);
  if (hit1->criterion > hit2->criterion) return(1);
  return(0);
}

int CompareHitsByIJ(const void *c1, const void *c2) {
  const besthit_t *hit1 = (besthit_t*)c1;
  const besthit_t *hit2 = (besthit_t*)c2;
  return hit1->i != hit2->i ? hit1->i - hit2->i : hit1->j - hit2->j;
}

void SortSaveBestHits(int iNode, /*IN/SORT*/besthit_t *besthits,
		      int nIn, int nOut,
		      /*IN/OUT*/top_hits_t *tophits) {
  assert(nIn > 0);
  assert(nOut > 0);
  top_hits_list_t *l = &tophits->top_hits_lists[iNode];
  /*  */
  qsort(/*IN/OUT*/besthits,nIn,sizeof(besthit_t),CompareHitsByCriterion);

  /* First count how many we will save
     Not sure if removing duplicates is actually necessary.
   */
  int nSave = 0;
  int jLast = -1;
  int iBest;
  for (iBest = 0; iBest < nIn && nSave < nOut; iBest++) {
    if (besthits[iBest].i < 0)
      continue;
    assert(besthits[iBest].i == iNode);
    int j = besthits[iBest].j;
    if (j != iNode && j != jLast && j >= 0) {
      nSave++;
      jLast = j;
    }
  }

  assert(nSave > 0);

#ifdef OPENMP
  omp_set_lock(&tophits->locks[iNode]);
#endif
  if (l->hits != NULL) {
    l->hits = myfree(l->hits, l->nHits * sizeof(hit_t));
    l->nHits = 0;
  }
  l->hits = mymalloc(sizeof(hit_t) * nSave);
  l->nHits = nSave;
  int iSave = 0;
  jLast = -1;
  for (iBest = 0; iBest < nIn && iSave < nSave; iBest++) {
    int j = besthits[iBest].j;
    if (j != iNode && j != jLast && j >= 0) {
      l->hits[iSave].j = j;
      l->hits[iSave].dist = besthits[iBest].dist;
      iSave++;
      jLast = j;
    }
  }
#ifdef OPENMP
  omp_unset_lock(&tophits->locks[iNode]);
#endif
  assert(iSave == nSave);
}

void TransferBestHits(/*IN/UPDATE*/NJ_t *NJ,
		       int nActive,
		      int iNode,
		      /*IN*/besthit_t *oldhits,
		      int nOldHits,
		      /*OUT*/besthit_t *newhits,
		      bool updateDistances) {
  assert(iNode >= 0);
  assert(NJ->parent[iNode] < 0);

  int iBest;
  for(iBest = 0; iBest < nOldHits; iBest++) {
    besthit_t *old = &oldhits[iBest];
    besthit_t *new = &newhits[iBest];
    new->i = iNode;
    new->j = ActiveAncestor(/*IN*/NJ, old->j);
    new->dist = old->dist;	/* may get reset below */
    new->weight = old->weight;
    new->criterion = old->criterion;

    if(new->j < 0 || new->j == iNode) {
      new->weight = 0;
      new->dist = -1e20;
      new->criterion = 1e20;
    } else if (new->i != old->i || new->j != old->j) {
      if (updateDistances)
	SetDistCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/new);
      else {
	new->dist = -1e20;
	new->criterion = 1e20;
      }
    } else {
      if (updateDistances)
	SetCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/new);
      else
	new->criterion = 1e20;	/* leave dist alone */
    }
  }
}

void HitsToBestHits(/*IN*/hit_t *hits, int nHits, int iNode, /*OUT*/besthit_t *newhits) {
  int i;
  for (i = 0; i < nHits; i++) {
    hit_t *hit = &hits[i];
    besthit_t *bh = &newhits[i];
    bh->i = iNode;
    bh->j = hit->j;
    bh->dist = hit->dist;
    bh->criterion = 1e20;
    bh->weight = -1;		/* not the true value -- we compute these directly when needed */
  }
}

besthit_t HitToBestHit(int i, hit_t hit) {
  besthit_t bh;
  bh.i = i;
  bh.j = hit.j;
  bh.dist = hit.dist;
  bh.criterion = 1e20;
  bh.weight = -1;
  return(bh);
}

bool UpdateBestHit(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *hit,
		   bool bUpdateDist) {
  int i = ActiveAncestor(/*IN*/NJ, hit->i);
  int j = ActiveAncestor(/*IN*/NJ, hit->j);
  if (i < 0 || j < 0 || i == j) {
    hit->i = -1;
    hit->j = -1;
    hit->weight = 0;
    hit->dist = 1e20;
    hit->criterion = 1e20;
    return(false);
  }
  if (i != hit->i || j != hit->j) {
    hit->i = i;
    hit->j = j;
    if (bUpdateDist) {
      SetDistCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/hit);
    } else {
      hit->dist = -1e20;
      hit->criterion = 1e20;
    }
  }
  return(true);
}