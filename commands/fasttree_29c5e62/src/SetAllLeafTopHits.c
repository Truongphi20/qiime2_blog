#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Helper function for sorting in SetAllLeafTopHits,
   and the global variables it needs
*/
NJ_t *CompareSeedNJ = NULL;
int *CompareSeedGaps = NULL;
int CompareSeeds(const void *c1, const void *c2) {
  int seed1 = *(int *)c1;
  int seed2 = *(int *)c2;
  int gapdiff = CompareSeedGaps[seed1] - CompareSeedGaps[seed2];
  if (gapdiff != 0) return(gapdiff);	/* fewer gaps is better */
  double outdiff = CompareSeedNJ->outDistances[seed1] - CompareSeedNJ->outDistances[seed2];
  if(outdiff < 0) return(-1);	/* closer to more nodes is better */
  if(outdiff > 0) return(1);
  return(0);
}


/* Using the seed heuristic and the close global variable */
void SetAllLeafTopHits(/*IN/UPDATE*/NJ_t *NJ, /*IN/OUT*/top_hits_t *tophits) {
  double close = tophitsClose;
  if (close < 0) {
    if (fastest && NJ->nSeq >= 50000) {
      close = 0.99;
    } else {
      double logN = log((double)NJ->nSeq)/log(2.0);
      close = logN/(logN+2.0);
    }
  }
  /* Sort the potential seeds, by a combination of nGaps and NJ->outDistances
     We don't store nGaps so we need to compute that
  */
  int *nGaps = (int*)mymalloc(sizeof(int)*NJ->nSeq);
  int iNode;
  for(iNode=0; iNode<NJ->nSeq; iNode++) {
    nGaps[iNode] = (int)(0.5 + NJ->nPos - NJ->selfweight[iNode]);
  }
  int *seeds = (int*)mymalloc(sizeof(int)*NJ->nSeq);
  for (iNode=0; iNode<NJ->nSeq; iNode++) seeds[iNode] = iNode;
  CompareSeedNJ = NJ;
  CompareSeedGaps = nGaps;
  qsort(/*IN/OUT*/seeds, NJ->nSeq, sizeof(int), CompareSeeds);
  CompareSeedNJ = NULL;
  CompareSeedGaps = NULL;

  /* For each seed, save its top 2*m hits and then look for close neighbors */
  assert(2 * tophits->m <= NJ->nSeq);
  int iSeed;
  int nHasTopHits = 0;
#ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 50)
#endif
  for(iSeed=0; iSeed < NJ->nSeq; iSeed++) {
    int seed = seeds[iSeed];
    if (iSeed > 0 && (iSeed % 100) == 0) {
#ifdef OPENMP
      #pragma omp critical
#endif
      ProgressReport("Top hits for %6d of %6d seqs (at seed %6d)",
		     nHasTopHits, NJ->nSeq,
		     iSeed, 0);
    }
    if (tophits->top_hits_lists[seed].nHits > 0) {
      if(verbose>2) fprintf(stderr, "Skipping seed %d\n", seed);
      continue;
    }

    besthit_t *besthitsSeed = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->nSeq);
    besthit_t *besthitsNeighbor = (besthit_t*)mymalloc(sizeof(besthit_t) * 2 * tophits->m);
    besthit_t bestjoin;

    if(verbose>2) fprintf(stderr,"Trying seed %d\n", seed);
    SetBestHit(seed, NJ, /*nActive*/NJ->nSeq, /*OUT*/&bestjoin, /*OUT*/besthitsSeed);

    /* sort & save top hits of self. besthitsSeed is now sorted. */
    SortSaveBestHits(seed, /*IN/SORT*/besthitsSeed, /*IN-SIZE*/NJ->nSeq,
		     /*OUT-SIZE*/tophits->m, /*IN/OUT*/tophits);
    nHasTopHits++;

    /* find "close" neighbors and compute their top hits */
    double neardist = besthitsSeed[2 * tophits->m - 1].dist * close;
    /* must have at least average weight, rem higher is better
       and allow a bit more than average, e.g. if we are looking for within 30% away,
       20% more gaps than usual seems OK
       Alternatively, have a coverage requirement in case neighbor is short
       If fastest, consider the top q/2 hits to be close neighbors, regardless
    */
    double nearweight = 0;
    int iClose;
    for (iClose = 0; iClose < 2 * tophits->m; iClose++)
      nearweight += besthitsSeed[iClose].weight;
    nearweight = nearweight/(2.0 * tophits->m); /* average */
    nearweight *= (1.0-2.0*neardist/3.0);
    double nearcover = 1.0 - neardist/2.0;

    if(verbose>2) fprintf(stderr,"Distance limit for close neighbors %f weight %f ungapped %d\n",
			  neardist, nearweight, NJ->nPos-nGaps[seed]);
    for (iClose = 0; iClose < tophits->m; iClose++) {
      besthit_t *closehit = &besthitsSeed[iClose];
      int closeNode = closehit->j;
      if (tophits->top_hits_lists[closeNode].nHits > 0)
	continue;

      /* If within close-distance, or identical, use as close neighbor */
      bool close = closehit->dist <= neardist
	&& (closehit->weight >= nearweight
	    || closehit->weight >= (NJ->nPos-nGaps[closeNode])*nearcover);
      bool identical = closehit->dist < 1e-6
	&& fabs(closehit->weight - (NJ->nPos - nGaps[seed])) < 1e-5
	&& fabs(closehit->weight - (NJ->nPos - nGaps[closeNode])) < 1e-5;
      if (useTopHits2nd && iClose < tophits->q && (close || identical)) {
	nHasTopHits++;
	nClose2Used++;
	int nUse = MIN(tophits->q * tophits2Safety, 2 * tophits->m);
	besthit_t *besthitsClose = mymalloc(sizeof(besthit_t) * nUse);
	TransferBestHits(NJ, /*nActive*/NJ->nSeq,
			 closeNode,
			 /*IN*/besthitsSeed, /*SIZE*/nUse,
			 /*OUT*/besthitsClose,
			 /*updateDistance*/true);
	SortSaveBestHits(closeNode, /*IN/SORT*/besthitsClose,
			 /*IN-SIZE*/nUse, /*OUT-SIZE*/tophits->q,
			 /*IN/OUT*/tophits);
	tophits->top_hits_lists[closeNode].hitSource = seed;
	besthitsClose = myfree(besthitsClose, sizeof(besthit_t) * nUse);
      } else if (close || identical || (fastest && iClose < (tophits->q+1)/2)) {
	nHasTopHits++;
	nCloseUsed++;
	if(verbose>2) fprintf(stderr, "Near neighbor %d (rank %d weight %f ungapped %d %d)\n",
			      closeNode, iClose, besthitsSeed[iClose].weight,
			      NJ->nPos-nGaps[seed],
			      NJ->nPos-nGaps[closeNode]);

	/* compute top 2*m hits */
	TransferBestHits(NJ, /*nActive*/NJ->nSeq,
			 closeNode,
			 /*IN*/besthitsSeed, /*SIZE*/2 * tophits->m,
			 /*OUT*/besthitsNeighbor,
			 /*updateDistance*/true);
	SortSaveBestHits(closeNode, /*IN/SORT*/besthitsNeighbor,
			 /*IN-SIZE*/2 * tophits->m, /*OUT-SIZE*/tophits->m,
			 /*IN/OUT*/tophits);

	/* And then try for a second level of transfer. We assume we
	   are in a good area, because of the 1st
	   level of transfer, and in a small neighborhood, because q is
	   small (32 for 1 million sequences), so we do not make any close checks.
	 */
	int iClose2;
	for (iClose2 = 0; iClose2 < tophits->q && iClose2 < 2 * tophits->m; iClose2++) {
	  int closeNode2 = besthitsNeighbor[iClose2].j;
	  assert(closeNode2 >= 0);
	  if (tophits->top_hits_lists[closeNode2].hits == NULL) {
	    nClose2Used++;
	    nHasTopHits++;
	    int nUse = MIN(tophits->q * tophits2Safety, 2 * tophits->m);
	    besthit_t *besthitsClose2 = mymalloc(sizeof(besthit_t) * nUse);
	    TransferBestHits(NJ, /*nActive*/NJ->nSeq,
			     closeNode2,
			     /*IN*/besthitsNeighbor, /*SIZE*/nUse,
			     /*OUT*/besthitsClose2,
			     /*updateDistance*/true);
	    SortSaveBestHits(closeNode2, /*IN/SORT*/besthitsClose2,
			     /*IN-SIZE*/nUse, /*OUT-SIZE*/tophits->q,
			     /*IN/OUT*/tophits);
	    tophits->top_hits_lists[closeNode2].hitSource = closeNode;
	    besthitsClose2 = myfree(besthitsClose2, sizeof(besthit_t) * nUse);
	  } /* end if should do 2nd-level transfer */
	}
      }
    } /* end loop over close candidates */
    besthitsSeed = myfree(besthitsSeed, sizeof(besthit_t)*NJ->nSeq);
    besthitsNeighbor = myfree(besthitsNeighbor, sizeof(besthit_t) * 2 * tophits->m);
  } /* end loop over seeds */

  for (iNode=0; iNode<NJ->nSeq; iNode++) {
    top_hits_list_t *l = &tophits->top_hits_lists[iNode];
    assert(l->hits != NULL);
    assert(l->hits[0].j >= 0);
    assert(l->hits[0].j < NJ->nSeq);
    assert(l->hits[0].j != iNode);
    tophits->visible[iNode] = l->hits[0];
  }

  if (verbose >= 2) fprintf(stderr, "#Close neighbors among leaves: 1st-level %ld 2nd-level %ld seeds %ld\n",
			    nCloseUsed, nClose2Used, NJ->nSeq-nCloseUsed-nClose2Used);
  nGaps = myfree(nGaps, sizeof(int)*NJ->nSeq);
  seeds = myfree(seeds, sizeof(int)*NJ->nSeq);

  /* Now add a "checking phase" where we ensure that the q or 2*sqrt(m) hits
     of i are represented in j (if they should be)
   */
  long lReplace = 0;
  int nCheck = tophits->q > 0 ? tophits->q : (int)(0.5 + 2.0*sqrt(tophits->m));
  for (iNode = 0; iNode < NJ->nSeq; iNode++) {
    if ((iNode % 100) == 0)
      ProgressReport("Checking top hits for %6d of %6d seqs",
		     iNode+1, NJ->nSeq, 0, 0);
    top_hits_list_t *lNode = &tophits->top_hits_lists[iNode];
    int iHit;
    for (iHit = 0; iHit < nCheck && iHit < lNode->nHits; iHit++) {
      besthit_t bh = HitToBestHit(iNode, lNode->hits[iHit]);
      SetCriterion(NJ, /*nActive*/NJ->nSeq, /*IN/OUT*/&bh);
      top_hits_list_t *lTarget = &tophits->top_hits_lists[bh.j];

      /* If this criterion is worse than the nCheck-1 entry of the target,
	 then skip the check.
	 This logic is based on assuming that the list is sorted,
	 which is true initially but may not be true later.
	 Still, is a good heuristic.
      */
      assert(nCheck > 0);
      assert(nCheck <= lTarget->nHits);
      besthit_t bhCheck = HitToBestHit(bh.j, lTarget->hits[nCheck-1]);
      SetCriterion(NJ, /*nActive*/NJ->nSeq, /*IN/OUT*/&bhCheck);
      if (bhCheck.criterion < bh.criterion)
	continue;		/* no check needed */

      /* Check if this is present in the top-hit list */
      int iHit2;
      bool bFound = false;
      for (iHit2 = 0; iHit2 < lTarget->nHits && !bFound; iHit2++)
	if (lTarget->hits[iHit2].j == iNode)
	  bFound = true;
      if (!bFound) {
	/* Find the hit with the worst criterion and replace it with this one */
	int iWorst = -1;
	double dWorstCriterion = -1e20;
	for (iHit2 = 0; iHit2 < lTarget->nHits; iHit2++) {
	  besthit_t bh2 = HitToBestHit(bh.j, lTarget->hits[iHit2]);
	  SetCriterion(NJ, /*nActive*/NJ->nSeq, /*IN/OUT*/&bh2);
	  if (bh2.criterion > dWorstCriterion) {
	    iWorst = iHit2;
	    dWorstCriterion = bh2.criterion;
	  }
	}
	if (dWorstCriterion > bh.criterion) {
	  assert(iWorst >= 0);
	  lTarget->hits[iWorst].j = iNode;
	  lTarget->hits[iWorst].dist = bh.dist;
	  lReplace++;
	  /* and perhaps update visible */
	  besthit_t v;
	  bool bSuccess = GetVisible(NJ, /*nActive*/NJ->nSeq, tophits, bh.j, /*OUT*/&v);
	  assert(bSuccess);
	  if (bh.criterion < v.criterion)
	    tophits->visible[bh.j] = lTarget->hits[iWorst];
	}
      }
    }
  }

  if (verbose >= 2)
    fprintf(stderr, "Replaced %ld top hit entries\n", lReplace);
}