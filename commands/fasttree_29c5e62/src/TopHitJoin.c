#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/*
  Create a top hit list for the new node, either
  from children (if there are enough best hits left) or by a "refresh"
  Also set visible set for newnode
  Also update visible set for other nodes if we stumble across a "better" hit
*/
 
void TopHitJoin(int newnode,
		/*IN/UPDATE*/NJ_t *NJ,
		int nActive,
		/*IN/OUT*/top_hits_t *tophits) {
  long startProfileOps = profileOps;
  long startOutProfileOps = outprofileOps;
  assert(NJ->child[newnode].nChild == 2);
  top_hits_list_t *lNew = &tophits->top_hits_lists[newnode];
  assert(lNew->hits == NULL);

  /* Copy the hits */
  int i;
  top_hits_list_t *lChild[2];
  for (i = 0; i< 2; i++) {
    lChild[i] = &tophits->top_hits_lists[NJ->child[newnode].child[i]];
    assert(lChild[i]->hits != NULL && lChild[i]->nHits > 0);
  }
  int nCombined = lChild[0]->nHits + lChild[1]->nHits;
  besthit_t *combinedList = (besthit_t*)mymalloc(sizeof(besthit_t)*nCombined);
  HitsToBestHits(lChild[0]->hits, lChild[0]->nHits, NJ->child[newnode].child[0],
		 /*OUT*/combinedList);
  HitsToBestHits(lChild[1]->hits, lChild[1]->nHits, NJ->child[newnode].child[1],
		 /*OUT*/combinedList + lChild[0]->nHits);
  int nUnique;
  /* UniqueBestHits() replaces children (used in the calls to HitsToBestHits)
     with active ancestors, so all distances & criteria will be recomputed */
  besthit_t *uniqueList = UniqueBestHits(/*IN/UPDATE*/NJ, nActive,
					 /*IN/SORT*/combinedList,
					 nCombined,
					 /*OUT*/&nUnique);
  int nUniqueAlloc = nCombined;
  combinedList = myfree(combinedList, sizeof(besthit_t)*nCombined);

  /* Forget the top-hit lists of the joined nodes */
  for (i = 0; i < 2; i++) {
    lChild[i]->hits = myfree(lChild[i]->hits, sizeof(hit_t) * lChild[i]->nHits);
    lChild[i]->nHits = 0;
  }

  /* Use the average age, rounded up, by 1 Versions 2.0 and earlier
     used the maximum age, which leads to more refreshes without
     improving the accuracy of the NJ phase. Intuitively, if one of
     them was just refreshed then another refresh is unlikely to help.
   */
  lNew->age = (lChild[0]->age+lChild[1]->age+1)/2 + 1;

  /* If top hit ages always match (perfectly balanced), then a
     limit of log2(m) would mean a refresh after
     m joins, which is about what we want.
  */
  int tophitAgeLimit = MAX(1, (int)(0.5 + log((double)tophits->m)/log(2.0)));

  /* Either use the merged list as candidate top hits, or
     move from 2nd level to 1st level, or do a refresh
     UniqueBestHits eliminates hits to self, so if nUnique==nActive-1,
     we've already done the exhaustive search.

     Either way, we set tophits, visible(newnode), update visible of its top hits,
     and modify topvisible: if we do a refresh, then we reset it, otherwise we update
  */
  bool bSecondLevel = lChild[0]->hitSource >= 0 && lChild[1]->hitSource >= 0;
  bool bUseUnique = nUnique==nActive-1
    || (lNew->age <= tophitAgeLimit
	&& nUnique >= (bSecondLevel ? (int)(0.5 + tophits2Refresh * tophits->q)
		       : (int)(0.5 + tophits->m * tophitsRefresh) ));
  if (bUseUnique && verbose > 2)
    fprintf(stderr,"Top hits for %d from combined %d nActive=%d tophitsage %d %s\n",
	    newnode,nUnique,nActive,lNew->age,
	    bSecondLevel ? "2ndlevel" : "1stlevel");

  if (!bUseUnique
      && bSecondLevel
      && lNew->age <= tophitAgeLimit) {
    int source = ActiveAncestor(NJ, lChild[0]->hitSource);
    if (source == newnode)
      source = ActiveAncestor(NJ, lChild[1]->hitSource);
    /* In parallel mode, it is possible that we would select a node as the
       hit-source and then over-write that top hit with a short list.
       So we need this sanity check.
    */
    if (source != newnode
	&& source >= 0
	&& tophits->top_hits_lists[source].hitSource < 0) {

      /* switch from 2nd-level to 1st-level top hits -- compute top hits list
	 of node from what we have so far plus the active source plus its top hits */
      top_hits_list_t *lSource = &tophits->top_hits_lists[source];
      assert(lSource->hitSource < 0);
      assert(lSource->nHits > 0);
      int nMerge = 1 + lSource->nHits + nUnique;
      besthit_t *mergeList = mymalloc(sizeof(besthit_t) * nMerge);
      memcpy(/*to*/mergeList, /*from*/uniqueList, nUnique * sizeof(besthit_t));
      
      int iMerge = nUnique;
      mergeList[iMerge].i = newnode;
      mergeList[iMerge].j = source;
      SetDistCriterion(NJ, nActive, /*IN/OUT*/&mergeList[iMerge]);
      iMerge++;
      HitsToBestHits(lSource->hits, lSource->nHits, newnode, /*OUT*/mergeList+iMerge);
      for (i = 0; i < lSource->nHits; i++) {
	SetDistCriterion(NJ, nActive, /*IN/OUT*/&mergeList[iMerge]);
	iMerge++;
      }
      assert(iMerge == nMerge);
      
      uniqueList = myfree(uniqueList, nUniqueAlloc * sizeof(besthit_t));
      uniqueList = UniqueBestHits(/*IN/UPDATE*/NJ, nActive,
				  /*IN/SORT*/mergeList,
				  nMerge,
				  /*OUT*/&nUnique);
      nUniqueAlloc = nMerge;
      mergeList = myfree(mergeList, sizeof(besthit_t)*nMerge);
      
      assert(nUnique > 0);
      bUseUnique = nUnique >= (int)(0.5 + tophits->m * tophitsRefresh);
      bSecondLevel = false;
      
      if (bUseUnique && verbose > 2)
	fprintf(stderr, "Top hits for %d from children and source %d's %d hits, nUnique %d\n",
		newnode, source, lSource->nHits, nUnique);
    }
  }

  if (bUseUnique) {
    if (bSecondLevel) {
      /* pick arbitrarily */
      lNew->hitSource = lChild[0]->hitSource;
    }
    int nSave = MIN(nUnique, bSecondLevel ? tophits->q : tophits->m);
    assert(nSave>0);
    if (verbose > 2)
      fprintf(stderr, "Combined %d ops so far %ld\n", nUnique, profileOps - startProfileOps);
    SortSaveBestHits(newnode, /*IN/SORT*/uniqueList, /*nIn*/nUnique,
		     /*nOut*/nSave, /*IN/OUT*/tophits);
    assert(lNew->hits != NULL); /* set by sort/save */
    tophits->visible[newnode] = lNew->hits[0];
    UpdateTopVisible(/*IN*/NJ, nActive, newnode, &tophits->visible[newnode],
		     /*IN/OUT*/tophits);
    UpdateVisible(/*IN/UPDATE*/NJ, nActive, /*IN*/uniqueList, nSave, /*IN/OUT*/tophits);
  } else {
    /* need to refresh: set top hits for node and for its top hits */
    if(verbose > 2) fprintf(stderr,"Top hits for %d by refresh (%d unique age %d) nActive=%d\n",
			  newnode,nUnique,lNew->age,nActive);
    nRefreshTopHits++;
    lNew->age = 0;

    int iNode;
    /* ensure all out-distances are up to date ahead of time
       to avoid any data overwriting issues.
    */
#ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 50)
#endif
    for (iNode = 0; iNode < NJ->maxnode; iNode++) {
      if (NJ->parent[iNode] < 0) {
	if (fastest) {
	  besthit_t bh;
	  bh.i = iNode;
	  bh.j = iNode;
	  bh.dist = 0;
	  SetCriterion(/*IN/UPDATE*/NJ, nActive, &bh);
	} else {
	  SetOutDistance(/*IN/UDPATE*/NJ, iNode, nActive);
	}
      }
    }

    /* exhaustively get the best 2*m hits for newnode, set visible, and save the top m */
    besthit_t *allhits = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnode);
    assert(2 * tophits->m <= NJ->maxnode);
    besthit_t bh;
    SetBestHit(newnode, NJ, nActive, /*OUT*/&bh, /*OUT*/allhits);
    qsort(/*IN/OUT*/allhits, NJ->maxnode, sizeof(besthit_t), CompareHitsByCriterion);
    SortSaveBestHits(newnode, /*IN/SORT*/allhits, /*nIn*/NJ->maxnode,
		     /*nOut*/tophits->m, /*IN/OUT*/tophits);

    /* Do not need to call UpdateVisible because we set visible below */

    /* And use the top 2*m entries to expand other best-hit lists, but only for top m */
    int iHit;
#ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 50)
#endif
    for (iHit=0; iHit < tophits->m; iHit++) {
      if (allhits[iHit].i < 0) continue;
      int iNode = allhits[iHit].j;
      assert(iNode>=0);
      if (NJ->parent[iNode] >= 0) continue;
      top_hits_list_t *l = &tophits->top_hits_lists[iNode];
      int nHitsOld = l->nHits;
      assert(nHitsOld <= tophits->m);
      l->age = 0;

      /* Merge: old hits into 0->nHitsOld and hits from iNode above that */
      besthit_t *bothList = (besthit_t*)mymalloc(sizeof(besthit_t) * 3 * tophits->m);
      HitsToBestHits(/*IN*/l->hits, nHitsOld, iNode, /*OUT*/bothList); /* does not compute criterion */
      for (i = 0; i < nHitsOld; i++)
	SetCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/&bothList[i]);
      if (nActive <= 2 * tophits->m)
	l->hitSource = -1;	/* abandon the 2nd-level top-hits heuristic */
      int nNewHits = l->hitSource >= 0 ? tophits->q : tophits->m;
      assert(nNewHits > 0);

      TransferBestHits(/*IN/UPDATE*/NJ, nActive, iNode,
		       /*IN*/allhits, /*nOldHits*/2 * nNewHits,
		       /*OUT*/&bothList[nHitsOld],
		       /*updateDist*/false); /* rely on UniqueBestHits to update dist and/or criterion */
      int nUnique2;
      besthit_t *uniqueList2 = UniqueBestHits(/*IN/UPDATE*/NJ, nActive,
					      /*IN/SORT*/bothList, nHitsOld + 2 * nNewHits,
					      /*OUT*/&nUnique2);
      assert(nUnique2 > 0);
      bothList = myfree(bothList,3 * tophits->m * sizeof(besthit_t));

      /* Note this will overwrite l, but we saved nHitsOld */
      SortSaveBestHits(iNode, /*IN/SORT*/uniqueList2, /*nIn*/nUnique2,
		       /*nOut*/nNewHits, /*IN/OUT*/tophits);
      /* will update topvisible below */
      tophits->visible[iNode] = tophits->top_hits_lists[iNode].hits[0];
      uniqueList2 = myfree(uniqueList2, (nHitsOld + 2 * tophits->m) * sizeof(besthit_t));
    }

    ResetTopVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits); /* outside of the parallel phase */
    allhits = myfree(allhits,sizeof(besthit_t)*NJ->maxnode);
  }
  uniqueList = myfree(uniqueList, nUniqueAlloc * sizeof(besthit_t));
  if (verbose > 2) {
    fprintf(stderr, "New top-hit list for %d profile-ops %ld (out-ops %ld): source %d age %d members ",
	    newnode,
	    profileOps - startProfileOps,
	    outprofileOps - startOutProfileOps,
	    lNew->hitSource, lNew->age);

    int i;
    for (i = 0; i < lNew->nHits; i++)
      fprintf(stderr, " %d", lNew->hits[i].j);
    fprintf(stderr,"\n");
  }
}