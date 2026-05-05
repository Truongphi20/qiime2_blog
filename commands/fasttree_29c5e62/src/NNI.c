#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

int NNI(/*IN/OUT*/NJ_t *NJ, int iRound, int nRounds, bool useML,
	/*IN/OUT*/nni_stats_t *stats,
	/*OUT*/double *dMaxDelta) {
  /* For each non-root node N, with children A,B, sibling C, and uncle D,
     we compare the current topology AB|CD to the alternate topologies
     AC|BD and AD|BC, by using the 4 relevant profiles.

     If useML is true, it uses quartet maximum likelihood, and it
     updates branch lengths as it goes.

     If useML is false, it uses the minimum-evolution criterion with
     log-corrected distances on profiles.  (If logdist is false, then
     the log correction is not done.) If useML is false, then NNI()
     does NOT modify the branch lengths.

     Regardless of whether it changes the topology, it recomputes the
     profile for the node, using the pairwise distances and BIONJ-like
     weightings (if bionj is set). The parent's profile has changed,
     but recomputing it is not necessary because we will visit it
     before we need it (we use postorder, so we may visit the sibling
     and its children before we visit the parent, but we never
     consider an ancestor's profile, so that is OK). When we change
     the parent's profile, this alters the uncle's up-profile, so we
     remove that.  Finally, if the topology has changed, we remove the
     up-profiles of the nodes.

     If we do an NNI during post-order traversal, the result is a bit
     tricky. E.g. if we are at node N, and have visited its children A
     and B but not its uncle C, and we do an NNI that swaps B & C,
     then the post-order traversal will visit C, and its children, but
     then on the way back up, it will skip N, as it has already
     visited it.  So, the profile of N will not be recomputed: any
     changes beneath C will not be reflected in the profile of N, and
     the profile of N will be slightly stale. This will be corrected
     on the next round of NNIs.
  */
  double supportThreshold = useML ? treeLogLkDelta : MEMinDelta;
  int i;
  *dMaxDelta = 0.0;
  int nNNIThisRound = 0;

  if (NJ->nSeq <= 3)
    return(0);			/* nothing to do */
  if (verbose > 2) {
    fprintf(stderr, "Beginning round %d of NNIs with ml? %d\n", iRound, useML?1:0);
    PrintNJInternal(/*WRITE*/stderr, NJ, /*useLen*/useML && iRound > 0 ? 1 : 0);
  }
  /* For each node the upProfile or NULL */
  profile_t **upProfiles = UpProfiles(NJ);

  traversal_t traversal = InitTraversal(NJ);

  /* Identify nodes we can skip traversing into */
  int node;
  if (fastNNI) {
    for (node = 0; node < NJ->maxnode; node++) {
      if (node != NJ->root
	  && node >= NJ->nSeq
	  && stats[node].age >= 2
	  && stats[node].subtreeAge >= 2
	  && stats[node].support > supportThreshold) {
	int nodeABCD[4];
	SetupABCD(NJ, node, NULL, NULL, /*OUT*/nodeABCD, useML);
	for (i = 0; i < 4; i++)
	  if (stats[nodeABCD[i]].age == 0 && stats[nodeABCD[i]].support > supportThreshold)
	    break;
	if (i == 4) {
	  SkipTraversalInto(node, /*IN/OUT*/traversal);
	  if (verbose > 2)
	    fprintf(stderr, "Skipping subtree at %d: child %d %d parent %d age %d subtreeAge %d support %.3f\n",
		    node, nodeABCD[0], nodeABCD[1], NJ->parent[node],
		    stats[node].age, stats[node].subtreeAge, stats[node].support);
	}
      }
    }
  }

  int iDone = 0;
  bool bUp;
  node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, &bUp)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */
    if (bUp) {
      if(verbose > 2)
	fprintf(stderr, "Going up back to node %d\n", node);
      /* No longer needed */
      for (i = 0; i < NJ->child[node].nChild; i++)
	DeleteUpProfile(upProfiles, NJ, NJ->child[node].child[i]);
      DeleteUpProfile(upProfiles, NJ, node);
      RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node, useML);
      continue;
    }
    if ((iDone % 100) == 0) {
      char buf[100];
      sprintf(buf, "%s NNI round %%d of %%d, %%d of %%d splits", useML ? "ML" : "ME");
      if (iDone > 0)
	sprintf(buf+strlen(buf), ", %d changes", nNNIThisRound);
      if (nNNIThisRound > 0)
	sprintf(buf+strlen(buf), " (max delta %.3f)", *dMaxDelta);
      ProgressReport(buf, iRound+1, nRounds, iDone+1, NJ->maxnode - NJ->nSeq);
    }
    iDone++;

    profile_t *profiles[4];
    int nodeABCD[4];
    /* Note -- during the first round of ML NNIs, we use the min-evo-based branch lengths,
       which may be suboptimal */
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, useML);

    /* Given our 4 profiles, consider doing a swap */
    int nodeA = nodeABCD[0];
    int nodeB = nodeABCD[1];
    int nodeC = nodeABCD[2];
    int nodeD = nodeABCD[3];

    nni_t choice = ABvsCD;

    if (verbose > 2)
      fprintf(stderr,"Considering NNI around %d: Swap A=%d B=%d C=%d D=up(%d) or parent %d\n",
	      node, nodeA, nodeB, nodeC, nodeD, NJ->parent[node]);
    if (verbose > 3 && useML) {
      double len[5] = { NJ->branchlength[nodeA], NJ->branchlength[nodeB], NJ->branchlength[nodeC], NJ->branchlength[nodeD],
			NJ->branchlength[node] };
      for (i=0; i < 5; i++)
	if (len[i] < MLMinBranchLength)
	  len[i] = MLMinBranchLength;
      fprintf(stderr, "Starting quartet likelihood %.3f len %.3f %.3f %.3f %.3f %.3f\n",
	      MLQuartetLogLk(profiles[0],profiles[1],profiles[2],profiles[3],NJ->nPos,NJ->transmat,&NJ->rates,len, /*site_lk*/NULL),
	      len[0], len[1], len[2], len[3], len[4]);
    }

    numeric_t newlength[5];
    double criteria[3];
    if (useML) {
      for (i = 0; i < 4; i++)
	newlength[i] = NJ->branchlength[nodeABCD[i]];
      newlength[4] = NJ->branchlength[node];
      bool bFast = mlAccuracy < 2 && stats[node].age > 0;
      choice = MLQuartetNNI(profiles, NJ->transmat, &NJ->rates, NJ->nPos, NJ->nConstraints,
			    /*OUT*/criteria, /*IN/OUT*/newlength, bFast);
    } else {
      choice = ChooseNNI(profiles, NJ->distance_matrix, NJ->nPos, NJ->nConstraints,
			 /*OUT*/criteria);
      /* invert criteria so that higher is better, as in ML case, to simplify code below */
      for (i = 0; i < 3; i++)
	criteria[i] = -criteria[i];
    }
    
    if (choice == ACvsBD) {
      /* swap B and C */
      ReplaceChild(/*IN/OUT*/NJ, node, nodeB, nodeC);
      ReplaceChild(/*IN/OUT*/NJ, NJ->parent[node], nodeC, nodeB);
    } else if (choice == ADvsBC) {
      /* swap A and C */
      ReplaceChild(/*IN/OUT*/NJ, node, nodeA, nodeC);
      ReplaceChild(/*IN/OUT*/NJ, NJ->parent[node], nodeC, nodeA);
    }
    
    if (useML) {
      /* update branch length for the internal branch, and of any
	 branches that lead to leaves, b/c those will not are not
	 the internal branch for NNI and would not otherwise be set.
      */
      if (choice == ADvsBC) {
	/* For ADvsBC, MLQuartetNNI swaps B with D, but we swap A with C */
	double length2[5] = { newlength[LEN_C], newlength[LEN_D],
			      newlength[LEN_A], newlength[LEN_B],
			      newlength[LEN_I] };
	int i;
	for (i = 0; i < 5; i++) newlength[i] = length2[i];
	/* and swap A and C */
	double tmp = newlength[LEN_A];
	newlength[LEN_A] = newlength[LEN_C];
	newlength[LEN_C] = tmp;
      } else if (choice == ACvsBD) {
	/* swap B and C */
	double tmp = newlength[LEN_B];
	newlength[LEN_B] = newlength[LEN_C];
	newlength[LEN_C] = tmp;
      }
      
      NJ->branchlength[node] = newlength[LEN_I];
      NJ->branchlength[nodeA] = newlength[LEN_A];
      NJ->branchlength[nodeB] = newlength[LEN_B];
      NJ->branchlength[nodeC] = newlength[LEN_C];
      NJ->branchlength[nodeD] = newlength[LEN_D];
    }
    
    if (verbose>2 && (choice != ABvsCD || verbose > 2))
      fprintf(stderr,"NNI around %d: Swap A=%d B=%d C=%d D=out(C) -- choose %s %s %.4f\n",
	      node, nodeA, nodeB, nodeC,
	      choice == ACvsBD ? "AC|BD" : (choice == ABvsCD ? "AB|CD" : "AD|BC"),
	      useML ? "delta-loglk" : "-deltaLen",
	      criteria[choice] - criteria[ABvsCD]);
    if(verbose >= 3 && slow && useML)
      fprintf(stderr, "Old tree lk -- %.4f\n", TreeLogLk(NJ, /*site_likelihoods*/NULL));
    
    /* update stats, *dMaxDelta, etc. */
    if (choice == ABvsCD) {
      stats[node].age++;
    } else {
      if (useML)
	nML_NNI++;
      else
	nNNI++;
      nNNIThisRound++;
      stats[node].age = 0;
      stats[nodeA].age = 0;
      stats[nodeB].age = 0;
      stats[nodeC].age = 0;
      stats[nodeD].age = 0;
    }
    stats[node].delta = criteria[choice] - criteria[ABvsCD]; /* 0 if ABvsCD */
    if (stats[node].delta > *dMaxDelta)
      *dMaxDelta = stats[node].delta;
    
    /* support is improvement of score for self over better of alternatives */
    stats[node].support = 1e20;
    for (i = 0; i < 3; i++)
      if (choice != i && criteria[choice]-criteria[i] < stats[node].support)
	stats[node].support = criteria[choice]-criteria[i];
    
    /* subtreeAge is the number of rounds since self or descendent had a significant improvement */
    if (stats[node].delta > supportThreshold)
      stats[node].subtreeAge = 0;
    else {
      stats[node].subtreeAge++;
      for (i = 0; i < 2; i++) {
	int child = NJ->child[node].child[i];
	if (stats[node].subtreeAge > stats[child].subtreeAge)
	  stats[node].subtreeAge = stats[child].subtreeAge;
      }
    }

    /* update profiles and free up unneeded up-profiles */
    if (choice == ABvsCD) {
      /* No longer needed */
      DeleteUpProfile(upProfiles, NJ, nodeA);
      DeleteUpProfile(upProfiles, NJ, nodeB);
      DeleteUpProfile(upProfiles, NJ, nodeC);
      RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node, useML);
      if(slow && useML)
	UpdateForNNI(NJ, node, upProfiles, useML);
    } else {
      UpdateForNNI(NJ, node, upProfiles, useML);
    }
    if(verbose > 2 && slow && useML) {
      /* Note we recomputed profiles back up to root already if slow */
      PrintNJInternal(/*WRITE*/stderr, NJ, /*useLen*/true);
      fprintf(stderr, "New tree lk -- %.4f\n", TreeLogLk(NJ, /*site_likelihoods*/NULL));
    }
  } /* end postorder traversal */
  traversal = FreeTraversal(traversal,NJ);
  if (verbose>=2) {
    int nUp = 0;
    for (i = 0; i < NJ->maxnodes; i++)
      if (upProfiles[i] != NULL)
	nUp++;
    fprintf(stderr, "N up profiles at end of NNI:  %d\n", nUp);
  }
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  return(nNNIThisRound);
}