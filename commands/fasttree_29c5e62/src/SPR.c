#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void SPR(/*IN/OUT*/NJ_t *NJ, int maxSPRLength, int iRound, int nRounds) {
  /* Given a non-root node N with children A,B, sibling C, and uncle D,
     we can try to move A by doing three types of moves (4 choices):
     "down" -- swap A with a child of B (if B is not a leaf) [2 choices]
     "over" -- swap B with C
     "up" -- swap A with D
     We follow down moves with down moves, over moves with down moves, and
     up moves with either up or over moves. (Other choices are just backing
     up and hence useless.)

     As with NNIs, we keep track of up-profiles as we go. However, some of the regular
     profiles may also become "stale" so it is a bit trickier.

     We store the traversal before we do SPRs to avoid any possible infinite loop
  */
  double last_tot_len = 0.0;
  if (NJ->nSeq <= 3 || maxSPRLength < 1)
    return;
  if (slow)
    last_tot_len = TreeLength(NJ, /*recomputeLengths*/true);
  int *nodeList = mymalloc(sizeof(int) * NJ->maxnodes);
  int nodeListLen = 0;
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    nodeList[nodeListLen++] = node;
  }
  assert(nodeListLen == NJ->maxnode);
  traversal = FreeTraversal(traversal,NJ);

  profile_t **upProfiles = UpProfiles(NJ);
  spr_step_t *steps = mymalloc(sizeof(spr_step_t) * maxSPRLength); /* current chain of SPRs */

  int i;
  for (i = 0; i < nodeListLen; i++) {
    node = nodeList[i];
    if ((i % 100) == 0)
      ProgressReport("SPR round %3d of %3d, %d of %d nodes",
		     iRound+1, nRounds, i+1, nodeListLen);
    if (node == NJ->root)
      continue; /* nothing to do for root */
    /* The nodes to NNI around */
    int nodeAround[2] = { NJ->parent[node], Sibling(NJ, node) };
    if (NJ->parent[node] == NJ->root) {
      /* NNI around both siblings instead */
      RootSiblings(NJ, node, /*OUT*/nodeAround);
    }
    bool bChanged = false;
    int iAround;
    for (iAround = 0; iAround < 2 && bChanged == false; iAround++) {
      int ACFirst;
      for (ACFirst = 0; ACFirst < 2 && bChanged == false; ACFirst++) {
	if(verbose > 3)
	  PrintNJInternal(stderr, NJ, /*useLen*/false);
	int chainLength = FindSPRSteps(/*IN/OUT*/NJ, node, nodeAround[iAround],
				       upProfiles, /*OUT*/steps, maxSPRLength, (bool)ACFirst);
	double dMinDelta = 0.0;
	int iCBest = -1;
	double dTotDelta = 0.0;
	int iC;
	for (iC = 0; iC < chainLength; iC++) {
	  dTotDelta += steps[iC].deltaLength;
	  if (dTotDelta < dMinDelta) {
	    dMinDelta = dTotDelta;
	    iCBest = iC;
	  }
	}
      
	if (verbose>3) {
	  fprintf(stderr, "SPR %s %d around %d chainLength %d of %d deltaLength %.5f swaps:",
		  iCBest >= 0 ? "move" : "abandoned",
		  node,nodeAround[iAround],iCBest+1,chainLength,dMinDelta);
	  for (iC = 0; iC < chainLength; iC++)
	    fprintf(stderr, " (%d,%d)%.4f", steps[iC].nodes[0], steps[iC].nodes[1], steps[iC].deltaLength);
	  fprintf(stderr,"\n");
	}
	for (iC = chainLength - 1; iC > iCBest; iC--)
	  UnwindSPRStep(/*IN/OUT*/NJ, /*IN*/&steps[iC], /*IN/OUT*/upProfiles);
	if(verbose > 3)
	  PrintNJInternal(stderr, NJ, /*useLen*/false);
	while (slow && iCBest >= 0) {
	  double expected_tot_len = last_tot_len + dMinDelta;
	  double new_tot_len = TreeLength(NJ, /*recompute*/true);
	  if (verbose > 2)
	    fprintf(stderr, "Total branch-length is now %.4f was %.4f expected %.4f\n",
		    new_tot_len, last_tot_len, expected_tot_len);
	  if (new_tot_len < last_tot_len) {
	    last_tot_len = new_tot_len;
	    break;		/* no rewinding necessary */
	  }
	  if (verbose > 2)
	    fprintf(stderr, "Rewinding SPR to %d\n",iCBest);
	  UnwindSPRStep(/*IN/OUT*/NJ, /*IN*/&steps[iCBest], /*IN/OUT*/upProfiles);
	  dMinDelta -= steps[iCBest].deltaLength;
	  iCBest--;
	}
	if (iCBest >= 0)
	  bChanged = true;
      }	/* loop over which step to take at 1st NNI */
    } /* loop over which node to pivot around */

    if (bChanged) {
      nSPR++;		/* the SPR move is OK */
      /* make sure all the profiles are OK */
      int j;
      for (j = 0; j < NJ->maxnodes; j++)
	DeleteUpProfile(upProfiles, NJ, j);
      int ancestor;
      for (ancestor = NJ->parent[node]; ancestor >= 0; ancestor = NJ->parent[ancestor])
	RecomputeProfile(/*IN/OUT*/NJ, upProfiles, ancestor, /*useML*/false);
    }
  } /* end loop over subtrees to prune & regraft */
  steps = myfree(steps, sizeof(spr_step_t) * maxSPRLength);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  nodeList = myfree(nodeList, sizeof(int) * NJ->maxnodes);
}