#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void OptimizeAllBranchLengths(/*IN/OUT*/NJ_t *NJ) {
  if (NJ->nSeq < 2)
    return;
  if (NJ->nSeq == 2) {
    int parent = NJ->root;
    assert(NJ->child[parent].nChild==2);
    int nodes[2] = { NJ->child[parent].child[0], NJ->child[parent].child[1] };
    double length = 1.0;
    (void)MLPairOptimize(NJ->profiles[nodes[0]], NJ->profiles[nodes[1]],
			 NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/&length);
    NJ->branchlength[nodes[0]] = length/2.0;
    NJ->branchlength[nodes[1]] = length/2.0;
    return;
  };

  traversal_t traversal = InitTraversal(NJ);
  profile_t **upProfiles = UpProfiles(NJ);
  int node = NJ->root;
  int iDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    int nChild = NJ->child[node].nChild;
    if (nChild > 0) {
      if ((iDone % 100) == 0)
	ProgressReport("ML Lengths %d of %d splits", iDone+1, NJ->maxnode - NJ->nSeq, 0, 0);
      iDone++;

      /* optimize the branch lengths between self, parent, and children,
         with two iterations
      */
      assert(nChild == 2 || nChild == 3);
      int nodes[3] = { NJ->child[node].child[0],
		       NJ->child[node].child[1],
		       nChild == 3 ? NJ->child[node].child[2] : node };
      profile_t *profiles[3] = { NJ->profiles[nodes[0]],
			   NJ->profiles[nodes[1]], 
			   nChild == 3 ? NJ->profiles[nodes[2]]
			   : GetUpProfile(/*IN/OUT*/upProfiles, NJ, node, /*useML*/true) };
      int iter;
      for (iter = 0; iter < 2; iter++) {
	int i;
	for (i = 0; i < 3; i++) {
	  profile_t *pA = profiles[i];
	  int b1 = (i+1) % 3;
	  int b2 = (i+2) % 3;
	  profile_t *pB = PosteriorProfile(profiles[b1], profiles[b2],
					   NJ->branchlength[nodes[b1]],
					   NJ->branchlength[nodes[b2]],
					   NJ->transmat, &NJ->rates, NJ->nPos, /*nConstraints*/0);
	  double len = NJ->branchlength[nodes[i]];
	  if (len < MLMinBranchLength)
	    len = MLMinBranchLength;
	  (void)MLPairOptimize(pA, pB, NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/&len);
	  NJ->branchlength[nodes[i]] = len;
	  pB = FreeProfile(pB, NJ->nPos, /*nConstraints*/0);
	  if (verbose>3)
	    fprintf(stderr, "Optimize length for %d to %.3f\n",
		    nodes[i], NJ->branchlength[nodes[i]]);
	}
      }
      if (node != NJ->root) {
	RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node, /*useML*/true);
	DeleteUpProfile(upProfiles, NJ, node);
      }
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}