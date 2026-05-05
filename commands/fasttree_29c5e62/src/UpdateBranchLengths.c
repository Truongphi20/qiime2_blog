#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Recomputes all branch lengths

   For internal branches such as (A,B) vs. (C,D), uses the formula 

   length(AB|CD) = (d(A,C)+d(A,D)+d(B,C)+d(B,D))/4 - d(A,B)/2 - d(C,D)/2

   (where all distances are profile distances - diameters).

   For external branches (e.g. to leaves) A vs. (B,C), use the formula

   length(A|BC) = (d(A,B)+d(A,C)-d(B,C))/2
*/
void UpdateBranchLengths(/*IN/OUT*/NJ_t *NJ) {
  if (NJ->nSeq < 2)
    return;
  else if (NJ->nSeq == 2) {
    int root = NJ->root;
    int nodeA = NJ->child[root].child[0];
    int nodeB = NJ->child[root].child[1];
    besthit_t h;
    ProfileDist(NJ->profiles[nodeA],NJ->profiles[nodeB],
		NJ->nPos, NJ->distance_matrix, /*OUT*/&h);
    if (logdist)
      h.dist = LogCorrect(h.dist);
    NJ->branchlength[nodeA] = h.dist/2.0;
    NJ->branchlength[nodeB] = h.dist/2.0;
    return;
  }

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    /* reset branch length of node (distance to its parent) */
    if (node == NJ->root)
      continue; /* no branch length to set */
    if (node < NJ->nSeq) { /* a leaf */
      profile_t *profileA = NJ->profiles[node];
      profile_t *profileB = NULL;
      profile_t *profileC = NULL;

      int sib = Sibling(NJ,node);
      if (sib == -1) { /* at root, have 2 siblings */
	int sibs[2];
	RootSiblings(NJ, node, /*OUT*/sibs);
	profileB = NJ->profiles[sibs[0]];
	profileC = NJ->profiles[sibs[1]];
      } else {
	profileB = NJ->profiles[sib];
	profileC = GetUpProfile(/*IN/OUT*/upProfiles, NJ, NJ->parent[node], /*useML*/false);
      }
      profile_t *profiles[3] = {profileA,profileB,profileC};
      double d[3]; /*AB,AC,BC*/
      CorrectedPairDistances(profiles, 3, NJ->distance_matrix, NJ->nPos, /*OUT*/d);
      /* d(A,BC) = (dAB+dAC-dBC)/2 */
      NJ->branchlength[node] = (d[0]+d[1]-d[2])/2.0;
    } else {
      profile_t *profiles[4];
      int nodeABCD[4];
      SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);
      double d[6];
      CorrectedPairDistances(profiles, 4, NJ->distance_matrix, NJ->nPos, /*OUT*/d);
      NJ->branchlength[node] = (d[qAC]+d[qAD]+d[qBC]+d[qBD])/4.0 - (d[qAB]+d[qCD])/2.0;
      
      /* no longer needed */
      DeleteUpProfile(upProfiles, NJ, nodeABCD[0]);
      DeleteUpProfile(upProfiles, NJ, nodeABCD[1]);
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}