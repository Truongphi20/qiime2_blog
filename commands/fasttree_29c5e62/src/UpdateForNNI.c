#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Update the profile of node and its ancestor, and delete nearby out-profiles */
void UpdateForNNI(/*IN/OUT*/NJ_t *NJ, int node, /*IN/OUT*/profile_t **upProfiles,
		  bool useML) {
  int i;
  if (slow) {
    /* exhaustive update */
    for (i = 0; i < NJ->maxnodes; i++)
      DeleteUpProfile(upProfiles, NJ, i);

    /* update profiles back to root */
    int ancestor;
    for (ancestor = node; ancestor >= 0; ancestor = NJ->parent[ancestor])
      RecomputeProfile(/*IN/OUT*/NJ, upProfiles, ancestor, useML);

    /* remove any up-profiles made while doing that*/
    for (i = 0; i < NJ->maxnodes; i++)
      DeleteUpProfile(upProfiles, NJ, i);
  } else {
    /* if fast, only update around self
       note that upProfile(parent) is still OK after an NNI, but
       up-profiles of uncles may not be
    */
    DeleteUpProfile(upProfiles, NJ, node);
    for (i = 0; i < NJ->child[node].nChild; i++)
      DeleteUpProfile(upProfiles, NJ, NJ->child[node].child[i]);
    assert(node != NJ->root);
    int parent = NJ->parent[node];
    int neighbors[2] = { parent, Sibling(NJ, node) };
    if (parent == NJ->root)
      RootSiblings(NJ, node, /*OUT*/neighbors);
    DeleteUpProfile(upProfiles, NJ, neighbors[0]);
    DeleteUpProfile(upProfiles, NJ, neighbors[1]);
    int uncle = Sibling(NJ, parent);
    if (uncle >= 0)
      DeleteUpProfile(upProfiles, NJ, uncle);
    RecomputeProfile(/*IN/OUT*/NJ, upProfiles, node, useML);
    RecomputeProfile(/*IN/OUT*/NJ, upProfiles, parent, useML);
  }
}