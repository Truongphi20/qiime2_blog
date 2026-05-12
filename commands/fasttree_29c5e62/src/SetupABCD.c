#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void SetupABCD(NJ_t *NJ, int node,
	       /* the 4 profiles; the last one is an outprofile */
	       /*OPTIONAL OUT*/profile_t *profiles[4], 
	       /*OPTIONAL IN/OUT*/profile_t **upProfiles,
	       /*OUT*/int nodeABCD[4],
	       bool useML) {
  int parent = NJ->parent[node];
  assert(parent >= 0);
  assert(NJ->child[node].nChild == 2);
  nodeABCD[0] = NJ->child[node].child[0]; /*A*/
  nodeABCD[1] = NJ->child[node].child[1]; /*B*/

  profile_t *profile4 = NULL;
  if (parent == NJ->root) {
    int sibs[2];
    RootSiblings(NJ, node, /*OUT*/sibs);
    nodeABCD[2] = sibs[0];
    nodeABCD[3] = sibs[1];
    if (profiles == NULL)
      return;
    profile4 = NJ->profiles[sibs[1]];
  } else {
    nodeABCD[2] = Sibling(NJ,node);
    assert(nodeABCD[2] >= 0);
    nodeABCD[3] = parent;
    if (profiles == NULL)
      return;
    profile4 = GetUpProfile(upProfiles,NJ,parent,useML);
  }
  assert(upProfiles != NULL);
  int i;
  for (i = 0; i < 3; i++)
    profiles[i] = NJ->profiles[nodeABCD[i]];
  profiles[3] = profile4;
}