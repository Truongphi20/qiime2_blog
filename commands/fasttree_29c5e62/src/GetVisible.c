#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

bool GetVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		/*IN/OUT*/top_hits_t *tophits,
		int iNode, /*OUT*/besthit_t *visible) {
  if (iNode < 0 || NJ->parent[iNode] >= 0)
    return(false);
  hit_t *v = &tophits->visible[iNode];
  if (v->j < 0 || NJ->parent[v->j] >= 0)
    return(false);
  *visible = HitToBestHit(iNode, *v);
  SetCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/visible);  
  return(true);
}