#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

nni_stats_t *InitNNIStats(NJ_t *NJ) {
  nni_stats_t *stats = mymalloc(sizeof(nni_stats_t)*NJ->maxnode);
  const int LargeAge = 1000000;
  int i;
  for (i = 0; i < NJ->maxnode; i++) {
    stats[i].delta = 0;
    stats[i].support = 0;
    if (i == NJ->root || i < NJ->nSeq) {
      stats[i].age = LargeAge;
      stats[i].subtreeAge = LargeAge;
    } else {
      stats[i].age = 0;
      stats[i].subtreeAge = 0;
    }
  }
  return(stats);
}