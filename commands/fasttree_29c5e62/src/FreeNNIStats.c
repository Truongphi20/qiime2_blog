#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

nni_stats_t *FreeNNIStats(nni_stats_t *stats, NJ_t *NJ) {
  return(myfree(stats, sizeof(nni_stats_t)*NJ->maxnode));
}