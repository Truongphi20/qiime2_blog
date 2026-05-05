#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double PairNegLogLk(double x, void *data) {
  quartet_opt_t *qo = (quartet_opt_t *)data;
  assert(qo != NULL);
  assert(qo->pair1 != NULL && qo->pair2 != NULL);
  qo->nEval++;
  double loglk = PairLogLk(qo->pair1, qo->pair2, x, qo->nPos, qo->transmat, qo->rates, /*site_lk*/NULL);
  assert(loglk < 1e100);
  if (verbose > 5)
    fprintf(stderr, "PairLogLk(%.4f) =  %.4f\n", x, loglk);
  return(-loglk);
}