#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double *PSameVector(double length, rates_t *rates) {
  double *pSame = mymalloc(sizeof(double) * rates->nRateCategories);
  int iRate;
  for (iRate = 0; iRate < rates->nRateCategories; iRate++)
    pSame[iRate] = 0.25 + 0.75 * exp((-4.0/3.0) * fabs(length*rates->rates[iRate]));
  return(pSame);
}