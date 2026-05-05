#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double *PDiffVector(double *pSame, rates_t *rates) {
  double *pDiff = mymalloc(sizeof(double) * rates->nRateCategories);
  int iRate;
  for (iRate = 0; iRate < rates->nRateCategories; iRate++)
    pDiff[iRate] = (1.0 - pSame[iRate])/3.0;
  return(pDiff);
}