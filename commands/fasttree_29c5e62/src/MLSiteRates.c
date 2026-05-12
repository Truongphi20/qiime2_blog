#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Caller must free the resulting vector of n rates */
numeric_t *MLSiteRates(int nRateCategories) {
  /* Even spacing from 1/nRate to nRate */
  double logNCat = log((double)nRateCategories);
  double logMinRate = -logNCat;
  double logMaxRate = logNCat;
  double logd = (logMaxRate-logMinRate)/(double)(nRateCategories-1);

  numeric_t *rates = mymalloc(sizeof(numeric_t)*nRateCategories);
  int i;
  for (i = 0; i < nRateCategories; i++)
    rates[i] = exp(logMinRate + logd*(double)i);
  return(rates);
}