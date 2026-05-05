#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

numeric_t *ExpEigenRates(double length, transition_matrix_t *transmat, rates_t *rates) {
  numeric_t *expeigen = mymalloc(sizeof(numeric_t) * nCodes * rates->nRateCategories);
  int iRate, j;
  for (iRate = 0; iRate < rates->nRateCategories; iRate++) {
    for (j = 0; j < nCodes; j++) {
      double relLen = length * rates->rates[iRate];
      /* very short branch lengths lead to numerical problems so prevent them */
      if (relLen < MLMinRelBranchLength)
	relLen  = MLMinRelBranchLength;
      expeigen[iRate*nCodes + j] = exp(relLen * transmat->eigenval[j]);
    }
  }
  return(expeigen);
}