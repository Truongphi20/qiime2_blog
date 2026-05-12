#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double *MLSiteLikelihoodsByRate(/*IN*/NJ_t *NJ, /*IN*/numeric_t *rates, int nRateCategories) {
  double *site_loglk = mymalloc(sizeof(double)*NJ->nPos*nRateCategories);

  /* save the original rates */
  assert(NJ->rates.nRateCategories > 0);
  numeric_t *oldRates = NJ->rates.rates;
  NJ->rates.rates = mymalloc(sizeof(numeric_t) * NJ->rates.nRateCategories);

  /* Compute site likelihood for each rate */
  int iPos;
  int iRate;
  for (iRate = 0; iRate  < nRateCategories; iRate++) {
    int i;
    for (i = 0; i < NJ->rates.nRateCategories; i++)
      NJ->rates.rates[i] = rates[iRate];
    RecomputeMLProfiles(/*IN/OUT*/NJ);
    double loglk = TreeLogLk(NJ, /*OUT*/&site_loglk[NJ->nPos*iRate]);
    ProgressReport("Site likelihoods with rate category %d of %d", iRate+1, nRateCategories, 0, 0);
    if(verbose > 2) {
      fprintf(stderr, "Rate %.3f Loglk %.3f SiteLogLk", rates[iRate], loglk);
      for (iPos = 0; iPos < NJ->nPos; iPos++)
	fprintf(stderr,"\t%.3f", site_loglk[NJ->nPos*iRate + iPos]);
      fprintf(stderr,"\n");
    }
  }

  /* restore original rates and profiles */
  myfree(NJ->rates.rates, sizeof(numeric_t) * NJ->rates.nRateCategories);
  NJ->rates.rates = oldRates;
  RecomputeMLProfiles(/*IN/OUT*/NJ);

  return(site_loglk);
}