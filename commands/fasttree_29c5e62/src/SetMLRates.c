#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void SetMLRates(/*IN/OUT*/NJ_t *NJ, int nRateCategories) {
  assert(nRateCategories > 0);
  AllocRateCategories(/*IN/OUT*/&NJ->rates, 1, NJ->nPos); /* set to 1 category of rate 1 */
  if (nRateCategories == 1) {
    RecomputeMLProfiles(/*IN/OUT*/NJ);
    return;
  }
  numeric_t *rates = MLSiteRates(nRateCategories);
  double *site_loglk = MLSiteLikelihoodsByRate(/*IN*/NJ, /*IN*/rates, nRateCategories);

  /* Select best rate for each site, correcting for the prior
     For a prior, use a gamma distribution with shape parameter 3, scale 1/3, so
     Prior(rate) ~ rate**2 * exp(-3*rate)
     log Prior(rate) = C + 2 * log(rate) - 3 * rate
  */
  double sumRates = 0;
  int iPos;
  int iRate;
  for (iPos = 0; iPos < NJ->nPos; iPos++) {
    int iBest = -1;
    double dBest = -1e20;
    for (iRate = 0; iRate < nRateCategories; iRate++) {
      double site_loglk_with_prior = site_loglk[NJ->nPos*iRate + iPos]
	+ 2.0 * log(rates[iRate]) - 3.0 * rates[iRate];
      if (site_loglk_with_prior > dBest) {
	iBest = iRate;
	dBest = site_loglk_with_prior;
      }
    }
    if (verbose > 2)
      fprintf(stderr, "Selected rate category %d rate %.3f for position %d\n",
	      iBest, rates[iBest], iPos+1);
    NJ->rates.ratecat[iPos] = iBest;
    sumRates += rates[iBest];
  }
  site_loglk = myfree(site_loglk, sizeof(double)*NJ->nPos*nRateCategories);

  /* Force the rates to average to 1 */
  double avgRate = sumRates/NJ->nPos;
  for (iRate = 0; iRate < nRateCategories; iRate++)
    rates[iRate] /= avgRate;
  
  /* Save the rates */
  NJ->rates.rates = myfree(NJ->rates.rates, sizeof(numeric_t) * NJ->rates.nRateCategories);
  NJ->rates.rates = rates;
  NJ->rates.nRateCategories = nRateCategories;

  /* Update profiles based on rates */
  RecomputeMLProfiles(/*IN/OUT*/NJ);

  if (verbose) {
    fprintf(stderr, "Switched to using %d rate categories (CAT approximation)\n", nRateCategories);
    fprintf(stderr, "Rate categories were divided by %.3f so that average rate = 1.0\n", avgRate);
    fprintf(stderr, "CAT-based log-likelihoods may not be comparable across runs\n");
    if (!gammaLogLk)
      fprintf(stderr, "Use -gamma for approximate but comparable Gamma(20) log-likelihoods\n");
  }
}