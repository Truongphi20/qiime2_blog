#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double GammaLogLk(/*IN*/siteratelk_t *s, /*OPTIONAL OUT*/double *gamma_loglk_sites) {
  int iRate, iPos;
  double *dRate = mymalloc(sizeof(double) * s->nRateCats);
  for (iRate = 0; iRate < s->nRateCats; iRate++) {
    /* The probability density for each rate is approximated by the total
       density between the midpoints */
    double pMin = iRate == 0 ? 0.0 :
      PGamma(s->mult * (s->rates[iRate-1] + s->rates[iRate])/2.0, s->alpha);
    double pMax = iRate == s->nRateCats-1 ? 1.0 :
      PGamma(s->mult * (s->rates[iRate]+s->rates[iRate+1])/2.0, s->alpha);
    dRate[iRate] = pMax-pMin;
  }

  double loglk = 0.0;
  for (iPos = 0; iPos < s->nPos; iPos++) {
    /* Prevent underflow on large trees by comparing to maximum loglk */
    double maxloglk = -1e20;
    for (iRate = 0; iRate < s->nRateCats; iRate++) {
      double site_loglk = s->site_loglk[s->nPos*iRate + iPos];
      if (site_loglk > maxloglk)
	maxloglk = site_loglk;
    }
    double rellk = 0; /* likelihood scaled by exp(maxloglk) */
    for (iRate = 0; iRate < s->nRateCats; iRate++) {
      double lk = exp(s->site_loglk[s->nPos*iRate + iPos] - maxloglk);
      rellk += lk * dRate[iRate];
    }
    double loglk_site = maxloglk + log(rellk);
    loglk += loglk_site;
    if (gamma_loglk_sites != NULL)
      gamma_loglk_sites[iPos] = loglk_site;
  }
  dRate = myfree(dRate, sizeof(double)*s->nRateCats);
  return(loglk);
}
