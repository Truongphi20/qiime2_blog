#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double OptAlpha(double alpha, void *data) {
  siteratelk_t *s = (siteratelk_t *)data;
  s->alpha = alpha;
  return(-GammaLogLk(s, NULL));
}

double OptMult(double mult, void *data) {
  siteratelk_t *s = (siteratelk_t *)data;
  s->mult = mult;
  return(-GammaLogLk(s, NULL));
}

double RescaleGammaLogLk(int nPos, int nRateCats, /*IN*/numeric_t *rates, /*IN*/double *site_loglk,
			 /*OPTIONAL*/FILE *fpLog) {
  siteratelk_t s = { /*mult*/1.0, /*alpha*/1.0, nPos, nRateCats, rates, site_loglk };
  double fx, f2x;
  int i;
  fx = -GammaLogLk(&s, NULL);
  if (verbose>2)
    fprintf(stderr, "Optimizing alpha, starting at loglk %.3f\n", -fx);
  for (i = 0; i < 10; i++) {
    ProgressReport("Optimizing alpha round %d", i+1, 0, 0, 0);
    double start = fx;
    s.alpha = onedimenmin(0.01, s.alpha, 10.0, OptAlpha, &s, 0.001, 0.001, &fx, &f2x);
    if (verbose>2)
      fprintf(stderr, "Optimize alpha round %d to %.3f lk %.3f\n", i+1, s.alpha, -fx);
    s.mult = onedimenmin(0.01, s.mult, 10.0, OptMult, &s, 0.001, 0.001, &fx, &f2x);
    if (verbose>2)
      fprintf(stderr, "Optimize mult round %d to %.3f lk %.3f\n", i+1, s.mult, -fx);
    if (fx > start - 0.001) {
      if (verbose>2)
	fprintf(stderr, "Optimizing alpha & mult converged\n");
      break;
    }
  }

  double *gamma_loglk_sites = mymalloc(sizeof(double) * nPos);
  double gammaLogLk = GammaLogLk(&s, /*OUT*/gamma_loglk_sites);
  if (verbose > 0)
    fprintf(stderr, "Gamma(%d) LogLk = %.3f alpha = %.3f rescaling lengths by %.3f\n",
	    nRateCats, gammaLogLk, s.alpha, 1/s.mult);
  if (fpLog) {
    int iPos;
    int iRate;
    fprintf(fpLog, "Gamma%dLogLk\t%.3f\tApproximate\tAlpha\t%.3f\tRescale\t%.3f\n",
	    nRateCats, gammaLogLk, s.alpha, 1/s.mult);
    fprintf(fpLog, "Gamma%d\tSite\tLogLk", nRateCats);
    for (iRate = 0; iRate < nRateCats; iRate++)
      fprintf(fpLog, "\tr=%.3f", rates[iRate]/s.mult);
    fprintf(fpLog,"\n");
    for (iPos = 0; iPos < nPos; iPos++) {
      fprintf(fpLog, "Gamma%d\t%d\t%.3f", nRateCats, iPos, gamma_loglk_sites[iPos]);
      for (iRate = 0; iRate < nRateCats; iRate++)
	fprintf(fpLog, "\t%.3f", site_loglk[nPos*iRate + iPos]);
      fprintf(fpLog,"\n");
    }
  }
  gamma_loglk_sites = myfree(gamma_loglk_sites, sizeof(double) * nPos);
  return(1.0/s.mult);
}