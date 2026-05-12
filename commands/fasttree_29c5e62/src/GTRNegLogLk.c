#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double GTRNegLogLk(double x, void *data) {
  
  gtr_opt_t *gtr = (gtr_opt_t*)data;
  assert(nCodes == 4);
  assert(gtr->NJ != NULL);
  assert(gtr->iRate >= 0 && gtr->iRate < 6);
  assert(x > 0);
  transition_matrix_t *old = gtr->NJ->transmat;
  double rates[6];
  int i;
  for (i = 0; i < 6; i++)
    rates[i] = gtr->rates[i];
  rates[gtr->iRate] = x;

  FILE *fpLog = gtr->fpLog;
  if (fpLog)
    fprintf(fpLog, "GTR_Opt\tfreq %.5f %.5f %.5f %.5f rates %.5f %.5f %.5f %.5f %.5f %.5f\n",
          gtr->freq[0], gtr->freq[1], gtr->freq[2], gtr->freq[3],
          rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]);

  gtr->NJ->transmat = CreateGTR(rates, gtr->freq);
  RecomputeMLProfiles(/*IN/OUT*/gtr->NJ);
  double loglk = TreeLogLk(gtr->NJ, /*site_loglk*/NULL);
  myfree(gtr->NJ->transmat, sizeof(transition_matrix_t));
  gtr->NJ->transmat = old;
  /* Do not recompute profiles -- assume the caller will do that */
  if (verbose > 2)
    fprintf(stderr, "GTR LogLk(%.5f %.5f %.5f %.5f %.5f %.5f) = %f\n",
	    rates[0], rates[1], rates[2], rates[3], rates[4], rates[5], loglk);
  if (fpLog)
    fprintf(fpLog, "GTR_Opt\tGTR LogLk(%.5f %.5f %.5f %.5f %.5f %.5f) = %f\n",
	    rates[0], rates[1], rates[2], rates[3], rates[4], rates[5], loglk);
  return(-loglk);
}