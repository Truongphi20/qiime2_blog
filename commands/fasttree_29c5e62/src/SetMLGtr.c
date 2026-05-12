#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void SetMLGtr(/*IN/OUT*/NJ_t *NJ, /*OPTIONAL IN*/double *freq_in, /*OPTIONAL WRITE*/FILE *fpLog) {
  int i;
  assert(nCodes==4);
  gtr_opt_t gtr;
  gtr.NJ = NJ;
  gtr.fpLog = fpLog;
  if (freq_in != NULL) {
    for (i=0; i<4; i++)
      gtr.freq[i]=freq_in[i];
  } else {
    /* n[] and sum were int in FastTree 2.1.9 and earlier -- this
       caused gtr analyses to fail on analyses with >2e9 positions */
    long n[4] = {1,1,1,1};	/* pseudocounts */
    for (i=0; i<NJ->nSeq; i++) {
      unsigned char *codes = NJ->profiles[i]->codes;
      int iPos;
      for (iPos=0; iPos<NJ->nPos; iPos++)
	if (codes[iPos] < 4)
	  n[codes[iPos]]++;
    }
    long sum = n[0]+n[1]+n[2]+n[3];
    for (i=0; i<4; i++)
      gtr.freq[i] = n[i]/(double)sum;
  }
  for (i=0; i<6; i++)
    gtr.rates[i] = 1.0;
  int nRounds = mlAccuracy < 2 ? 2 : mlAccuracy;
  for (i = 0; i < nRounds; i++) {
    for (gtr.iRate = 0; gtr.iRate < 6; gtr.iRate++) {
      ProgressReport("Optimizing GTR model, step %d of %d", i*6+gtr.iRate+1, 12, 0, 0);
      double negloglk, f2x;
      gtr.rates[gtr.iRate] = onedimenmin(/*xmin*/0.05,
					 /*xguess*/gtr.rates[gtr.iRate],
					 /*xmax*/20.0,
					 GTRNegLogLk,
					 /*data*/&gtr,
					 /*ftol*/0.001,
					 /*atol*/0.0001,
					 /*OUT*/&negloglk,
					 /*OUT*/&f2x);
    }
  }
  /* normalize gtr so last rate is 1 -- specifying that rate separately is useful for optimization only */
  for (i = 0; i < 5; i++)
    gtr.rates[i] /= gtr.rates[5];
  gtr.rates[5] = 1.0;
  if (verbose) {
    fprintf(stderr, "GTR Frequencies: %.4f %.4f %.4f %.4f\n", gtr.freq[0], gtr.freq[1], gtr.freq[2], gtr.freq[3]);
    fprintf(stderr, "GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f\n",
	    gtr.rates[0],gtr.rates[1],gtr.rates[2],gtr.rates[3],gtr.rates[4],gtr.rates[5]);
  }
  if (fpLog != NULL) {
    fprintf(fpLog, "GTRFreq\t%.4f\t%.4f\t%.4f\t%.4f\n", gtr.freq[0], gtr.freq[1], gtr.freq[2], gtr.freq[3]);
    fprintf(fpLog, "GTRRates\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	    gtr.rates[0],gtr.rates[1],gtr.rates[2],gtr.rates[3],gtr.rates[4],gtr.rates[5]);
  }
  myfree(NJ->transmat, sizeof(transition_matrix_t));
  NJ->transmat = CreateGTR(gtr.rates, gtr.freq);
  RecomputeMLProfiles(/*IN/OUT*/NJ);
  OptimizeAllBranchLengths(/*IN/OUT*/NJ);
}