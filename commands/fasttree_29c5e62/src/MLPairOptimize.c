#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Input site_loglk must be for each rate */
double MLPairOptimize(profile_t *pA, profile_t *pB,
		      int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		      /*IN/OUT*/double *branch_length) {
  quartet_opt_t qopt = { nPos, transmat, rates,
			 /*nEval*/0, /*pair1*/pA, /*pair2*/pB };
  double f2x,negloglk;
  *branch_length = onedimenmin(/*xmin*/MLMinBranchLength,
			       /*xguess*/*branch_length,
			       /*xmax*/6.0,
			       PairNegLogLk,
			       /*data*/&qopt,
			       /*ftol*/MLFTolBranchLength,
			       /*atol*/MLMinBranchLengthTolerance,
			       /*OUT*/&negloglk,
			       /*OUT*/&f2x);
  return(-negloglk);		/* the log likelihood */
}