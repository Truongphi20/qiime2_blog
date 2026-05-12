#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double MLQuartetLogLk(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
		      int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		      /*IN*/double branch_lengths[5],
		      /*OPTIONAL OUT*/double *site_likelihoods) {
  profile_t *pAB = PosteriorProfile(pA, pB,
				    branch_lengths[0], branch_lengths[1],
				    transmat,
				    rates,
				    nPos, /*nConstraints*/0);
  profile_t *pCD = PosteriorProfile(pC, pD,
				    branch_lengths[2], branch_lengths[3],
				    transmat,
				    rates,
				    nPos, /*nConstraints*/0);
  if (site_likelihoods != NULL) {
    int i;
    for (i = 0; i < nPos; i++)
      site_likelihoods[i] = 1.0;
  }
  /* Roughly, P(A,B,C,D) = P(A) P(B|A) P(D|C) P(AB | CD) */
  double loglk = PairLogLk(pA, pB, branch_lengths[0]+branch_lengths[1],
			   nPos, transmat, rates, /*OPTIONAL IN/OUT*/site_likelihoods)
    + PairLogLk(pC, pD, branch_lengths[2]+branch_lengths[3],
		nPos, transmat, rates, /*OPTIONAL IN/OUT*/site_likelihoods)
    + PairLogLk(pAB, pCD, branch_lengths[4],
		nPos, transmat, rates, /*OPTIONAL IN/OUT*/site_likelihoods);
  pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);
  pCD = FreeProfile(pCD, nPos, /*nConstraints*/0);
  return(loglk);
}