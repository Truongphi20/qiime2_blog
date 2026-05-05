#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


double MLQuartetOptimize(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
			 int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
			 /*IN/OUT*/double branch_lengths[5],
			 /*OPTIONAL OUT*/bool *pStarTest,
			 /*OPTIONAL OUT*/double *site_likelihoods) {
  int j;
  double start_length[5];
  for (j = 0; j < 5; j++) {
    start_length[j] = branch_lengths[j];
    if (branch_lengths[j] < MLMinBranchLength)
      branch_lengths[j] = MLMinBranchLength;
  }
  quartet_opt_t qopt = { nPos, transmat, rates, /*nEval*/0,
			 /*pair1*/NULL, /*pair2*/NULL };
  double f2x, negloglk;

  if (pStarTest != NULL)
    *pStarTest = false;

  /* First optimize internal branch, then branch to A, B, C, D, in turn
     May use star test to quit after internal branch
   */
  profile_t *pAB = PosteriorProfile(pA, pB,
				    branch_lengths[LEN_A], branch_lengths[LEN_B],
				    transmat, rates, nPos, /*nConstraints*/0);
  profile_t *pCD = PosteriorProfile(pC, pD,
				    branch_lengths[LEN_C], branch_lengths[LEN_D],
				    transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pAB;
  qopt.pair2 = pCD;
  branch_lengths[LEN_I] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_I],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);

  if (pStarTest != NULL) {
    assert(site_likelihoods == NULL);
    double loglkStar = -PairNegLogLk(MLMinBranchLength, &qopt);
    if (loglkStar < -negloglk - closeLogLkLimit) {
      *pStarTest = true;
      double off = PairLogLk(pA, pB,
			     branch_lengths[LEN_A] + branch_lengths[LEN_B],
			     qopt.nPos, qopt.transmat, qopt.rates, /*site_lk*/NULL)
	+ PairLogLk(pC, pD,
		    branch_lengths[LEN_C] + branch_lengths[LEN_D],
		    qopt.nPos, qopt.transmat, qopt.rates, /*site_lk*/NULL);
      pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);
      pCD = FreeProfile(pCD, nPos, /*nConstraints*/0);
      return (-negloglk + off);
    }
  }
  pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);
  profile_t *pBCD = PosteriorProfile(pB, pCD,
				     branch_lengths[LEN_B], branch_lengths[LEN_I],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pA;
  qopt.pair2 = pBCD;
  branch_lengths[LEN_A] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_A],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);
  pBCD = FreeProfile(pBCD, nPos, /*nConstraints*/0);
  profile_t *pACD = PosteriorProfile(pA, pCD,
				     branch_lengths[LEN_A], branch_lengths[LEN_I],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pB;
  qopt.pair2 = pACD;
  branch_lengths[LEN_B] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_B],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);
  pACD = FreeProfile(pACD, nPos, /*nConstraints*/0);
  pCD = FreeProfile(pCD, nPos, /*nConstraints*/0);
  pAB = PosteriorProfile(pA, pB,
			 branch_lengths[LEN_A], branch_lengths[LEN_B],
			 transmat, rates, nPos, /*nConstraints*/0);
  profile_t *pABD = PosteriorProfile(pAB, pD,
				     branch_lengths[LEN_I], branch_lengths[LEN_D],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pC;
  qopt.pair2 = pABD;
  branch_lengths[LEN_C] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_C],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);
  pABD = FreeProfile(pABD, nPos, /*nConstraints*/0);
  profile_t *pABC = PosteriorProfile(pAB, pC,
				     branch_lengths[LEN_I], branch_lengths[LEN_C],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pD;
  qopt.pair2 = pABC;
  branch_lengths[LEN_D] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_D],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);

  /* Compute the total quartet likelihood
     PairLogLk(ABC,D) + PairLogLk(AB,C) + PairLogLk(A,B)
   */
  double loglkABCvsD = -negloglk;
  if (site_likelihoods) {
    for (j = 0; j < nPos; j++)
      site_likelihoods[j] = 1.0;
    PairLogLk(pABC, pD, branch_lengths[LEN_D],
	      qopt.nPos, qopt.transmat, qopt.rates, /*IN/OUT*/site_likelihoods);
  }
  double quartetloglk = loglkABCvsD
    + PairLogLk(pAB, pC, branch_lengths[LEN_I] + branch_lengths[LEN_C],
		qopt.nPos, qopt.transmat, qopt.rates,
		/*IN/OUT*/site_likelihoods)
    + PairLogLk(pA, pB, branch_lengths[LEN_A] + branch_lengths[LEN_B],
		qopt.nPos, qopt.transmat, qopt.rates,
		/*IN/OUT*/site_likelihoods);

  pABC = FreeProfile(pABC, nPos, /*nConstraints*/0);
  pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);

  if (verbose > 3) {
    double loglkStart = MLQuartetLogLk(pA, pB, pC, pD, nPos, transmat, rates, start_length, /*site_lk*/NULL);
    fprintf(stderr, "Optimize loglk from %.5f to %.5f eval %d lengths from\n"
	    "   %.5f %.5f %.5f %.5f %.5f to\n"
	    "   %.5f %.5f %.5f %.5f %.5f\n",
	    loglkStart, quartetloglk, qopt.nEval,
	    start_length[0], start_length[1], start_length[2], start_length[3], start_length[4],
	    branch_lengths[0], branch_lengths[1], branch_lengths[2], branch_lengths[3], branch_lengths[4]);
  }
  return(quartetloglk);
}