#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

nni_t MLQuartetNNI(profile_t *profiles[4],
		   /*OPTIONAL*/transition_matrix_t *transmat,
		   rates_t *rates,
		   int nPos, int nConstraints,
		   /*OUT*/double criteria[3], /* The three potential quartet log-likelihoods */
		   /*IN/OUT*/numeric_t len[5],
		   bool bFast)
{
  int i;
  double lenABvsCD[5] = {len[LEN_A], len[LEN_B], len[LEN_C], len[LEN_D], len[LEN_I]};
  double lenACvsBD[5] = {len[LEN_A], len[LEN_C], len[LEN_B], len[LEN_D], len[LEN_I]};   /* Swap B & C */
  double lenADvsBC[5] = {len[LEN_A], len[LEN_D], len[LEN_C], len[LEN_B], len[LEN_I]};   /* Swap B & D */
  bool bConsiderAC = true;
  bool bConsiderAD = true;
  int iRound;
  int nRounds = mlAccuracy < 2 ? 2 : mlAccuracy;
  double penalty[3];
  QuartetConstraintPenalties(profiles, nConstraints, /*OUT*/penalty);
  if (penalty[ABvsCD] > penalty[ACvsBD] || penalty[ABvsCD] > penalty[ADvsBC])
    bFast = false;
#ifdef OPENMP
      bFast = false;		/* turn off star topology test */
#endif

  for (iRound = 0; iRound < nRounds; iRound++) {
    bool bStarTest = false;
    {
#ifdef OPENMP
      #pragma omp parallel
      #pragma omp sections
#endif
      {
#ifdef OPENMP
        #pragma omp section
#endif
	{
	  criteria[ABvsCD] = MLQuartetOptimize(profiles[0], profiles[1], profiles[2], profiles[3],
					       nPos, transmat, rates,
					       /*IN/OUT*/lenABvsCD,
					       bFast ? &bStarTest : NULL,
					       /*site_likelihoods*/NULL)
	    - penalty[ABvsCD];	/* subtract penalty b/c we are trying to maximize log lk */
	}

#ifdef OPENMP
        #pragma omp section
#else
	if (bStarTest) {
	  nStarTests++;
	  criteria[ACvsBD] = -1e20;
	  criteria[ADvsBC] = -1e20;
	  len[LEN_I] = lenABvsCD[LEN_I];
	  return(ABvsCD);
	}
#endif
	{
	  if (bConsiderAC)
	    criteria[ACvsBD] = MLQuartetOptimize(profiles[0], profiles[2], profiles[1], profiles[3],
						 nPos, transmat, rates,
						 /*IN/OUT*/lenACvsBD, NULL, /*site_likelihoods*/NULL)
	      - penalty[ACvsBD];
	}
	
#ifdef OPENMP
        #pragma omp section
#endif
	{
	  if (bConsiderAD)
	    criteria[ADvsBC] = MLQuartetOptimize(profiles[0], profiles[3], profiles[2], profiles[1],
						 nPos, transmat, rates,
						 /*IN/OUT*/lenADvsBC, NULL, /*site_likelihoods*/NULL)
	      - penalty[ADvsBC];
	}
      }
    } /* end parallel sections */
    if (mlAccuracy < 2) {
      /* If clearly worse then ABvsCD, or have short internal branch length and worse, then
         give up */
      if (criteria[ACvsBD] < criteria[ABvsCD] - closeLogLkLimit
	  || (lenACvsBD[LEN_I] <= 2.0*MLMinBranchLength && criteria[ACvsBD] < criteria[ABvsCD]))
	bConsiderAC = false;
      if (criteria[ADvsBC] < criteria[ABvsCD] - closeLogLkLimit
	  || (lenADvsBC[LEN_I] <= 2.0*MLMinBranchLength && criteria[ADvsBC] < criteria[ABvsCD]))
	bConsiderAD = false;
      if (!bConsiderAC && !bConsiderAD)
	break;
      /* If clearly better than either alternative, then give up
         (Comparison is probably biased in favor of ABvsCD anyway) */
      if (criteria[ACvsBD] > criteria[ABvsCD] + closeLogLkLimit
	  && criteria[ACvsBD] > criteria[ADvsBC] + closeLogLkLimit)
	break;
      if (criteria[ADvsBC] > criteria[ABvsCD] + closeLogLkLimit
	  && criteria[ADvsBC] > criteria[ACvsBD] + closeLogLkLimit)
	break;
    }
  } /* end loop over rounds */

  if (verbose > 2) {
    fprintf(stderr, "Optimized quartet for %d rounds: ABvsCD %.5f ACvsBD %.5f ADvsBC %.5f\n",
	    iRound, criteria[ABvsCD], criteria[ACvsBD], criteria[ADvsBC]);
  }
  if (criteria[ACvsBD] > criteria[ABvsCD] && criteria[ACvsBD] > criteria[ADvsBC]) {
    for (i = 0; i < 5; i++) len[i] = lenACvsBD[i];
    return(ACvsBD);
  } else if (criteria[ADvsBC] > criteria[ABvsCD] && criteria[ADvsBC] > criteria[ACvsBD]) {
    for (i = 0; i < 5; i++) len[i] = lenADvsBC[i];
    return(ADvsBC);
  } else {
    for (i = 0; i < 5; i++) len[i] = lenABvsCD[i];
    return(ABvsCD);
  }
}