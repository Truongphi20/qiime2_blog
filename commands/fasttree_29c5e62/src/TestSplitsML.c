#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void TestSplitsML(/*IN/OUT*/NJ_t *NJ, /*OUT*/SplitCount_t *splitcount, int nBootstrap) {
  const double tolerance = 1e-6;
  splitcount->nBadSplits = 0;
  splitcount->nConstraintViolations = 0;
  splitcount->nBadBoth = 0;
  splitcount->nSplits = 0;
  splitcount->dWorstDeltaUnconstrained = 0;
  splitcount->dWorstDeltaConstrained = 0;

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  int *col = nBootstrap > 0 ? ResampleColumns(NJ->nPos, nBootstrap) : NULL;
  double *site_likelihoods[3];
  int choice;
  for (choice = 0; choice < 3; choice++)
    site_likelihoods[choice] = mymalloc(sizeof(double)*NJ->nPos);

  int iNodesDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */
    
    if(iNodesDone > 0 && (iNodesDone % 100) == 0)
      ProgressReport("ML split tests for %6d of %6d internal splits", iNodesDone, NJ->nSeq-3, 0, 0);
    iNodesDone++;

    profile_t *profiles[4];
    int nodeABCD[4];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/true);
    double loglk[3];
    double len[5];
    int i;
    for (i = 0; i < 4; i++)
      len[i] = NJ->branchlength[nodeABCD[i]];
    len[4] = NJ->branchlength[node];
    double lenABvsCD[5] = {len[LEN_A], len[LEN_B], len[LEN_C], len[LEN_D], len[LEN_I]};
    double lenACvsBD[5] = {len[LEN_A], len[LEN_C], len[LEN_B], len[LEN_D], len[LEN_I]};   /* Swap B & C */
    double lenADvsBC[5] = {len[LEN_A], len[LEN_D], len[LEN_C], len[LEN_B], len[LEN_I]};   /* Swap B & D */

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
	  /* Lengths are already optimized for ABvsCD */
	  loglk[ABvsCD] = MLQuartetLogLk(profiles[0], profiles[1], profiles[2], profiles[3],
					 NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenABvsCD,
					 /*OUT*/site_likelihoods[ABvsCD]);
	}

#ifdef OPENMP
      #pragma omp section
#endif
	{
	  loglk[ACvsBD] = MLQuartetOptimize(profiles[0], profiles[2], profiles[1], profiles[3],
					    NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenACvsBD, /*pStarTest*/NULL,
					    /*OUT*/site_likelihoods[ACvsBD]);
	}

#ifdef OPENMP
      #pragma omp section
#endif
	{
	  loglk[ADvsBC] = MLQuartetOptimize(profiles[0], profiles[3], profiles[2], profiles[1],
					    NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenADvsBC, /*pStarTest*/NULL,
					    /*OUT*/site_likelihoods[ADvsBC]);
	}
      }
    }

    /* do a second pass on the better alternative if it is close */
    if (loglk[ACvsBD] > loglk[ADvsBC]) {
      if (mlAccuracy > 1 || loglk[ACvsBD] > loglk[ABvsCD] - closeLogLkLimit) {
	loglk[ACvsBD] = MLQuartetOptimize(profiles[0], profiles[2], profiles[1], profiles[3],
					  NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenACvsBD, /*pStarTest*/NULL,
					  /*OUT*/site_likelihoods[ACvsBD]);
      }
    } else {
      if (mlAccuracy > 1 || loglk[ADvsBC] > loglk[ABvsCD] - closeLogLkLimit) {
	loglk[ADvsBC] = MLQuartetOptimize(profiles[0], profiles[3], profiles[2], profiles[1],
					  NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenADvsBC, /*pStarTest*/NULL,
					  /*OUT*/site_likelihoods[ADvsBC]);
      }
    }

    if (loglk[ABvsCD] >= loglk[ACvsBD] && loglk[ABvsCD] >= loglk[ADvsBC])
      choice = ABvsCD;
    else if (loglk[ACvsBD] >= loglk[ABvsCD] && loglk[ACvsBD] >= loglk[ADvsBC])
      choice = ACvsBD;
    else
      choice = ADvsBC;
    bool badSplit = loglk[choice] > loglk[ABvsCD] + treeLogLkDelta; /* ignore small changes in likelihood */

    /* constraint penalties, indexed by nni_t (lower is better) */
    double p[3];
    QuartetConstraintPenalties(profiles, NJ->nConstraints, /*OUT*/p);
    bool bBadConstr = p[ABvsCD] > p[ACvsBD] + tolerance || p[ABvsCD] > p[ADvsBC] + tolerance;
    bool violateConstraint = false;
    int iC;
    for (iC=0; iC < NJ->nConstraints; iC++) {
      if (SplitViolatesConstraint(profiles, iC)) {
	violateConstraint = true;
	break;
      }
    }
    splitcount->nSplits++;
    if (violateConstraint)
      splitcount->nConstraintViolations++;
    if (badSplit)
      splitcount->nBadSplits++;
    if (badSplit && bBadConstr)
      splitcount->nBadBoth++;
    if (badSplit) {
      double delta = loglk[choice] - loglk[ABvsCD];
      /* If ABvsCD is favored over the more likely NNI by constraints,
	 then this is probably a bad split because of the constraint */
      if (p[choice] > p[ABvsCD] + tolerance)
	splitcount->dWorstDeltaConstrained = MAX(delta, splitcount->dWorstDeltaConstrained);
      else
	splitcount->dWorstDeltaUnconstrained = MAX(delta, splitcount->dWorstDeltaUnconstrained);
    }
    if (nBootstrap>0)
      NJ->support[node] = badSplit ? 0.0 : SHSupport(NJ->nPos, nBootstrap, col, loglk, site_likelihoods);

    /* No longer needed */
    DeleteUpProfile(upProfiles, NJ, nodeABCD[0]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[1]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[2]);
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  if (nBootstrap>0)
    col = myfree(col, sizeof(int)*((size_t)NJ->nPos)*nBootstrap);
  for (choice = 0; choice < 3; choice++)
    site_likelihoods[choice] = myfree(site_likelihoods[choice], sizeof(double)*NJ->nPos);
}