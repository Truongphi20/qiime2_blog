#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double TreeLogLk(/*IN*/NJ_t *NJ, /*OPTIONAL OUT*/double *site_loglk) {
  int i;
  if (NJ->nSeq < 2)
    return(0.0);
  double loglk = 0.0;
  double *site_likelihood = NULL;
  if (site_loglk != NULL) {
    site_likelihood = mymalloc(sizeof(double)*NJ->nPos);
    for (i = 0; i < NJ->nPos; i++) {
      site_likelihood[i] = 1.0;
      site_loglk[i] = 0.0;
    }
  }
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    int nChild = NJ->child[node].nChild;
    if (nChild == 0)
      continue;
    assert(nChild >= 2);
    int *children = NJ->child[node].child;
    double loglkchild = PairLogLk(NJ->profiles[children[0]], NJ->profiles[children[1]],
				  NJ->branchlength[children[0]]+NJ->branchlength[children[1]],
				  NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/site_likelihood);
    loglk += loglkchild;
    if (site_likelihood != NULL) {
      /* prevent underflows */
      for (i = 0; i < NJ->nPos; i++) {
	while(site_likelihood[i] < LkUnderflow) {
	  site_likelihood[i] *= LkUnderflowInv;
	  site_loglk[i] -= LogLkUnderflow;
	}
      }
    }
    if (verbose > 2)
      fprintf(stderr, "At %d: LogLk(%d:%.4f,%d:%.4f) = %.3f\n",
	      node,
	      children[0], NJ->branchlength[children[0]],
	      children[1], NJ->branchlength[children[1]],
	      loglkchild);
    if (NJ->child[node].nChild == 3) {
      assert(node == NJ->root);
      /* Infer the common parent of the 1st two to define the third... */
      profile_t *pAB = PosteriorProfile(NJ->profiles[children[0]],
					NJ->profiles[children[1]],
					NJ->branchlength[children[0]],
					NJ->branchlength[children[1]],
					NJ->transmat, &NJ->rates,
					NJ->nPos, /*nConstraints*/0);
      double loglkup = PairLogLk(pAB, NJ->profiles[children[2]],
				 NJ->branchlength[children[2]],
				 NJ->nPos, NJ->transmat, &NJ->rates,
				 /*IN/OUT*/site_likelihood);
      loglk += loglkup;
      if (verbose > 2)
	fprintf(stderr, "At root %d: LogLk((%d/%d),%d:%.3f) = %.3f\n",
		node, children[0], children[1], children[2],
		NJ->branchlength[children[2]],
		loglkup);
      pAB = FreeProfile(pAB, NJ->nPos, NJ->nConstraints);
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  if (site_likelihood != NULL) {
    for (i = 0; i < NJ->nPos; i++) {
      site_loglk[i] += log(site_likelihood[i]);
    }
    site_likelihood = myfree(site_likelihood, sizeof(double)*NJ->nPos);
  }

  /* For Jukes-Cantor, with a tree of size 4, if the children of the root are
     (A,B), C, and D, then
     P(ABCD) = P(A) P(B|A) P(C|AB) P(D|ABC)
     
     Above we compute P(B|A) P(C|AB) P(D|ABC) -- note P(B|A) is at the child of root
     and P(C|AB) P(D|ABC) is at root.

     Similarly if the children of the root are C, D, and (A,B), then
     P(ABCD) = P(C|D) P(A|B) P(AB|CD) P(D), and above we compute that except for P(D)

     So we need to multiply by P(A) = 0.25, so we pay log(4) at each position
     (if ungapped). Each gapped position in any sequence reduces the payment by log(4)

     For JTT or GTR, we are computing P(A & B) and the posterior profiles are scaled to take
     the prior into account, so we do not need any correction.
     codeFreq[NOCODE] is scaled x higher so that P(-) = 1 not P(-)=1/nCodes, so gaps
     do not need to be corrected either.
   */

  if (nCodes == 4 && NJ->transmat == NULL) {
    int nGaps = 0;
    double logNCodes = log((double)nCodes);
    for (i = 0; i < NJ->nPos; i++) {
      int nGapsThisPos = 0;
      for (node = 0; node < NJ->nSeq; node++) {
	unsigned char *codes = NJ->profiles[node]->codes;
	if (codes[i] == NOCODE)
	  nGapsThisPos++;
      }
      nGaps += nGapsThisPos;
      if (site_loglk != NULL) {
	site_loglk[i] += nGapsThisPos * logNCodes;
	if (nCodes == 4 && NJ->transmat == NULL)
	  site_loglk[i] -= logNCodes;
      }
    }
    loglk -= NJ->nPos * logNCodes;
    loglk += nGaps * logNCodes;	/* do not pay for gaps -- only Jukes-Cantor */
  }
  return(loglk);
}