#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


double PairLogLk(profile_t *pA, profile_t *pB, double length, int nPos,
		 /*OPTIONAL*/transition_matrix_t *transmat,
		 rates_t *rates,
		 /*OPTIONAL IN/OUT*/double *site_likelihoods) {
  double lk = 1.0;
  double loglk = 0.0;		/* stores underflow of lk during the loop over positions */
  int i,j;
  assert(rates != NULL && rates->nRateCategories > 0);
  numeric_t *expeigenRates = NULL;
  if (transmat != NULL)
    expeigenRates = ExpEigenRates(length, transmat, rates);

  if (transmat == NULL) {	/* Jukes-Cantor */
    assert (nCodes == 4);
    double *pSame = PSameVector(length, rates);
    double *pDiff = PDiffVector(pSame, rates);
    
    int iFreqA = 0;
    int iFreqB = 0;
    for (i = 0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      double wA = pA->weights[i];
      double wB = pB->weights[i];
      int codeA = pA->codes[i];
      int codeB = pB->codes[i];
      numeric_t *fA = GET_FREQ(pA,i,/*IN/OUT*/iFreqA);
      numeric_t *fB = GET_FREQ(pB,i,/*IN/OUT*/iFreqB);
      double lkAB = 0;

      if (fA == NULL && fB == NULL) {
	if (codeA == NOCODE) {	/* A is all gaps */
	  /* gap to gap is sum(j) 0.25 * (0.25 * pSame + 0.75 * pDiff) = sum(i) 0.25*0.25 = 0.25
	     gap to any character gives the same result
	  */
	  lkAB = 0.25;
	} else if (codeB == NOCODE) { /* B is all gaps */
	  lkAB = 0.25;
	} else if (codeA == codeB) { /* A and B match */
	  lkAB = pSame[iRate] * wA*wB + 0.25 * (1-wA*wB);
	} else {		/* codeA != codeB */
	  lkAB = pDiff[iRate] * wA*wB + 0.25 * (1-wA*wB);
	}
      } else if (fA == NULL) {
	/* Compare codeA to profile of B */
	if (codeA == NOCODE)
	  lkAB = 0.25;
	else
	  lkAB = wA * (pDiff[iRate] + fB[codeA] * (pSame[iRate]-pDiff[iRate])) + (1.0-wA) * 0.25;
	/* because lkAB = wA * P(codeA->B) + (1-wA) * 0.25 
	   P(codeA -> B) = sum(j) P(B==j) * (j==codeA ? pSame : pDiff)
	   = sum(j) P(B==j) * pDiff + 
	   = pDiff + P(B==codeA) * (pSame-pDiff)
	*/
      } else if (fB == NULL) { /* Compare codeB to profile of A */
	if (codeB == NOCODE)
	  lkAB = 0.25;
	else
	  lkAB = wB * (pDiff[iRate] + fA[codeB] * (pSame[iRate]-pDiff[iRate])) + (1.0-wB) * 0.25;
      } else { /* both are full profiles */
	for (j = 0; j < 4; j++)
	  lkAB += fB[j] * (fA[j] * pSame[iRate] + (1-fA[j])* pDiff[iRate]); /* P(A|B) */
      }
      assert(lkAB > 0);
      lk *= lkAB;
      while (lk < LkUnderflow) {
	lk *= LkUnderflowInv;
	loglk -= LogLkUnderflow;
      }
      if (site_likelihoods != NULL)
	site_likelihoods[i] *= lkAB;
    }
    pSame = myfree(pSame, sizeof(double) * rates->nRateCategories);
    pDiff = myfree(pDiff, sizeof(double) * rates->nRateCategories);
  } else if (nCodes == 4) {	/* matrix model on nucleotides */
    int iFreqA = 0;
    int iFreqB = 0;
    numeric_t fAmix[4], fBmix[4];
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];

    for (i = 0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      numeric_t *expeigen = &expeigenRates[iRate*4];
      double wA = pA->weights[i];
      double wB = pB->weights[i];
      if (wA == 0 && wB == 0 && pA->codes[i] == NOCODE && pB->codes[i] == NOCODE) {
	/* Likelihood of A vs B is 1, so nothing changes
	   Do not need to advance iFreqA or iFreqB */
	continue;		
      }
      numeric_t *fA = GET_FREQ(pA,i,/*IN/OUT*/iFreqA);
      numeric_t *fB = GET_FREQ(pB,i,/*IN/OUT*/iFreqB);
      if (fA == NULL)
	fA = &transmat->codeFreq[pA->codes[i]][0];
      if (wA > 0.0 && wA < 1.0) {
	for (j  = 0; j < 4; j++)
	  fAmix[j] = wA*fA[j] + (1.0-wA)*fGap[j];
	fA = fAmix;
      }
      if (fB == NULL)
	fB = &transmat->codeFreq[pB->codes[i]][0];
      if (wB > 0.0 && wB < 1.0) {
	for (j  = 0; j < 4; j++)
	  fBmix[j] = wB*fB[j] + (1.0-wB)*fGap[j];
	fB = fBmix;
      }
      /* SSE3 instructions do not speed this step up:
	 numeric_t lkAB = vector_multiply3_sum(expeigen, fA, fB); */
		// dsp this is where check for <=0 was added in 2.1.1.LG
      double lkAB = 0;
      for (j = 0; j < 4; j++)
	lkAB += expeigen[j]*fA[j]*fB[j];
      assert(lkAB > 0);
      if (site_likelihoods != NULL)
	site_likelihoods[i] *= lkAB;
      lk *= lkAB;
      while (lk < LkUnderflow) {
	lk *= LkUnderflowInv;
	loglk -= LogLkUnderflow;
      }
      while (lk > LkUnderflowInv) {
	lk *= LkUnderflow;
	loglk += LogLkUnderflow;
      }
    }
  } else if (nCodes == 20) {	/* matrix model on amino acids */
    int iFreqA = 0;
    int iFreqB = 0;
    numeric_t fAmix[20], fBmix[20];
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];

    for (i = 0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      numeric_t *expeigen = &expeigenRates[iRate*20];
      double wA = pA->weights[i];
      double wB = pB->weights[i];
      if (wA == 0 && wB == 0 && pA->codes[i] == NOCODE && pB->codes[i] == NOCODE) {
	/* Likelihood of A vs B is 1, so nothing changes
	   Do not need to advance iFreqA or iFreqB */
	continue;		
      }
      numeric_t *fA = GET_FREQ(pA,i,/*IN/OUT*/iFreqA);
      numeric_t *fB = GET_FREQ(pB,i,/*IN/OUT*/iFreqB);
      if (fA == NULL)
	fA = &transmat->codeFreq[pA->codes[i]][0];
      if (wA > 0.0 && wA < 1.0) {
	for (j  = 0; j < 20; j++)
	  fAmix[j] = wA*fA[j] + (1.0-wA)*fGap[j];
	fA = fAmix;
      }
      if (fB == NULL)
	fB = &transmat->codeFreq[pB->codes[i]][0];
      if (wB > 0.0 && wB < 1.0) {
	for (j  = 0; j < 20; j++)
	  fBmix[j] = wB*fB[j] + (1.0-wB)*fGap[j];
	fB = fBmix;
      }
      numeric_t lkAB = vector_multiply3_sum(expeigen, fA, fB, 20);
      if (!(lkAB > 0)) {
	/* If this happens, it indicates a numerical problem that needs to be addressed elsewhere,
	   so report all the details */
	fprintf(stderr, "# FastTree.c::PairLogLk -- numerical problem!\n");
	fprintf(stderr, "# This block is intended for loading into R\n");

	fprintf(stderr, "lkAB = %.8g\n", lkAB);
	fprintf(stderr, "Branch_length= %.8g\nalignment_position=%d\nnCodes=%d\nrate_category=%d\nrate=%.8g\n",
		length, i, nCodes, iRate, rates->rates[iRate]);
	fprintf(stderr, "wA=%.8g\nwB=%.8g\n", wA, wB);
	fprintf(stderr, "codeA = %d\ncodeB = %d\n", pA->codes[i], pB->codes[i]);

	fprintf(stderr, "fA = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", fA[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "fB = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", fB[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "stat = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", transmat->stat[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "eigenval = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", transmat->eigenval[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "expeigen = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", expeigen[j]);
	fprintf(stderr,")\n");

	int k;
	fprintf(stderr, "codeFreq = c(");
	for (j = 0; j < nCodes; j++) for(k = 0; k < nCodes; k++) fprintf(stderr, "%s %.8g", j==0 && k==0?"":",",
									     transmat->codeFreq[j][k]);
	fprintf(stderr,")\n");

	fprintf(stderr, "eigeninv = c(");
	for (j = 0; j < nCodes; j++) for(k = 0; k < nCodes; k++) fprintf(stderr, "%s %.8g", j==0 && k==0?"":",",
									     transmat->eigeninv[j][k]);
	fprintf(stderr,")\n");

	fprintf(stderr, "# Transform into matrices and compute un-rotated vectors for profiles A and B\n");
	fprintf(stderr, "codeFreq = matrix(codeFreq,nrow=20);\n");
	fprintf(stderr, "eigeninv = matrix(eigeninv,nrow=20);\n");
	fputs("unrotA = stat * (eigeninv %*% fA)\n", stderr);
	fputs("unrotB = stat * (eigeninv %*% fB)\n", stderr);
	fprintf(stderr,"# End of R block\n");
      }
      assert(lkAB > 0);
      if (site_likelihoods != NULL)
	site_likelihoods[i] *= lkAB;
      lk *= lkAB;
      while (lk < LkUnderflow) {
	lk *= LkUnderflowInv;
	loglk -= LogLkUnderflow;
      }
      while (lk > LkUnderflowInv) {
	lk *= LkUnderflow;
	loglk += LogLkUnderflow;
      }
    }
  } else {
    assert(0);			/* illegal nCodes */
  }
  if (transmat != NULL)
    expeigenRates = myfree(expeigenRates, sizeof(numeric_t) * rates->nRateCategories * 20);
  loglk += log(lk);
  nLkCompute++;
  return(loglk);
}