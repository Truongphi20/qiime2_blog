#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


profile_t *PosteriorProfile(profile_t *p1, profile_t *p2,
			    double len1, double len2,
			    /*OPTIONAL*/transition_matrix_t *transmat,
			    rates_t *rates,
			    int nPos, int nConstraints) {
  if (len1 < MLMinBranchLength)
    len1 = MLMinBranchLength;
  if (len2 < MLMinBranchLength)
    len2 = MLMinBranchLength;

  int i,j,k;
  profile_t *out = NewProfile(nPos, nConstraints);
  for (i = 0; i < nPos; i++) {
    out->codes[i] = NOCODE;
    out->weights[i] = 1.0;
  }
  out->nVectors = nPos;
  out->vectors = (numeric_t*)mymalloc(sizeof(numeric_t)*nCodes*out->nVectors);
  for (i = 0; i < nCodes * out->nVectors; i++) out->vectors[i] = 0;
  int iFreqOut = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  numeric_t *expeigenRates1 = NULL, *expeigenRates2 = NULL;

  if (transmat != NULL) {
    expeigenRates1 = ExpEigenRates(len1, transmat, rates);
    expeigenRates2 = ExpEigenRates(len2, transmat, rates);
  }

  if (transmat == NULL) {	/* Jukes-Cantor */
    assert(nCodes == 4);

    double *PSame1 = PSameVector(len1, rates);
    double *PDiff1 = PDiffVector(PSame1, rates);
    double *PSame2 = PSameVector(len2, rates);
    double *PDiff2 = PDiffVector(PSame2, rates);

    numeric_t mix1[4], mix2[4];

    for (i=0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      double w1 = p1->weights[i];
      double w2 = p2->weights[i];
      int code1 = p1->codes[i];
      int code2 = p2->codes[i];
      numeric_t *f1 = GET_FREQ(p1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(p2,i,/*IN/OUT*/iFreq2);

      /* First try to store a simple profile */
      if (f1 == NULL && f2 == NULL) {
	if (code1 == NOCODE && code2 == NOCODE) {
	  out->codes[i] = NOCODE;
	  out->weights[i] = 0.0;
	  continue;
	} else if (code1 == NOCODE) {
	  /* Posterior(parent | character & gap, len1, len2) = Posterior(parent | character, len1)
	     = PSame() for matching characters and 1-PSame() for the rest
	     = (pSame - pDiff) * character + (1-(pSame-pDiff)) * gap
	  */
	  out->codes[i] = code2;
	  out->weights[i] = w2 * (PSame2[iRate] - PDiff2[iRate]);
	  continue;
	} else if (code2 == NOCODE) {
	  out->codes[i] = code1;
	  out->weights[i] = w1 * (PSame1[iRate] - PDiff1[iRate]);
	  continue;
	} else if (code1 == code2) {
	  out->codes[i] = code1;
	  double f12code = (w1*PSame1[iRate] + (1-w1)*0.25) * (w2*PSame2[iRate] + (1-w2)*0.25);
	  double f12other = (w1*PDiff1[iRate] + (1-w1)*0.25) * (w2*PDiff2[iRate] + (1-w2)*0.25);
	  /* posterior probability of code1/code2 after scaling */
	  double pcode = f12code/(f12code+3*f12other);
	  /* Now f = w * (code ? 1 : 0) + (1-w) * 0.25, so to get pcode we need
	     fcode = 1/4 + w1*3/4 or w = (f-1/4)*4/3
	   */
	  out->weights[i] = (pcode - 0.25) * 4.0/3.0;
	  /* This can be zero because of numerical problems, I think */
	  if (out->weights[i] < 1e-6) {
	    if (verbose > 1)
	      fprintf(stderr, "Replaced weight %f with %f from w1 %f w2 %f PSame %f %f f12code %f f12other %f\n",
		      out->weights[i], 1e-6,
		      w1, w2,
		      PSame1[iRate], PSame2[iRate],
		      f12code, f12other);
	    out->weights[i] = 1e-6;
	  }
	  continue;
	}
      }
      /* if we did not compute a simple profile, then do the full computation and
         store the full vector
      */
      if (f1 == NULL) {
	for (j = 0; j < 4; j++)
	  mix1[j] = (1-w1)*0.25;
	if(code1 != NOCODE)
	  mix1[code1] += w1;
	f1 = mix1;
      }
      if (f2 == NULL) {
	for (j = 0; j < 4; j++)
	  mix2[j] = (1-w2)*0.25;
	if(code2 != NOCODE)
	  mix2[code2] += w2;
	f2 = mix2;
      }
      out->codes[i] = NOCODE;
      out->weights[i] = 1.0;
      numeric_t *f = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      double lkAB = 0;
      for (j = 0; j < 4; j++) {
	f[j] = (f1[j] * PSame1[iRate] + (1.0-f1[j]) * PDiff1[iRate])
	  * (f2[j] * PSame2[iRate] + (1.0-f2[j]) * PDiff2[iRate]);
	lkAB += f[j];
      }
      double lkABInv = 1.0/lkAB;
      for (j = 0; j < 4; j++)
	f[j] *= lkABInv;
    }
    PSame1 = myfree(PSame1, sizeof(double) * rates->nRateCategories);
    PSame2 = myfree(PSame2, sizeof(double) * rates->nRateCategories);
    PDiff1 = myfree(PDiff1, sizeof(double) * rates->nRateCategories);
    PDiff2 = myfree(PDiff2, sizeof(double) * rates->nRateCategories);
  } else if (nCodes == 4) {	/* matrix model on nucleotides */
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];
    numeric_t f1mix[4], f2mix[4];
    
    for (i=0; i < nPos; i++) {
      if (p1->codes[i] == NOCODE && p2->codes[i] == NOCODE
	  && p1->weights[i] == 0 && p2->weights[i] == 0) {
	/* aligning gap with gap -- just output a gap
	   out->codes[i] is already set to NOCODE so need not set that */
	out->weights[i] = 0;
	continue;
      }
      int iRate = rates->ratecat[i];
      numeric_t *expeigen1 = &expeigenRates1[iRate*4];
      numeric_t *expeigen2 = &expeigenRates2[iRate*4];
      numeric_t *f1 = GET_FREQ(p1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(p2,i,/*IN/OUT*/iFreq2);
      numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      assert(fOut != NULL);

      if (f1 == NULL) {
	f1 = &transmat->codeFreq[p1->codes[i]][0]; /* codeFreq includes an entry for NOCODE */
	double w = p1->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 4; j++)
	    f1mix[j] = w * f1[j] + (1.0-w) * fGap[j];
	  f1 = f1mix;
	}
      }
      if (f2 == NULL) {
	f2 = &transmat->codeFreq[p2->codes[i]][0];
	double w = p2->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 4; j++)
	    f2mix[j] = w * f2[j] + (1.0-w) * fGap[j];
	  f2 = f2mix;
	}
      }
      numeric_t fMult1[4] ALIGNED;	/* rotated1 * expeigen1 */
      numeric_t fMult2[4] ALIGNED;	/* rotated2 * expeigen2 */
#if 0 /* SSE3 is slower */
      vector_multiply(f1, expeigen1, 4, /*OUT*/fMult1);
      vector_multiply(f2, expeigen2, 4, /*OUT*/fMult2);
#else
      for (j = 0; j < 4; j++) {
	fMult1[j] = f1[j]*expeigen1[j];
	fMult2[j] = f2[j]*expeigen2[j];
      }
#endif
      numeric_t fPost[4] ALIGNED;		/* in  unrotated space */
      for (j = 0; j < 4; j++) {
#if 0 /* SSE3 is slower */
	fPost[j] = vector_dot_product_rot(fMult1, fMult2, &transmat->codeFreq[j][0], 4)
	  * transmat->statinv[j]; */
#else
	double out1 = 0;
	double out2 = 0;
	for (k = 0; k < 4; k++) {
	  out1 += fMult1[k] * transmat->codeFreq[j][k];
	  out2 += fMult2[k] * transmat->codeFreq[j][k];
	}
	fPost[j] = out1*out2*transmat->statinv[j];
#endif
      }
      double fPostTot = 0;
      for (j = 0; j < 4; j++)
	fPostTot += fPost[j];
      assert(fPostTot > fPostTotalTolerance);
      double fPostInv = 1.0/fPostTot;
#if 0 /* SSE3 is slower */
      vector_multiply_by(fPost, fPostInv, 4);
#else
      for (j = 0; j < 4; j++)
	fPost[j] *= fPostInv;
#endif

      /* and finally, divide by stat again & rotate to give the new frequencies */
      matrixt_by_vector4(transmat->eigeninvT, fPost, /*OUT*/fOut);
    }  /* end loop over position i */
  } else if (nCodes == 20) {	/* matrix model on amino acids */
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];
    numeric_t f1mix[20] ALIGNED;
    numeric_t f2mix[20] ALIGNED;
    
    for (i=0; i < nPos; i++) {
      if (p1->codes[i] == NOCODE && p2->codes[i] == NOCODE
	  && p1->weights[i] == 0 && p2->weights[i] == 0) {
	/* aligning gap with gap -- just output a gap
	   out->codes[i] is already set to NOCODE so need not set that */
	out->weights[i] = 0;
	continue;
      }
      int iRate = rates->ratecat[i];
      numeric_t *expeigen1 = &expeigenRates1[iRate*20];
      numeric_t *expeigen2 = &expeigenRates2[iRate*20];
      numeric_t *f1 = GET_FREQ(p1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(p2,i,/*IN/OUT*/iFreq2);
      numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      assert(fOut != NULL);

      if (f1 == NULL) {
	f1 = &transmat->codeFreq[p1->codes[i]][0]; /* codeFreq includes an entry for NOCODE */
	double w = p1->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 20; j++)
	    f1mix[j] = w * f1[j] + (1.0-w) * fGap[j];
	  f1 = f1mix;
	}
      }
      if (f2 == NULL) {
	f2 = &transmat->codeFreq[p2->codes[i]][0];
	double w = p2->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 20; j++)
	    f2mix[j] = w * f2[j] + (1.0-w) * fGap[j];
	  f2 = f2mix;
	}
      }
      numeric_t fMult1[20] ALIGNED;	/* rotated1 * expeigen1 */
      numeric_t fMult2[20] ALIGNED;	/* rotated2 * expeigen2 */
      vector_multiply(f1, expeigen1, 20, /*OUT*/fMult1);
      vector_multiply(f2, expeigen2, 20, /*OUT*/fMult2);
      numeric_t fPost[20] ALIGNED;		/* in  unrotated space */
      for (j = 0; j < 20; j++) {
	numeric_t value = vector_dot_product_rot(fMult1, fMult2, &transmat->codeFreq[j][0], 20)
	  * transmat->statinv[j];
	/* Added this logic try to avoid rare numerical problems */
	fPost[j] = value >= 0 ? value : 0;
      }
      double fPostTot = vector_sum(fPost, 20);
      assert(fPostTot > fPostTotalTolerance);
      double fPostInv = 1.0/fPostTot;
      vector_multiply_by(/*IN/OUT*/fPost, fPostInv, 20);
      int ch = -1;		/* the dominant character, if any */
      if (!exactML) {
	for (j = 0; j < 20; j++) {
	  if (fPost[j] >= approxMLminf) {
	    ch = j;
	    break;
	  }
	}
      }

      /* now, see if we can use the approximation 
	 fPost ~= (1 or 0) * w + nearP * (1-w)
	 to avoid rotating */
      double w = 0;
      if (ch >= 0) {
	w = (fPost[ch] - transmat->nearP[ch][ch]) / (1.0 - transmat->nearP[ch][ch]);
	for (j = 0; j < 20; j++) {
	  if (j != ch) {
	    double fRough = (1.0-w) * transmat->nearP[ch][j];
	    if (fRough < fPost[j]  * approxMLminratio) {
	      ch = -1;		/* give up on the approximation */
	      break;
	    }
	  }
	}
      }
      if (ch >= 0) {
	nAAPosteriorRough++;
	double wInvStat = w * transmat->statinv[ch];
	for (j = 0; j < 20; j++)
	  fOut[j] = wInvStat * transmat->codeFreq[ch][j] + (1.0-w) * transmat->nearFreq[ch][j];
      } else {
	/* and finally, divide by stat again & rotate to give the new frequencies */
	nAAPosteriorExact++;
	for (j = 0; j < 20; j++)
	  fOut[j] = vector_multiply_sum(fPost, &transmat->eigeninv[j][0], 20);
      }
    } /* end loop over position i */
  } else {
    assert(0);			/* illegal nCodes */
  }

  if (transmat != NULL) {
    expeigenRates1 = myfree(expeigenRates1, sizeof(numeric_t) * rates->nRateCategories * nCodes);
    expeigenRates2 = myfree(expeigenRates2, sizeof(numeric_t) * rates->nRateCategories * nCodes);
  }

  /* Reallocate out->vectors to be the right size */
  out->nVectors = iFreqOut;
  if (out->nVectors == 0)
    out->vectors = (numeric_t*)myfree(out->vectors, sizeof(numeric_t)*nCodes*nPos);
  else
    out->vectors = (numeric_t*)myrealloc(out->vectors,
				     /*OLDSIZE*/sizeof(numeric_t)*nCodes*nPos,
				     /*NEWSIZE*/sizeof(numeric_t)*nCodes*out->nVectors,
				     /*copy*/true); /* try to save space */
  nProfileFreqAlloc += out->nVectors;
  nProfileFreqAvoid += nPos - out->nVectors;

  /* compute total constraints */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] = p1->nOn[i] + p2->nOn[i];
    out->nOff[i] = p1->nOff[i] + p2->nOff[i];
  }
  nPosteriorCompute++;
  return(out);
}