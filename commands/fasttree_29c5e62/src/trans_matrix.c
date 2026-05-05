#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

transition_matrix_t *CreateTransitionMatrix(/*IN*/double matrix[MAXCODES][MAXCODES],
					    /*IN*/double stat[MAXCODES]) {
  int i,j,k;
  transition_matrix_t *transmat = mymalloc(sizeof(transition_matrix_t));
  double sqrtstat[20];
  for (i = 0; i < nCodes; i++) {
    transmat->stat[i] = stat[i];
    transmat->statinv[i] = 1.0/stat[i];
    sqrtstat[i] = sqrt(stat[i]);
  }

  double sym[20*20];		/* symmetrized matrix M' */
  /* set diagonals so columns sums are 0 before symmetrization */
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      sym[nCodes*i+j] = matrix[i][j];
  for (j = 0; j < nCodes; j++) {
    double sum = 0;
    sym[nCodes*j+j] = 0;
    for (i = 0; i < nCodes; i++)
      sum += sym[nCodes*i+j];
    sym[nCodes*j+j] = -sum;
  }
  /* M' = S**-1 M S */
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      sym[nCodes*i+j] *= sqrtstat[j]/sqrtstat[i];

  /* eigen decomposition of M' -- note that eigenW is the transpose of what we want,
     which is eigenvectors in columns */
  double eigenW[20*20], eval[20], e[20];
  for (i = 0; i < nCodes*nCodes; i++)
    eigenW[i] = sym[i];
  tred2(eigenW, nCodes, nCodes, eval, e);       
  tqli(eval, e, nCodes , nCodes, eigenW);

  /* save eigenvalues */
  for (i = 0; i < nCodes; i++)
    transmat->eigenval[i] = eval[i];

  /* compute eigen decomposition of M into t(codeFreq): V = S*W */
  /* compute inverse of V in eigeninv: V**-1 = t(W) S**-1  */
  for (i = 0; i < nCodes; i++) {
    for (j = 0; j < nCodes; j++) {
      transmat->eigeninv[i][j] = eigenW[nCodes*i+j] / sqrtstat[j];
      transmat->eigeninvT[j][i] = transmat->eigeninv[i][j];
    }
  }
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      transmat->codeFreq[i][j] = eigenW[j*nCodes+i] * sqrtstat[i];
  /* codeFreq[NOCODE] is the rotation of (1,1,...) not (1/nCodes,1/nCodes,...), which
     gives correct posterior probabilities
  */
  for (j = 0; j < nCodes; j++) {
    transmat->codeFreq[NOCODE][j] = 0.0;
    for (i = 0; i < nCodes; i++)
      transmat->codeFreq[NOCODE][j] += transmat->codeFreq[i][j];
  }
  /* save some posterior probabilities for approximating later:
     first, we compute P(B | A, t) for t = approxMLnearT, by using
     V * exp(L*t) * V**-1 */
  double expvalues[MAXCODES];
  for (i = 0; i < nCodes; i++)
    expvalues[i] = exp(approxMLnearT * transmat->eigenval[i]);
  double LVinv[MAXCODES][MAXCODES]; /* exp(L*t) * V**-1 */
  for (i = 0; i < nCodes; i++) {
    for (j = 0; j < nCodes; j++)
      LVinv[i][j] = transmat->eigeninv[i][j] * expvalues[i];
  }
  /* matrix transform for converting A -> B given t: transt[i][j] = P(j->i | t) */
  double transt[MAXCODES][MAXCODES];
  for (i = 0; i < nCodes; i++) {
    for (j = 0; j < nCodes; j++) {
      transt[i][j] = 0;
      for (k = 0; k < nCodes; k++)
	transt[i][j] += transmat->codeFreq[i][k] * LVinv[k][j];
    }
  }
  /* nearP[i][j] = P(parent = j | both children are i) = P(j | i,i) ~ stat(j) * P(j->i | t)**2 */
  for (i = 0; i < nCodes; i++) {
    double nearP[MAXCODES];
    double tot = 0;
    for (j = 0; j < nCodes; j++) {
      assert(transt[j][i] > 0);
      assert(transmat->stat[j] > 0);
      nearP[j] = transmat->stat[j] * transt[i][j] * transt[i][j];
      tot += nearP[j];
    }
    assert(tot > 0);
    for (j = 0; j < nCodes; j++)
      nearP[j] *= 1.0/tot;
    /* save nearP in transmat->nearP[i][] */
    for (j = 0; j < nCodes; j++)
      transmat->nearP[i][j] = nearP[j];
    /* multiply by 1/stat and rotate nearP */
    for (j = 0; j < nCodes; j++)
      nearP[j] /= transmat->stat[j];
    for (j = 0; j < nCodes; j++) {
      double rot = 0;
      for (k = 0; k < nCodes; k++)
	rot += nearP[k] * transmat->codeFreq[i][j];
      transmat->nearFreq[i][j] = rot;
    }
  }
  return(transmat);
  assert(0);
}

distance_matrix_t *TransMatToDistanceMat(transition_matrix_t *transmat) {
  if (transmat == NULL)
    return(NULL);
  distance_matrix_t *dmat = mymalloc(sizeof(distance_matrix_t));
  int i, j;
  for (i=0; i<nCodes; i++) {
    for (j=0; j<nCodes; j++) {
      dmat->distances[i][j] = 0;	/* never actually used */
      dmat->eigeninv[i][j] = transmat->eigeninv[i][j];
      dmat->codeFreq[i][j] = transmat->codeFreq[i][j];
    }
  }
  /* eigentot . rotated-vector is the total frequency of the unrotated vector
     (used to normalize in NormalizeFreq()
     For transition matrices, we rotate by transpose of eigenvectors, so
     we need to multiply by the inverse matrix by 1....1 to get this vector,
     or in other words, sum the columns
  */
  for(i = 0; i<nCodes; i++) {
      dmat->eigentot[i] = 0.0;
      for (j = 0; j<nCodes; j++)
	dmat->eigentot[i] += transmat->eigeninv[i][j];
  }
  return(dmat);
}