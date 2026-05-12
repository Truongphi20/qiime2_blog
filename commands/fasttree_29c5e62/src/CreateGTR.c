#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


transition_matrix_t *CreateGTR(double *r/*ac ag at cg ct gt*/, double *f/*acgt*/) {
  double matrix[MAXCODES][MAXCODES];
  assert(nCodes==4);
  int i, j;
  /* Place rates onto a symmetric matrix, but correct by f(target), so that
     stationary distribution f[] is maintained
     Leave diagonals as 0 (CreateTransitionMatrix will fix them)
  */
  int imat = 0;
  for (i = 0; i < nCodes; i++) {
    matrix[i][i] = 0;
    for (j = i+1; j < nCodes; j++) {
      double rate = r[imat++];
      assert(rate > 0);
      /* Want t(matrix) * f to be 0 */
      matrix[i][j] = rate * f[i];
      matrix[j][i] = rate * f[j];
    }
  }
  /* Compute average mutation rate */
  double total_rate = 0;
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      total_rate += f[i] * matrix[i][j];
  assert(total_rate > 1e-6);
  double inv = 1.0/total_rate;
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      matrix[i][j] *= inv;
  return(CreateTransitionMatrix(matrix,f));
}