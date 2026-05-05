#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


/* Make the (unrotated) frequencies sum to 1
   Simply dividing by total_weight is not ideal because of roundoff error
   So compute total_freq instead
*/
void NormalizeFreq(/*IN/OUT*/numeric_t *freq, distance_matrix_t *dmat) {
  double total_freq = 0;
  int k;
  if (dmat != NULL) {
    /* The total frequency is dot_product(true_frequencies, 1)
       So we rotate the 1 vector by eigeninv (stored in eigentot)
    */
    total_freq = vector_multiply_sum(freq, dmat->eigentot, nCodes);
  } else {
    for (k = 0; k < nCodes; k++)
      total_freq += freq[k];
  }
  if (total_freq > fPostTotalTolerance) {
    numeric_t inverse_weight = 1.0/total_freq;
    vector_multiply_by(/*IN/OUT*/freq, inverse_weight, nCodes);
  } else {
    /* This can happen if we are in a very low-weight region, e.g. if a mostly-gap position gets weighted down
       repeatedly; just set them all to arbitrary but legal values */
    if (dmat == NULL) {
      for (k = 0; k < nCodes; k++)
	freq[k] = 1.0/nCodes;
    } else {
      for (k = 0; k < nCodes; k++)
	freq[k] = dmat->codeFreq[0][k];
    }
  }
}