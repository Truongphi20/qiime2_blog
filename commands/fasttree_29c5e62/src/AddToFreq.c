#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* This should not be called if the update weight is 0, as
   in that case code==NOCODE and in=NULL is possible, and then
   it will fail.
*/
void AddToFreq(/*IN/OUT*/numeric_t *fOut,
	       double weight,
	       unsigned int codeIn, /*OPTIONAL*/numeric_t *fIn,
	       /*OPTIONAL*/distance_matrix_t *dmat) {
  assert(fOut != NULL);
  if (fIn != NULL) {
    vector_add_mult(fOut, fIn, weight, nCodes);
  } else if (dmat) {
    assert(codeIn != NOCODE);
    vector_add_mult(fOut, dmat->codeFreq[codeIn], weight, nCodes);
  } else {
    assert(codeIn != NOCODE);
    fOut[codeIn] += weight;
  }
}