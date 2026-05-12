#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Pick columns for resampling, stored as returned_vector[iBoot*nPos + j] */
int *ResampleColumns(int nPos, int nBootstrap) {
  long lPos = nPos; /* to prevent overflow on very long alignments when multiplying nPos * nBootstrap */
  int *col = (int*)mymalloc(sizeof(int)*lPos*(size_t)nBootstrap);
  int i;
  for (i = 0; i < nBootstrap; i++) {
    int j;
    for (j = 0; j < nPos; j++) {
      int pos   = (int)(knuth_rand() * nPos);
      if (pos<0)
	pos = 0;
      else if (pos == nPos)
	pos = nPos-1;
      col[i*lPos + j] = pos;
    }
  }
  if (verbose > 5) {
    for (i=0; i < 3 && i < nBootstrap; i++) {
      fprintf(stderr,"Boot%d",i);
      int j;
      for (j = 0; j < nPos; j++) {
	fprintf(stderr,"\t%d",col[i*lPos+j]);
      }
      fprintf(stderr,"\n");
    }
  }
  return(col);
}