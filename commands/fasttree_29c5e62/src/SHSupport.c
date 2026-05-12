#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

double SHSupport(int nPos, int nBootstrap, int *col, double loglk[3], double *site_likelihoods[3]) {
  long lPos = nPos;		/* to avoid overflow when multiplying */
  assert(nBootstrap>0);
  double delta1 = loglk[0]-loglk[1];
  double delta2 = loglk[0]-loglk[2];
  double delta = delta1 < delta2 ? delta1 : delta2;

  double *siteloglk[3];
  int i,j;
  for (i = 0; i < 3; i++) {
    siteloglk[i] = mymalloc(sizeof(double)*nPos);
    for (j = 0; j < nPos; j++)
      siteloglk[i][j] = log(site_likelihoods[i][j]);
  }

  int nSupport = 0;
  int iBoot;
  for (iBoot = 0; iBoot < nBootstrap; iBoot++) {
    double resampled[3];
    for (i = 0; i < 3; i++)
      resampled[i] = -loglk[i];
    for (j = 0; j < nPos; j++) {
      int pos = col[iBoot*lPos+j];
      for (i = 0; i < 3; i++)
	resampled[i] += siteloglk[i][pos];
    }
    int iBest = 0;
    for (i = 1; i < 3; i++)
      if (resampled[i] > resampled[iBest])
	iBest = i;
    double resample1 = resampled[iBest] - resampled[(iBest+1)%3];
    double resample2 = resampled[iBest] - resampled[(iBest+2)%3];
    double resampleDelta = resample1 < resample2 ? resample1 : resample2;
    if (resampleDelta < delta)
      nSupport++;
  }
  for (i=0;i<3;i++)
    siteloglk[i] = myfree(siteloglk[i], sizeof(double)*nPos);
  return(nSupport/(double)nBootstrap);
}