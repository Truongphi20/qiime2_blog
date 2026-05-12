#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* The BIONJ-like formula for the weight of A when building a profile for AB is
     1/2 + (avgD(B,CD) - avgD(A,CD))/(2*d(A,B))
*/
double QuartetWeight(profile_t *profiles[4], distance_matrix_t *dmat, int nPos) {
  if (!bionj)
    return(-1.0); /* even weighting */
  double d[6];
  CorrectedPairDistances(profiles, 4, dmat, nPos, /*OUT*/d);
  if (d[qAB] < 0.01)
    return -1.0;
  double weight = 0.5 + ((d[qBC]+d[qBD])-(d[qAC]+d[qAD]))/(4*d[qAB]);
  if (weight < 0)
    weight = 0;
  if (weight > 1)
    weight = 1;
  return (weight);
}