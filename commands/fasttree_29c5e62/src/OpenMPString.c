#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

char *OpenMPString(void) {
#ifdef OPENMP
  static char buf[100];
  sprintf(buf, ", OpenMP (%d threads)", omp_get_max_threads());
  return(buf);
#else
  return("");
#endif
}