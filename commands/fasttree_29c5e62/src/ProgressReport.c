#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void ProgressReport(char *format, int i1, int i2, int i3, int i4) {
  static bool time_set = false;
  static struct timeval time_last;
  static struct timeval time_begin;

  if (!showProgress)
    return;

  static struct timeval time_now;
  gettimeofday(&time_now,NULL);
  if (!time_set) {
    time_begin = time_last = time_now;
    time_set = true;
  }
  static struct timeval elapsed;
  timeval_subtract(&elapsed,&time_now,&time_last);
  
  if (elapsed.tv_sec > 1 || elapsed.tv_usec > 100*1000 || verbose > 1) {
    timeval_subtract(&elapsed,&time_now,&time_begin);
    fprintf(stderr, "%7i.%2.2i seconds: ", (int)elapsed.tv_sec, (int)(elapsed.tv_usec/10000));
    fprintf(stderr, format, i1, i2, i3, i4);
    if (verbose > 1 || !isatty(STDERR_FILENO)) {
      fprintf(stderr, "\n");
    } else {
      fprintf(stderr, "   \r");
    }
    fflush(stderr);
    time_last = time_now;
  }
}