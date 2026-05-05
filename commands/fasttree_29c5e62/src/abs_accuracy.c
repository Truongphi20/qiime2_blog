#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Algorithm 26.2.17 from Abromowitz and Stegun, Handbook of Mathematical Functions
   Absolute accuracy of only about 1e-7, which is enough for us
*/
double pnorm(double x)
{
  double b1 =  0.319381530;
  double b2 = -0.356563782;
  double b3 =  1.781477937;
  double b4 = -1.821255978;
  double b5 =  1.330274429;
  double p  =  0.2316419;
  double c  =  0.39894228;

  if(x >= 0.0) {
    double t = 1.0 / ( 1.0 + p * x );
    return (1.0 - c * exp( -x * x / 2.0 ) * t *
	    ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  /*else*/
  double t = 1.0 / ( 1.0 - p * x );
  return ( c * exp( -x * x / 2.0 ) * t *
	   ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
}

void *mymalloc(size_t sz) {
  if (sz == 0) return(NULL);
  void *new = malloc(sz);
  if (new == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }
  szAllAlloc += sz;
  mymallocUsed += sz;
#ifdef TRACK_MEMORY
  struct mallinfo mi = mallinfo();
  if (mi.arena+mi.hblkhd > maxmallocHeap)
    maxmallocHeap = mi.arena+mi.hblkhd;
#endif
  /* gcc malloc should always return 16-byte-aligned values... */
  assert(IS_ALIGNED(new));
  return (new);
}

void *mymemdup(void *data, size_t sz) {
  if(data==NULL) return(NULL);
  void *new = mymalloc(sz);
  memcpy(/*to*/new, /*from*/data, sz);
  return(new);
}

void *myrealloc(void *data, size_t szOld, size_t szNew, bool bCopy) {
  if (data == NULL && szOld == 0)
    return(mymalloc(szNew));
  if (data == NULL || szOld == 0 || szNew == 0) {
    fprintf(stderr,"Empty myrealloc\n");
    exit(1);
  }
  if (szOld == szNew)
    return(data);
  void *new = NULL;
  if (bCopy) {
    /* Try to reduce memory fragmentation by allocating anew and copying
       Seems to help in practice */
    new = mymemdup(data, szNew);
    myfree(data, szOld);
  } else {
    new = realloc(data,szNew);
    if (new == NULL) {
      fprintf(stderr, "Out of memory\n");
      exit(1);
    }
    assert(IS_ALIGNED(new));
    szAllAlloc += (szNew-szOld);
    mymallocUsed += (szNew-szOld);
#ifdef TRACK_MEMORY
    struct mallinfo mi = mallinfo();
    if (mi.arena+mi.hblkhd > maxmallocHeap)
      maxmallocHeap = mi.arena+mi.hblkhd;
#endif
  }
  return(new);
}

void *myfree(void *p, size_t sz) {
  if(p==NULL) return(NULL);
  free(p);
  mymallocUsed -= sz;
  return(NULL);
}