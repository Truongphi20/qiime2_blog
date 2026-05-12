#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


#ifdef USE_SSE3
inline float mm_sum(register __m128 sum) {
#if 1
  /* stupider but faster */
  float f[4] ALIGNED;
  _mm_store_ps(f,sum);
  return(f[0]+f[1]+f[2]+f[3]);
#else
  /* first we get sum[0]+sum[1], sum[2]+sum[3] by selecting 0/1 and 2/3 */
  sum = _mm_add_ps(sum,_mm_shuffle_ps(sum,sum,_MM_SHUFFLE(0,1,2,3)));
  /* then get sum[0]+sum[1]+sum[2]+sum[3] by selecting 0/1 and 0/1 */
  sum = _mm_add_ps(sum,_mm_shuffle_ps(sum,sum,_MM_SHUFFLE(0,1,0,1)));
  float f;
  _mm_store_ss(&f, sum);	/* save the lowest word */
  return(f);
#endif
}
#endif

void vector_multiply(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, int n, /*OUT*/numeric_t *fOut) {
#ifdef USE_SSE3
  int i;
  for (i = 0; i < n; i += 4) {
    __m128 a, b, c;
    a = _mm_load_ps(f1+i);
    b = _mm_load_ps(f2+i);
    c = _mm_mul_ps(a, b);
    _mm_store_ps(fOut+i,c);
  }
#else
  int i;
  for (i = 0; i < n; i++)
    fOut[i] = f1[i]*f2[i];
#endif
}

numeric_t vector_multiply_sum(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, int n) {
#ifdef USE_SSE3
  if (n == 4)
    return(f1[0]*f2[0]+f1[1]*f2[1]+f1[2]*f2[2]+f1[3]*f2[3]);
  __m128 sum = _mm_setzero_ps();
  int i;
  for (i = 0; i < n; i += 4) {
    __m128 a, b, c;
    a = _mm_load_ps(f1+i);
    b = _mm_load_ps(f2+i);
    c = _mm_mul_ps(a, b);
    sum = _mm_add_ps(c, sum);
  }
  return(mm_sum(sum));
#else
  int i;
  numeric_t out = 0.0;
  for (i=0; i < n; i++)
    out += f1[i]*f2[i];
  return(out);
#endif
}

/* sum(f1*f2*f3) */
numeric_t vector_multiply3_sum(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, /*IN*/numeric_t* f3, int n) {
#ifdef USE_SSE3
  __m128 sum = _mm_setzero_ps();
  int i;
  for (i = 0; i < n; i += 4) {
    __m128 a1, a2, a3;
    a1 = _mm_load_ps(f1+i);
    a2 = _mm_load_ps(f2+i);
    a3 = _mm_load_ps(f3+i);
    sum = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(a1,a2),a3),sum);
  }
  return(mm_sum(sum));
#else
  int i;
  numeric_t sum = 0.0;
  for (i = 0; i < n; i++)
    sum += f1[i]*f2[i]*f3[i];
  return(sum);
#endif
}

numeric_t vector_dot_product_rot(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, /*IN*/numeric_t *fBy, int n) {
#ifdef USE_SSE3
  __m128 sum1 = _mm_setzero_ps();
  __m128 sum2 = _mm_setzero_ps();
  int i;
  for (i = 0; i < n; i += 4) {
    __m128 a1, a2, aBy;
    a1 = _mm_load_ps(f1+i);
    a2 = _mm_load_ps(f2+i);
    aBy = _mm_load_ps(fBy+i);
    sum1 = _mm_add_ps(_mm_mul_ps(a1, aBy), sum1);
    sum2 = _mm_add_ps(_mm_mul_ps(a2, aBy), sum2);
  }
  return(mm_sum(sum1)*mm_sum(sum2));
#else
  int i;
  numeric_t out1 = 0.0;
  numeric_t out2 = 0.0;
  for (i=0; i < n; i++) {
    out1 += f1[i]*fBy[i];
    out2 += f2[i]*fBy[i];
  }
  return(out1*out2);
#endif
}

numeric_t vector_sum(/*IN*/numeric_t *f1, int n) {
#ifdef USE_SSE3
  if (n==4)
    return(f1[0]+f1[1]+f1[2]+f1[3]);
  __m128 sum = _mm_setzero_ps();
  int i;
  for (i = 0; i < n; i+=4) {
    __m128 a;
    a = _mm_load_ps(f1+i);
    sum = _mm_add_ps(a, sum);
  }
  return(mm_sum(sum));
#else
  numeric_t out = 0.0;
  int i;
  for (i = 0; i < n; i++)
    out += f1[i];
  return(out);
#endif
}

void vector_multiply_by(/*IN/OUT*/numeric_t *f, /*IN*/numeric_t fBy, int n) {
  int i;
#ifdef USE_SSE3
  __m128 c = _mm_set1_ps(fBy);
  for (i = 0; i < n; i += 4) {
    __m128 a, b;
    a = _mm_load_ps(f+i);
    b = _mm_mul_ps(a,c);
    _mm_store_ps(f+i,b);
  }
#else
  for (i = 0; i < n; i++)
    f[i] *= fBy;
#endif
}

void vector_add_mult(/*IN/OUT*/numeric_t *fTot, /*IN*/numeric_t *fAdd, numeric_t weight, int n) {
#ifdef USE_SSE3
  int i;
  __m128 w = _mm_set1_ps(weight);
  for (i = 0; i < n; i += 4) {
    __m128 tot, add;
    tot = _mm_load_ps(fTot+i);
    add = _mm_load_ps(fAdd+i);
    _mm_store_ps(fTot+i, _mm_add_ps(tot, _mm_mul_ps(add,w)));
  }
#else
  int i;
  for (i = 0; i < n; i++)
    fTot[i] += fAdd[i] * weight;
#endif
}

void matrixt_by_vector4(/*IN*/numeric_t mat[4][MAXCODES], /*IN*/numeric_t vec[4], /*OUT*/numeric_t out[4]) {
#ifdef USE_SSE3
  /*__m128 v = _mm_load_ps(vec);*/
  __m128 o = _mm_setzero_ps();
  int j;
  /* result is a sum of vectors: sum(k) v[k] * mat[k][] */
  for (j = 0; j < 4; j++) {
    __m128 m = _mm_load_ps(&mat[j][0]);
    __m128 vj = _mm_load1_ps(&vec[j]);	/* is it faster to shuffle v? */
    o = _mm_add_ps(o, _mm_mul_ps(vj,m));
  }
  _mm_store_ps(out, o);
#else
  int j,k;
  for (j = 0; j < 4; j++) {
    double sum = 0;
    for (k = 0; k < 4; k++)
      sum += vec[k] * mat[k][j];
    out[j] = sum;
  }
#endif
}