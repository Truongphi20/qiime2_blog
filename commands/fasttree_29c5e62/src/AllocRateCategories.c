#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

/* Allocate or reallocate the rate categories, and set every position
   to category 0 and every category's rate to 1.0
   If nRateCategories=0, just deallocate
*/
void AllocRateCategories(/*IN/OUT*/rates_t *rates, int nRateCategories, int nPos) {
  assert(nRateCategories >= 0);
  rates->rates = myfree(rates->rates, sizeof(numeric_t)*rates->nRateCategories);
  rates->ratecat = myfree(rates->ratecat, sizeof(unsigned int)*nPos);
  rates->nRateCategories = nRateCategories;
  if (rates->nRateCategories > 0) {
    rates->rates = (numeric_t*)mymalloc(sizeof(numeric_t)*rates->nRateCategories);
    int i;
    for (i = 0; i < nRateCategories; i++)
      rates->rates[i] = 1.0;
    rates->ratecat = (unsigned int *)mymalloc(sizeof(unsigned int)*nPos);
    for (i = 0; i < nPos; i++)
      rates->ratecat[i] = 0;
  }
}