#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


void LogMLRates(/*OPTIONAL WRITE*/FILE *fpLog, NJ_t *NJ) {
  if (fpLog != NULL) {
    rates_t *rates = &NJ->rates;
    fprintf(fpLog, "NCategories\t%d\nRates",rates->nRateCategories);
    assert(rates->nRateCategories > 0);
    int iRate;
    for (iRate = 0; iRate < rates->nRateCategories; iRate++)
      fprintf(fpLog, " %f", rates->rates[iRate]);
    fprintf(fpLog,"\nSiteCategories");
    int iPos;
    for (iPos = 0; iPos < NJ->nPos; iPos++) {
      iRate = rates->ratecat[iPos];
      fprintf(fpLog," %d",iRate+1);
    }
    fprintf(fpLog,"\n");
    fflush(fpLog);
  }
}

void LogTree(char *format, int i, /*OPTIONAL WRITE*/FILE *fpLog, NJ_t *NJ, char **names, uniquify_t *unique, bool bQuote) {
  if(fpLog != NULL) {
    fprintf(fpLog, format, i);
    fprintf(fpLog, "\t");
    PrintNJ(fpLog, NJ, names, unique, /*support*/false, bQuote);
    fflush(fpLog);
  }
}

void ReadTreeError(char *err, char *token) {
  fprintf(stderr, "Tree parse error: unexpected token '%s' -- %s\n",
	  token == NULL ? "(End of file)" : token,
	  err);
  exit(1);
}

double LogCorrect(double dist) {
  const double maxscore = 3.0;
  if (nCodes == 4 && !useMatrix) { /* Jukes-Cantor */
    dist = dist < 0.74 ? -0.75*log(1.0 - dist * 4.0/3.0) : maxscore;
  } else {			/* scoredist-like */
    dist = dist < 0.99 ? -1.3*log(1.0 - dist) : maxscore;
  }
  return (dist < maxscore ? dist : maxscore);
}