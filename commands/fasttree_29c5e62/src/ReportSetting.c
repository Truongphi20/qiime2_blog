#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void ReportSetting(FastTreeOptions_t opt, int nFPs, FILE *fps[2])
{
    char tophitString[100] = "no";
    char tophitsCloseStr[100] = "default";
    if(tophitsClose > 0) sprintf(tophitsCloseStr,"%.2f",tophitsClose);
    if(tophitsMult>0) sprintf(tophitString,"%.2f*sqrtN close=%s refresh=%.2f",
			      tophitsMult, tophitsCloseStr, tophitsRefresh);
    char supportString[100] = "none";
    if (opt.nBootstrap>0) {
      if (opt.MLnni != 0 || opt.MLlen)
	sprintf(supportString, "SH-like %d", opt.nBootstrap);
      else
	sprintf(supportString,"Local boot %d",opt.nBootstrap);
    }
    char nniString[100] = "(no NNI)";
    if (opt.nni > 0)
      sprintf(nniString, "+NNI (%d rounds)", opt.nni);
    if (opt.nni == -1)
      strcpy(nniString, "+NNI");
    char sprString[100] = "(no SPR)";
    if (opt.spr > 0)
      sprintf(sprString, "+SPR (%d rounds range %d)", opt.spr, opt.maxSPRLength);
    char mlnniString[100] = "(no ML-NNI)";
    if(opt.MLnni > 0)
      sprintf(mlnniString, "+ML-NNI (%d rounds)", opt.MLnni);
    else if (opt.MLnni == -1)
      sprintf(mlnniString, "+ML-NNI");
    else if (opt.MLlen)
      sprintf(mlnniString, "+ML branch lengths");
    if ((opt.MLlen || opt.MLnni != 0) && !exactML)
      strcat(mlnniString, " approx");
    if (opt.MLnni != 0)
      sprintf(mlnniString+strlen(mlnniString), " opt-each=%d",mlAccuracy);

    for (int i = 0; i < nFPs; i++) {
      FILE *fp = fps[i];
      fprintf(fp,"FastTree Version %s %s%s\nAlignment: %s",
	      FT_VERSION, SSE_STRING, OpenMPString(), opt.fileName != NULL ? opt.fileName : "standard input");
      if (opt.nAlign>1)
	fprintf(fp, " (%d alignments)", opt.nAlign);
      fprintf(fp,"\n%s distances: %s Joins: %s Support: %s\n",
	      nCodes == 20 ? "Amino acid" : "Nucleotide",
	      opt.matrixPrefix ? opt.matrixPrefix : (useMatrix? "BLOSUM45"
					     : (nCodes==4 && logdist ? "Jukes-Cantor" : "%different")),
	      bionj ? "weighted" : "balanced" ,
	      supportString);
      if (opt.intreeFile == NULL)
	fprintf(fp, "Search: %s%s %s %s %s\nTopHits: %s\n",
		slow?"Exhaustive (slow)" : (fastest ? "Fastest" : "Normal"),
		useTopHits2nd ? "+2nd" : "",
		nniString, sprString, mlnniString,
		tophitString);
      else
	fprintf(fp, "Start at tree from %s %s %s\n", opt.intreeFile, nniString, sprString);
      
      if (opt.MLnni != 0 || opt.MLlen) {
	fprintf(fp, "ML Model: %s,",
		(nCodes == 4) ? 
                (opt.bUseGtr ? "Generalized Time-Reversible" : "Jukes-Cantor") : 
                (opt.transitionFile ? opt.transitionFile :
                 (opt.bUseLg ? "Le-Gascuel 2008" : (opt.bUseWag ? "Whelan-And-Goldman" : "Jones-Taylor-Thorton"))));
	if (opt.nRateCats == 1)
	  fprintf(fp, " No rate variation across sites");
	else
	  fprintf(fp, " CAT approximation with %d rate categories", opt.nRateCats);
	fprintf(fp, "\n");
	if (nCodes == 4 && opt.bUseGtrRates)
	  fprintf(fp, "GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f\n",
		  opt.gtrrates[0],opt.gtrrates[1],opt.gtrrates[2],opt.gtrrates[3],opt.gtrrates[4],opt.gtrrates[5]);
	if (nCodes == 4 && opt.bUseGtrFreq)
	  fprintf(fp, "GTR frequencies(A C G T) %.4f %.4f %.4f %.4f\n",
		  opt.gtrfreq[0],opt.gtrfreq[1],opt.gtrfreq[2],opt.gtrfreq[3]);
      }
      if (opt.constraintsFile != NULL)
	fprintf(fp, "Constraints: %s Weight: %.3f\n", opt.constraintsFile, constraintWeight);
      if (pseudoWeight > 0)
	fprintf(fp, "Pseudocount weight for comparing sequences with little overlap: %.3lf\n",pseudoWeight);
      fflush(fp);
    }
}