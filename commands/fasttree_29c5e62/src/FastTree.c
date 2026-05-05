#include "include_stuff.h"
#include "datastructs.h"
#include "hyper_parameters.h"
#include "support_functions.h"

int main(int argc, char **argv) {
  FastTreeOptions_t opt;
  InitOptions(&opt);
  ParseCommandLine(argc, argv, &opt);

  int nAlign = opt.nAlign;
  char *matrixPrefix = opt.matrixPrefix;
  char *transitionFile = opt.transitionFile;
  distance_matrix_t *distance_matrix = NULL;
  bool make_matrix = opt.make_matrix;
  char *constraintsFile = opt.constraintsFile;
  char *intreeFile = opt.intreeFile;
  bool intree1 = opt.intree1;
  int nni = opt.nni;
  int spr = opt.spr;
  int maxSPRLength = opt.maxSPRLength;
  int MLnni = opt.MLnni;
  bool MLlen = opt.MLlen;
  int nBootstrap = opt.nBootstrap;
  int nRateCats = opt.nRateCats;
  char *logfile = opt.logfile;
  bool bUseGtr = opt.bUseGtr;
  bool bUseLg = opt.bUseLg;
  bool bUseWag = opt.bUseWag;
  bool bUseGtrRates = opt.bUseGtrRates;
  double *gtrrates = opt.gtrrates;
  bool bUseGtrFreq = opt.bUseGtrFreq;
  double *gtrfreq = opt.gtrfreq;
  bool bQuote = opt.bQuote;
  FILE *fpOut = opt.fpOut;

  codesString = nCodes == 20 ? codesStringAA : codesStringNT;
  if (nCodes == 4 && matrixPrefix == NULL)
    useMatrix = false; 		/* no default nucleotide matrix */
  if (transitionFile && nCodes != 20) {
    fprintf(stderr, "The -trans option is only supported for amino acid alignments\n");
    exit(1);
  }
#ifndef USE_DOUBLE
  if (transitionFile)
    fprintf(stderr,
            "Warning: custom matrices may create numerical problems for single-precision FastTree.\n"
            "You may want to recompile with -DUSE_DOUBLE\n");
#endif

  char *fileName = opt.fileName;

  if (slow && fastest) {
    fprintf(stderr,"Cannot be both slow and fastest\n");
    exit(1);
  }
  if (slow && tophitsMult > 0) {
    tophitsMult = 0.0;
  }

  FILE *fpLog = NULL;
  if (logfile != NULL) {
    fpLog = fopen(logfile, "w");
    if (fpLog == NULL) {
      fprintf(stderr, "Cannot write to: %s\n", logfile);
      exit(1);
    }
    fprintf(fpLog, "Command:");
    int i;
    for (i=0; i < argc; i++)
      fprintf(fpLog, " %s", argv[i]);
    fprintf(fpLog,"\n");
    fflush(fpLog);
  }

    int i;
  FILE *fps[2] = {NULL,NULL};
  int nFPs = 0;
  if (verbose)
    fps[nFPs++] = stderr;
  if (fpLog != NULL)
    fps[nFPs++] = fpLog;
  
  if (!make_matrix) {		/* Report settings */
    char tophitString[100] = "no";
    char tophitsCloseStr[100] = "default";
    if(tophitsClose > 0) sprintf(tophitsCloseStr,"%.2f",tophitsClose);
    if(tophitsMult>0) sprintf(tophitString,"%.2f*sqrtN close=%s refresh=%.2f",
			      tophitsMult, tophitsCloseStr, tophitsRefresh);
    char supportString[100] = "none";
    if (nBootstrap>0) {
      if (MLnni != 0 || MLlen)
	sprintf(supportString, "SH-like %d", nBootstrap);
      else
	sprintf(supportString,"Local boot %d",nBootstrap);
    }
    char nniString[100] = "(no NNI)";
    if (nni > 0)
      sprintf(nniString, "+NNI (%d rounds)", nni);
    if (nni == -1)
      strcpy(nniString, "+NNI");
    char sprString[100] = "(no SPR)";
    if (spr > 0)
      sprintf(sprString, "+SPR (%d rounds range %d)", spr, maxSPRLength);
    char mlnniString[100] = "(no ML-NNI)";
    if(MLnni > 0)
      sprintf(mlnniString, "+ML-NNI (%d rounds)", MLnni);
    else if (MLnni == -1)
      sprintf(mlnniString, "+ML-NNI");
    else if (MLlen)
      sprintf(mlnniString, "+ML branch lengths");
    if ((MLlen || MLnni != 0) && !exactML)
      strcat(mlnniString, " approx");
    if (MLnni != 0)
      sprintf(mlnniString+strlen(mlnniString), " opt-each=%d",mlAccuracy);

    for (i = 0; i < nFPs; i++) {
      FILE *fp = fps[i];
      fprintf(fp,"FastTree Version %s %s%s\nAlignment: %s",
	      FT_VERSION, SSE_STRING, OpenMPString(), fileName != NULL ? fileName : "standard input");
      if (nAlign>1)
	fprintf(fp, " (%d alignments)", nAlign);
      fprintf(fp,"\n%s distances: %s Joins: %s Support: %s\n",
	      nCodes == 20 ? "Amino acid" : "Nucleotide",
	      matrixPrefix ? matrixPrefix : (useMatrix? "BLOSUM45"
					     : (nCodes==4 && logdist ? "Jukes-Cantor" : "%different")),
	      bionj ? "weighted" : "balanced" ,
	      supportString);
      if (intreeFile == NULL)
	fprintf(fp, "Search: %s%s %s %s %s\nTopHits: %s\n",
		slow?"Exhaustive (slow)" : (fastest ? "Fastest" : "Normal"),
		useTopHits2nd ? "+2nd" : "",
		nniString, sprString, mlnniString,
		tophitString);
      else
	fprintf(fp, "Start at tree from %s %s %s\n", intreeFile, nniString, sprString);
      
      if (MLnni != 0 || MLlen) {
	fprintf(fp, "ML Model: %s,",
		(nCodes == 4) ? 
                (bUseGtr ? "Generalized Time-Reversible" : "Jukes-Cantor") : 
                (transitionFile ? transitionFile :
                 (bUseLg ? "Le-Gascuel 2008" : (bUseWag ? "Whelan-And-Goldman" : "Jones-Taylor-Thorton"))));
	if (nRateCats == 1)
	  fprintf(fp, " No rate variation across sites");
	else
	  fprintf(fp, " CAT approximation with %d rate categories", nRateCats);
	fprintf(fp, "\n");
	if (nCodes == 4 && bUseGtrRates)
	  fprintf(fp, "GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f\n",
		  gtrrates[0],gtrrates[1],gtrrates[2],gtrrates[3],gtrrates[4],gtrrates[5]);
	if (nCodes == 4 && bUseGtrFreq)
	  fprintf(fp, "GTR frequencies(A C G T) %.4f %.4f %.4f %.4f\n",
		  gtrfreq[0],gtrfreq[1],gtrfreq[2],gtrfreq[3]);
      }
      if (constraintsFile != NULL)
	fprintf(fp, "Constraints: %s Weight: %.3f\n", constraintsFile, constraintWeight);
      if (pseudoWeight > 0)
	fprintf(fp, "Pseudocount weight for comparing sequences with little overlap: %.3lf\n",pseudoWeight);
      fflush(fp);
    }
  }
  if (matrixPrefix != NULL) {
    if (!useMatrix) {
      fprintf(stderr,"Cannot use both -matrix and -nomatrix arguments!");
      exit(1);
    }
    distance_matrix = ReadDistanceMatrix(matrixPrefix);
  } else if (useMatrix) { 	/* use default matrix */
    assert(nCodes==20);
    distance_matrix = &matrixBLOSUM45;
    SetupDistanceMatrix(distance_matrix);
  } else {
    distance_matrix = NULL;
  }

  int iAln;
  FILE *fpIn = fileName != NULL ? fopen(fileName, "r") : stdin;
  if (fpIn == NULL) {
    fprintf(stderr, "Cannot read %s\n", fileName);
    exit(1);
  }
  FILE *fpConstraints = NULL;
  if (constraintsFile != NULL) {
    fpConstraints = fopen(constraintsFile, "r");
    if (fpConstraints == NULL) {
      fprintf(stderr, "Cannot read %s\n", constraintsFile);
      exit(1);
    }
  }

  FILE *fpInTree = NULL;
  if (intreeFile != NULL) {
    fpInTree = fopen(intreeFile,"r");
    if (fpInTree == NULL) {
      fprintf(stderr, "Cannot read %s\n", intreeFile);
      exit(1);
    }
  }

  for(iAln = 0; iAln < nAlign; iAln++) {
    alignment_t *aln = ReadAlignment(fpIn, bQuote);
    if (aln->nSeq < 1) {
      fprintf(stderr, "No alignment sequences\n");
      exit(1);
    }
    if (fpLog) {
      fprintf(fpLog, "Read %d sequences, %d positions\n", aln->nSeq, aln->nPos);
      fflush(fpLog);
    }

    struct timeval clock_start;
    gettimeofday(&clock_start,NULL);
    ProgressReport("Read alignment",0,0,0,0);

    /* Check that all names in alignment are unique */
    hashstrings_t *hashnames = MakeHashtable(aln->names, aln->nSeq);
    int i;
    for (i=0; i<aln->nSeq; i++) {
      hashiterator_t hi = FindMatch(hashnames,aln->names[i]);
      if (HashCount(hashnames,hi) != 1) {
	fprintf(stderr,"Non-unique name '%s' in the alignment\n",aln->names[i]);
	exit(1);
      }
    }

    /* Make a list of unique sequences -- note some lists are bigger than required */
    ProgressReport("Hashed the names",0,0,0,0);
    if (make_matrix) {
      NJ_t *NJ = InitNJ(aln->seqs, aln->nSeq, aln->nPos,
			/*constraintSeqs*/NULL, /*nConstraints*/0,
			distance_matrix, /*transmat*/NULL);
      printf("   %d\n",aln->nSeq);
      int i,j;
      for(i = 0; i < NJ->nSeq; i++) {
	printf("%s",aln->names[i]);
	for (j = 0; j < NJ->nSeq; j++) {
	  besthit_t hit;
	  SeqDist(NJ->profiles[i]->codes,NJ->profiles[j]->codes,NJ->nPos,NJ->distance_matrix,/*OUT*/&hit);
	  if (logdist)
	    hit.dist = LogCorrect(hit.dist);
	  /* Make sure -0 prints as 0 */
	  printf(" %f", hit.dist <= 0.0 ? 0.0 : hit.dist);
	}
	printf("\n");
      }
    } else {
      /* reset counters*/
      profileOps = 0;
      outprofileOps = 0;
      seqOps = 0;
      profileAvgOps = 0;
      nHillBetter = 0;
      nCloseUsed = 0;
      nClose2Used = 0;
      nRefreshTopHits = 0;
      nVisibleUpdate = 0;
      nNNI = 0;
      nML_NNI = 0;
      nProfileFreqAlloc = 0;
      nProfileFreqAvoid = 0;
      szAllAlloc = 0;
      mymallocUsed = 0;
      maxmallocHeap = 0;
      nLkCompute = 0;
      nPosteriorCompute = 0;
      nAAPosteriorExact = 0;
      nAAPosteriorRough = 0;
      nStarTests = 0;

      uniquify_t *unique = UniquifyAln(aln);
      ProgressReport("Identified unique sequences",0,0,0,0);

      /* read constraints */
      alignment_t *constraints = NULL;
      char **uniqConstraints = NULL;
      if (constraintsFile != NULL) {
	constraints = ReadAlignment(fpConstraints, bQuote);
	if (constraints->nSeq < 4) {
	  fprintf(stderr, "Warning: constraints file with less than 4 sequences ignored:\nalignment #%d in %s\n",
		  iAln+1, constraintsFile);
	  constraints = FreeAlignment(constraints);
	} else {
	  uniqConstraints = AlnToConstraints(constraints, unique, hashnames);
	  ProgressReport("Read the constraints",0,0,0,0);
	}
      }	/* end load constraints */

      transition_matrix_t *transmat = NULL;
      if (nCodes == 20) {
        transmat = transitionFile? ReadAATransitionMatrix(transitionFile) :
          (bUseLg? CreateTransitionMatrix(matrixLG08,statLG08) : 
           (bUseWag? CreateTransitionMatrix(matrixWAG01,statWAG01) :
            CreateTransitionMatrix(matrixJTT92,statJTT92)));
      } else if (nCodes == 4 && bUseGtr && (bUseGtrRates || bUseGtrFreq)) {
	transmat = CreateGTR(gtrrates,gtrfreq);
      }
      NJ_t *NJ = InitNJ(unique->uniqueSeq, unique->nUnique, aln->nPos,
			uniqConstraints,
			uniqConstraints != NULL ? constraints->nPos : 0, /* nConstraints */
			distance_matrix,
			transmat);
      if (verbose>2) fprintf(stderr, "read %s seqs %d (%d unique) positions %d nameLast %s seqLast %s\n",
			     fileName ? fileName : "standard input",
			     aln->nSeq, unique->nUnique, aln->nPos, aln->names[aln->nSeq-1], aln->seqs[aln->nSeq-1]);
      FreeAlignmentSeqs(/*IN/OUT*/aln); /*no longer needed*/
      if (fpInTree != NULL) {
	if (intree1)
	  fseek(fpInTree, 0L, SEEK_SET);
	ReadTree(/*IN/OUT*/NJ, /*IN*/unique, /*IN*/hashnames, /*READ*/fpInTree);
	if (verbose > 2)
	  fprintf(stderr, "Read tree from %s\n", intreeFile);
	if (verbose > 2)
	  PrintNJ(stderr, NJ, aln->names, unique, /*support*/false, bQuote);
      } else {
	FastNJ(NJ);
      }
      LogTree("NJ", 0, fpLog, NJ, aln->names, unique, bQuote);

      /* profile-frequencies for the "up-profiles" in ReliabilityNJ take only diameter(Tree)*L*a
	 space not N*L*a space, because we can free them as we go.
	 And up-profile by their nature tend to be complicated.
	 So save the profile-frequency memory allocation counters now to exclude later results.
      */
#ifdef TRACK_MEMORY
      long svProfileFreqAlloc = nProfileFreqAlloc;
      long svProfileFreqAvoid = nProfileFreqAvoid;
#endif
      int nniToDo = nni == -1 ? (int)(0.5 + 4.0 * log(NJ->nSeq)/log(2)) : nni;
      int sprRemaining = spr;
      int MLnniToDo = (MLnni != -1) ? MLnni : (int)(0.5 + 2.0*log(NJ->nSeq)/log(2));
      if(verbose>0) {
	if (fpInTree == NULL)
	  fprintf(stderr, "Initial topology in %.2f seconds\n", clockDiff(&clock_start));
	if (spr > 0 || nniToDo > 0 || MLnniToDo > 0)
	  fprintf(stderr,"Refining topology: %d rounds ME-NNIs, %d rounds ME-SPRs, %d rounds ML-NNIs\n", nniToDo, spr, MLnniToDo);
      }  

      if (nniToDo>0) {
	int i;
	bool bConverged = false;
	nni_stats_t *nni_stats = InitNNIStats(NJ);
	for (i=0; i < nniToDo; i++) {
	  double maxDelta;
	  if (!bConverged) {
	    int nChange = NNI(/*IN/OUT*/NJ, i, nniToDo, /*use ml*/false, /*IN/OUT*/nni_stats, /*OUT*/&maxDelta);
	    LogTree("ME_NNI%d",i+1, fpLog, NJ, aln->names, unique, bQuote);
	    if (nChange == 0) {
	      bConverged = true;
	      if (verbose>1)
		fprintf(stderr, "Min_evolution NNIs converged at round %d -- skipping some rounds\n", i+1);
	      if (fpLog)
		fprintf(fpLog, "Min_evolution NNIs converged at round %d -- skipping some rounds\n", i+1);
	    }
	  }

	  /* Interleave SPRs with NNIs (typically 1/3rd NNI, SPR, 1/3rd NNI, SPR, 1/3rd NNI */
	  if (sprRemaining > 0 && (nniToDo/(spr+1) > 0 && ((i+1) % (nniToDo/(spr+1))) == 0)) {
	    SPR(/*IN/OUT*/NJ, maxSPRLength, spr-sprRemaining, spr);
	    LogTree("ME_SPR%d",spr-sprRemaining+1, fpLog, NJ, aln->names, unique, bQuote);
	    sprRemaining--;
	    /* Restart the NNIs -- set all ages to 0, etc. */
	    bConverged = false;
	    nni_stats = FreeNNIStats(nni_stats, NJ);
	    nni_stats = InitNNIStats(NJ);
	  }
	}
	nni_stats = FreeNNIStats(nni_stats, NJ);
      }
      while(sprRemaining > 0) {	/* do any remaining SPR rounds */
	SPR(/*IN/OUT*/NJ, maxSPRLength, spr-sprRemaining, spr);
	LogTree("ME_SPR%d",spr-sprRemaining+1, fpLog, NJ, aln->names, unique, bQuote);
	sprRemaining--;
      }

      /* In minimum-evolution mode, update branch lengths, even if no NNIs or SPRs,
	 so that they are log-corrected, do not include penalties from constraints,
	 and avoid errors due to approximation of out-distances.
	 If doing maximum-likelihood NNIs, then we'll also use these
	 to get estimates of starting distances for quartets, etc.
	*/
      UpdateBranchLengths(/*IN/OUT*/NJ);
      LogTree("ME_Lengths",0, fpLog, NJ, aln->names, unique, bQuote);

      double total_len = 0;
      int iNode;
      for (iNode = 0; iNode < NJ->maxnode; iNode++)
	total_len += fabs(NJ->branchlength[iNode]);

      if (verbose>0) {
	fprintf(stderr, "Total branch-length %.3f after %.2f sec\n",
		total_len, clockDiff(&clock_start));
	fflush(stderr);
      }
      if (fpLog) {
	fprintf(fpLog, "Total branch-length %.3f after %.2f sec\n",
		total_len, clockDiff(&clock_start));
	fflush(stderr);
      }

#ifdef TRACK_MEMORY
  if (verbose>1) {
    struct mallinfo mi = mallinfo();
    fprintf(stderr, "Memory @ end of ME phase: %.2f MB (%.1f byte/pos) useful %.2f expected %.2f\n",
	    (mi.arena+mi.hblkhd)/1.0e6, (mi.arena+mi.hblkhd)/(double)(NJ->nSeq*(double)NJ->nPos),
	    mi.uordblks/1.0e6, mymallocUsed/1e6);
  }
#endif

      SplitCount_t splitcount = {0,0,0,0,0.0,0.0};

      if (MLnniToDo > 0 || MLlen) {
	bool warn_len = total_len/NJ->maxnode < 0.001 && MLMinBranchLengthTolerance > 1.0/aln->nPos;
	bool warn = warn_len || (total_len/NJ->maxnode < 0.001 && aln->nPos >= 10000);
	if (warn)
	  fprintf(stderr, "\nWARNING! This alignment consists of closely-related and very-long sequences.\n");
	if (warn_len)
	  fprintf(stderr,
		  "This version of FastTree may not report reasonable branch lengths!\n"
#ifdef USE_DOUBLE
		  "Consider changing MLMinBranchLengthTolerance.\n"
#else
		  "Consider recompiling FastTree with -DUSE_DOUBLE.\n"
#endif
		  "For more information, visit\n"
		  "http://www.microbesonline.org/fasttree/#BranchLen\n\n");
	if (warn)
	  fprintf(stderr, "WARNING! FastTree (or other standard maximum-likelihood tools)\n"
		  "may not be appropriate for aligments of very closely-related sequences\n"
		  "like this one, as FastTree does not account for recombination or gene conversion\n\n");

	/* Do maximum-likelihood computations */
	/* Convert profiles to use the transition matrix */
	distance_matrix_t *tmatAsDist = TransMatToDistanceMat(/*OPTIONAL*/NJ->transmat);
	RecomputeProfiles(NJ, /*OPTIONAL*/tmatAsDist);
	tmatAsDist = myfree(tmatAsDist, sizeof(distance_matrix_t));
	double lastloglk = -1e20;
	nni_stats_t *nni_stats = InitNNIStats(NJ);
	bool resetGtr = nCodes == 4 && bUseGtr && !bUseGtrRates;

	if (MLlen) {
	  int iRound;
	  int maxRound = (int)(0.5 + log(NJ->nSeq)/log(2));
	  double dLastLogLk = -1e20;
	  for (iRound = 1; iRound <= maxRound; iRound++) {
	    int node;
	    numeric_t *oldlength = (numeric_t*)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
	    for (node = 0; node < NJ->maxnode; node++)
	      oldlength[node] = NJ->branchlength[node];
	    OptimizeAllBranchLengths(/*IN/OUT*/NJ);
	    LogTree("ML_Lengths",iRound, fpLog, NJ, aln->names, unique, bQuote);
	    double dMaxChange = 0; /* biggest change in branch length */
	    for (node = 0; node < NJ->maxnode; node++) {
	      double d = fabs(oldlength[node] - NJ->branchlength[node]);
	      if (dMaxChange < d)
		dMaxChange = d;
	    }
	    oldlength = myfree(oldlength, sizeof(numeric_t)*NJ->maxnodes);
	    double loglk = TreeLogLk(NJ, /*site_likelihoods*/NULL);
	    bool bConverged = iRound > 1 && (dMaxChange < 0.001 || loglk < (dLastLogLk+treeLogLkDelta));
	    if (verbose)
	      fprintf(stderr, "%d rounds ML lengths: LogLk %s= %.3lf Max-change %.4lf%s Time %.2f\n",
		      iRound,
		      exactML || nCodes != 20 ? "" : "~",
		      loglk,
		      dMaxChange,
		      bConverged ? " (converged)" : "",
		      clockDiff(&clock_start));
	    if (fpLog)
	      fprintf(fpLog, "TreeLogLk\tLength%d\t%.4lf\tMaxChange\t%.4lf\n",
		      iRound, loglk, dMaxChange);
	    if (iRound == 1) {
	      if (resetGtr)
		SetMLGtr(/*IN/OUT*/NJ, bUseGtrFreq ? gtrfreq : NULL, fpLog);
	      SetMLRates(/*IN/OUT*/NJ, nRateCats);
	      LogMLRates(fpLog, NJ);
	    }
	    if (bConverged)
	      break;
	  }
	}

	if (MLnniToDo > 0) {
	  /* This may help us converge faster, and is fast */
	  OptimizeAllBranchLengths(/*IN/OUT*/NJ);
	  LogTree("ML_Lengths%d",1, fpLog, NJ, aln->names, unique, bQuote);
	}

	int iMLnni;
	double maxDelta;
	bool bConverged = false;
	for (iMLnni = 0; iMLnni < MLnniToDo; iMLnni++) {
	  int changes = NNI(/*IN/OUT*/NJ, iMLnni, MLnniToDo, /*use ml*/true, /*IN/OUT*/nni_stats, /*OUT*/&maxDelta);
	  LogTree("ML_NNI%d",iMLnni+1, fpLog, NJ, aln->names, unique, bQuote);
	  double loglk = TreeLogLk(NJ, /*site_likelihoods*/NULL);
	  bool bConvergedHere = (iMLnni > 0) && ((loglk < lastloglk + treeLogLkDelta) || maxDelta < treeLogLkDelta);
	  if (verbose)
	    fprintf(stderr, "ML-NNI round %d: LogLk %s= %.3f NNIs %d max delta %.2f Time %.2f%s\n",
		    iMLnni+1,
		    exactML || nCodes != 20 ? "" : "~",
		    loglk, changes, maxDelta,  clockDiff(&clock_start),
		    bConverged ? " (final)" : "");
	  if (fpLog)
	    fprintf(fpLog, "TreeLogLk\tML_NNI%d\t%.4lf\tMaxChange\t%.4lf\n", iMLnni+1, loglk, maxDelta);
	  if (bConverged)
	    break;		/* we did our extra round */
	  if (bConvergedHere)
	    bConverged = true;
	  if (bConverged || iMLnni == MLnniToDo-2) {
	    /* last round uses high-accuracy seettings -- reset NNI stats to tone down heuristics */
	    nni_stats = FreeNNIStats(nni_stats, NJ);
	    nni_stats = InitNNIStats(NJ);
	    if (verbose)
	      fprintf(stderr, "Turning off heuristics for final round of ML NNIs%s\n",
		      bConvergedHere? " (converged)" : "");
	    if (fpLog)
	      fprintf(fpLog, "Turning off heuristics for final round of ML NNIs%s\n",
		      bConvergedHere? " (converged)" : "");
	  }
	  lastloglk = loglk;
	  if (iMLnni == 0 && NJ->rates.nRateCategories == 1) {
	    if (resetGtr)
	      SetMLGtr(/*IN/OUT*/NJ, bUseGtrFreq ? gtrfreq : NULL, fpLog);
	    SetMLRates(/*IN/OUT*/NJ, nRateCats);
	    LogMLRates(fpLog, NJ);
	  }
	}
	nni_stats = FreeNNIStats(nni_stats, NJ);	

	/* This does not take long and improves the results */
	if (MLnniToDo > 0) {
	  OptimizeAllBranchLengths(/*IN/OUT*/NJ);
	  LogTree("ML_Lengths%d",2, fpLog, NJ, aln->names, unique, bQuote);
	  if (verbose || fpLog) {
	    double loglk = TreeLogLk(NJ, /*site_likelihoods*/NULL);
	    if (verbose)
	      fprintf(stderr, "Optimize all lengths: LogLk %s= %.3f Time %.2f\n",
		      exactML || nCodes != 20 ? "" : "~",
		      loglk, 
		      clockDiff(&clock_start));
	    if (fpLog) {
	      fprintf(fpLog, "TreeLogLk\tML_Lengths%d\t%.4f\n", 2, loglk);
	      fflush(fpLog);
	    }
	  }
	}

	/* Count bad splits and compute SH-like supports if desired */
	if ((MLnniToDo > 0 && !fastest) || nBootstrap > 0)
	  TestSplitsML(NJ, /*OUT*/&splitcount, nBootstrap);

	/* Compute gamma-based likelihood? */
	if (gammaLogLk && nRateCats > 1) {
	  numeric_t *rates = MLSiteRates(nRateCats);
	  double *site_loglk = MLSiteLikelihoodsByRate(NJ, rates, nRateCats);
	  double scale = RescaleGammaLogLk(NJ->nPos, nRateCats, rates, /*IN*/site_loglk, /*OPTIONAL*/fpLog);
	  rates = myfree(rates, sizeof(numeric_t) * nRateCats);
	  site_loglk = myfree(site_loglk, sizeof(double) * nRateCats * NJ->nPos);

	  for (i = 0; i < NJ->maxnodes; i++)
	    NJ->branchlength[i] *= scale;
	}
      } else {
	/* Minimum evolution supports */
	TestSplitsMinEvo(NJ, /*OUT*/&splitcount);
	if (nBootstrap > 0)
	  ReliabilityNJ(NJ, nBootstrap);
      }

      for (i = 0; i < nFPs; i++) {
	FILE *fp = fps[i];
	fprintf(fp, "Total time: %.2f seconds Unique: %d/%d Bad splits: %d/%d",
		clockDiff(&clock_start),
		NJ->nSeq, aln->nSeq,
		splitcount.nBadSplits, splitcount.nSplits);
	if (splitcount.dWorstDeltaUnconstrained >  0)
	  fprintf(fp, " Worst %sdelta-%s %.3f",
		  uniqConstraints != NULL ? "unconstrained " : "",
		  (MLnniToDo > 0 || MLlen) ? "LogLk" : "Len",
		  splitcount.dWorstDeltaUnconstrained);
	fprintf(fp,"\n");
	if (NJ->nSeq > 3 && NJ->nConstraints > 0) {
	    fprintf(fp, "Violating constraints: %d both bad: %d",
		    splitcount.nConstraintViolations, splitcount.nBadBoth);
	    if (splitcount.dWorstDeltaConstrained >  0)
	      fprintf(fp, " Worst delta-%s due to constraints: %.3f",
		      (MLnniToDo > 0 || MLlen) ? "LogLk" : "Len",
		      splitcount.dWorstDeltaConstrained);
	    fprintf(fp,"\n");
	}
	if (verbose > 1 || fp == fpLog) {
	  double dN2 = NJ->nSeq*(double)NJ->nSeq;
	  fprintf(fp, "Dist/N**2: by-profile %.3f (out %.3f) by-leaf %.3f avg-prof %.3f\n",
		  profileOps/dN2, outprofileOps/dN2, seqOps/dN2, profileAvgOps/dN2);
	  if (nCloseUsed>0 || nClose2Used > 0 || nRefreshTopHits>0)
	    fprintf(fp, "Top hits: close neighbors %ld/%d 2nd-level %ld refreshes %ld",
		    nCloseUsed, NJ->nSeq, nClose2Used, nRefreshTopHits);
	  if(!slow) fprintf(fp, " Hill-climb: %ld Update-best: %ld\n", nHillBetter, nVisibleUpdate);
	  if (nniToDo > 0 || spr > 0 || MLnniToDo > 0)
	    fprintf(fp, "NNI: %ld SPR: %ld ML-NNI: %ld\n", nNNI, nSPR, nML_NNI);
	  if (MLnniToDo > 0) {
	    fprintf(fp, "Max-lk operations: lk %ld posterior %ld", nLkCompute, nPosteriorCompute);
	    if (nAAPosteriorExact > 0 || nAAPosteriorRough > 0)
	      fprintf(fp, " approximate-posteriors %.2f%%",
		      (100.0*nAAPosteriorRough)/(double)(nAAPosteriorExact+nAAPosteriorRough));
	    if (mlAccuracy < 2)
	      fprintf(fp, " star-only %ld", nStarTests);
	    fprintf(fp, "\n");
	  }
	}
#ifdef TRACK_MEMORY
	fprintf(fp, "Memory: %.2f MB (%.1f byte/pos) ",
		maxmallocHeap/1.0e6, maxmallocHeap/(double)(aln->nSeq*(double)aln->nPos));
	/* Only report numbers from before we do reliability estimates */
	fprintf(fp, "profile-freq-alloc %ld avoided %.2f%%\n", 
		svProfileFreqAlloc,
		svProfileFreqAvoid > 0 ?
		100.0*svProfileFreqAvoid/(double)(svProfileFreqAlloc+svProfileFreqAvoid)
		: 0);
#endif
	fflush(fp);
      }
      PrintNJ(fpOut, NJ, aln->names, unique, /*support*/nBootstrap > 0, bQuote);
      fflush(fpOut);
      if (fpLog) {
	fprintf(fpLog,"TreeCompleted\n");
	fflush(fpLog);
      }
      FreeNJ(NJ);
      if (uniqConstraints != NULL)
	uniqConstraints = myfree(uniqConstraints, sizeof(char*) * unique->nUnique);
      constraints = FreeAlignment(constraints);
      unique = FreeUniquify(unique);
    } /* end build tree */
    hashnames = FreeHashtable(hashnames);
    aln = FreeAlignment(aln);
  } /* end loop over alignments */
  if (fpLog != NULL)
    fclose(fpLog);
  if (fpOut != stdout) fclose(fpOut);
  exit(0);
}

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

NJ_t *InitNJ(char **sequences, int nSeq, int nPos,
	     /*OPTIONAL*/char **constraintSeqs, int nConstraints,
	     /*OPTIONAL*/distance_matrix_t *distance_matrix,
	     /*OPTIONAL*/transition_matrix_t *transmat) {
  int iNode;

  NJ_t *NJ = (NJ_t*)mymalloc(sizeof(NJ_t));
  NJ->root = -1; 		/* set at end of FastNJ() */
  NJ->maxnode = NJ->nSeq = nSeq;
  NJ->nPos = nPos;
  NJ->maxnodes = 2*nSeq;
  NJ->seqs = sequences;
  NJ->distance_matrix = distance_matrix;
  NJ->transmat = transmat;
  NJ->nConstraints = nConstraints;
  NJ->constraintSeqs = constraintSeqs;

  NJ->profiles = (profile_t **)mymalloc(sizeof(profile_t*) * NJ->maxnodes);

  unsigned long counts[256];
  int i;
  for (i = 0; i < 256; i++)
    counts[i] = 0;
  for (iNode = 0; iNode < NJ->nSeq; iNode++) {
    NJ->profiles[iNode] = SeqToProfile(NJ, NJ->seqs[iNode], nPos,
				       constraintSeqs != NULL ? constraintSeqs[iNode] : NULL,
				       nConstraints,
				       iNode,
				       /*IN/OUT*/counts);
  }
  unsigned long totCount = 0;
  for (i = 0; i < 256; i++)
    totCount += counts[i];

  /* warnings about unknown characters */
  for (i = 0; i < 256; i++) {
    if (counts[i] == 0 || i == '.' || i == '-')
      continue;
    unsigned char *codesP;
    bool bMatched = false;
    for (codesP = codesString; *codesP != '\0'; codesP++) {
      if (*codesP == i || tolower(*codesP) == i) {
	bMatched = true;
	break;
      }
    }
    if (!bMatched)
      fprintf(stderr, "Ignored unknown character %c (seen %lu times)\n", i, counts[i]);
  }
    

  /* warnings about the counts */
  double fACGTUN = (counts['A'] + counts['C'] + counts['G'] + counts['T'] + counts['U'] + counts['N']
		    + counts['a'] + counts['c'] + counts['g'] + counts['t'] + counts['u'] + counts['n'])
    / (double)(totCount - counts['-'] - counts['.']);
  if (nCodes == 4 && fACGTUN < 0.9)
    fprintf(stderr, "WARNING! ONLY %.1f%% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A NUCLEOTIDE ALIGNMENT?\n",
	    100.0 * fACGTUN);
  else if (nCodes == 20 && fACGTUN >= 0.9)
    fprintf(stderr, "WARNING! %.1f%% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN ALIGNMENT?\n",
	    100.0 * fACGTUN);

  if(verbose>10) fprintf(stderr,"Made sequence profiles\n");
  for (iNode = NJ->nSeq; iNode < NJ->maxnodes; iNode++) 
    NJ->profiles[iNode] = NULL; /* not yet exists */

  NJ->outprofile = OutProfile(NJ->profiles, NJ->nSeq,
			      NJ->nPos, NJ->nConstraints,
			      NJ->distance_matrix);
  if(verbose>10) fprintf(stderr,"Made out-profile\n");

  NJ->totdiam = 0.0;

  NJ->diameter = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->diameter[iNode] = 0;

  NJ->varDiameter = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->varDiameter[iNode] = 0;

  NJ->selfdist = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->selfdist[iNode] = 0;

  NJ->selfweight = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->nSeq; iNode++)
    NJ->selfweight[iNode] = NJ->nPos - NGaps(NJ,iNode);

  NJ->outDistances = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
  NJ->nOutDistActive = (int *)mymalloc(sizeof(int)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++)
    NJ->nOutDistActive[iNode] = NJ->nSeq * 10; /* unreasonably high value */
  NJ->parent = NULL;		/* so SetOutDistance ignores it */
  for (iNode = 0; iNode < NJ->nSeq; iNode++)
    SetOutDistance(/*IN/UPDATE*/NJ, iNode, /*nActive*/NJ->nSeq);

  if (verbose>2) {
    for (iNode = 0; iNode < 4 && iNode < NJ->nSeq; iNode++)
      fprintf(stderr, "Node %d outdist %f\n", iNode, NJ->outDistances[iNode]);
  }

  NJ->parent = (int *)mymalloc(sizeof(int)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->parent[iNode] = -1;

  NJ->branchlength = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes); /* distance to parent */
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->branchlength[iNode] = 0;

  NJ->support = (numeric_t *)mymalloc(sizeof(numeric_t)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->support[iNode] = -1.0;

  NJ->child = (children_t*)mymalloc(sizeof(children_t)*NJ->maxnodes);
  for (iNode= 0; iNode < NJ->maxnode; iNode++) NJ->child[iNode].nChild = 0;

  NJ->rates.nRateCategories = 0;
  NJ->rates.rates = NULL;
  NJ->rates.ratecat = NULL;
  AllocRateCategories(&NJ->rates, 1, NJ->nPos);
  return(NJ);
}

NJ_t *FreeNJ(NJ_t *NJ) {
  if (NJ==NULL)
    return(NJ);

  int i;
  for (i=0; i < NJ->maxnode; i++)
    NJ->profiles[i] = FreeProfile(NJ->profiles[i], NJ->nPos, NJ->nConstraints);
  NJ->profiles = myfree(NJ->profiles, sizeof(profile_t*) * NJ->maxnodes);
  NJ->outprofile = FreeProfile(NJ->outprofile, NJ->nPos, NJ->nConstraints);
  NJ->diameter = myfree(NJ->diameter, sizeof(numeric_t)*NJ->maxnodes);
  NJ->varDiameter = myfree(NJ->varDiameter, sizeof(numeric_t)*NJ->maxnodes);
  NJ->selfdist = myfree(NJ->selfdist, sizeof(numeric_t)*NJ->maxnodes);
  NJ->selfweight = myfree(NJ->selfweight, sizeof(numeric_t)*NJ->maxnodes);
  NJ->outDistances = myfree(NJ->outDistances, sizeof(numeric_t)*NJ->maxnodes);
  NJ->nOutDistActive = myfree(NJ->nOutDistActive, sizeof(int)*NJ->maxnodes);
  NJ->parent = myfree(NJ->parent, sizeof(int)*NJ->maxnodes);
  NJ->branchlength = myfree(NJ->branchlength, sizeof(numeric_t)*NJ->maxnodes);
  NJ->support = myfree(NJ->support, sizeof(numeric_t)*NJ->maxnodes);
  NJ->child = myfree(NJ->child, sizeof(children_t)*NJ->maxnodes);
  NJ->transmat = myfree(NJ->transmat, sizeof(transition_matrix_t));
  AllocRateCategories(&NJ->rates, 0, NJ->nPos);
  return(myfree(NJ, sizeof(NJ_t)));
}

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

void ExhaustiveNJSearch(NJ_t *NJ, int nActive, /*OUT*/besthit_t *join) {
  join->i = -1;
  join->j = -1;
  join->weight = 0;
  join->dist = 1e20;
  join->criterion = 1e20;
  double bestCriterion = 1e20;

  int i, j;
  for (i = 0; i < NJ->maxnode-1; i++) {
    if (NJ->parent[i] < 0) {
      for (j = i+1; j < NJ->maxnode; j++) {
	if (NJ->parent[j] < 0) {
	  besthit_t hit;
	  hit.i = i;
	  hit.j = j;
	  SetDistCriterion(NJ, nActive, /*IN/OUT*/&hit);
	  if (hit.criterion < bestCriterion) {
	    *join = hit;
	    bestCriterion = hit.criterion;
	  }
	}
      }
    }
  }
  assert (join->i >= 0 && join->j >= 0);
}

void FastNJSearch(NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *besthits, /*OUT*/besthit_t *join) {
  join->i = -1;
  join->j = -1;
  join->dist = 1e20;
  join->weight = 0;
  join->criterion = 1e20;
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    int jNode = besthits[iNode].j;
    if (NJ->parent[iNode] < 0 && NJ->parent[jNode] < 0) { /* both i and j still active */
      /* recompute criterion to reflect the current out-distances */
      SetCriterion(NJ, nActive, /*IN/OUT*/&besthits[iNode]);
      if (besthits[iNode].criterion < join->criterion)
	*join = besthits[iNode];      
    }
  }

  if(!fastest) {
    int changed;
    do {
      changed = 0;
      assert(join->i >= 0 && join->j >= 0);
      SetBestHit(join->i, NJ, nActive, /*OUT*/&besthits[join->i], /*OUT IGNORED*/NULL);
      if (besthits[join->i].j != join->j) {
	changed = 1;
	if (verbose>2)
	  fprintf(stderr,"BetterI\t%d\t%d\t%d\t%d\t%f\t%f\n",
		  join->i,join->j,besthits[join->i].i,besthits[join->i].j,
		  join->criterion,besthits[join->i].criterion);
      }
      
      /* Save the best hit either way, because the out-distance has probably changed
	 since we started the computation. */
      join->j = besthits[join->i].j;
      join->weight = besthits[join->i].weight;
      join->dist = besthits[join->i].dist;
      join->criterion = besthits[join->i].criterion;
      
      SetBestHit(join->j, NJ, nActive, /*OUT*/&besthits[join->j], /*OUT IGNORE*/NULL);
      if (besthits[join->j].j != join->i) {
	changed = 1;
	if (verbose>2)
	  fprintf(stderr,"BetterJ\t%d\t%d\t%d\t%d\t%f\t%f\n",
		  join->i,join->j,besthits[join->j].i,besthits[join->j].j,
		  join->criterion,besthits[join->j].criterion);
	join->i = besthits[join->j].j;
	join->weight = besthits[join->j].weight;
	join->dist = besthits[join->j].dist;
	join->criterion = besthits[join->j].criterion;
      }
      if(changed) nHillBetter++;
    } while(changed);
  }
}

/* A token is one of ():;, or an alphanumeric string without whitespace
   Any whitespace between tokens is ignored */
char *ReadTreeToken(FILE *fp) {
  static char buf[BUFFER_SIZE];
  int len = 0;
  int c;
  for (c = fgetc(fp); c != EOF; c = fgetc(fp)) {
    if (c == '(' || c == ')' || c == ':' || c == ';' || c == ',') {
      /* standalone token */
      if (len == 0) {
	buf[len++] = c;
	buf[len] = '\0';
	return(buf);
      } else {
	ungetc(c, fp);
	buf[len] = '\0';
	return(buf);
      }
    } else if (isspace(c)) {
      if (len > 0) {
	buf[len] = '\0';
	return(buf);
      }
      /* else ignore whitespace at beginning of token */
    } else {
      /* not whitespace or standalone token */
      buf[len++] = c;
      if (len >= BUFFER_SIZE) {
	buf[BUFFER_SIZE-1] = '\0';
	fprintf(stderr, "Token too long in tree file, token begins with\n%s\n", buf);
	exit(1);
      }
    }
  }
  if (len > 0) {
    /* return the token we have so far */
    buf[len] = '\0';
    return(buf);
  }
  /* else */
  return(NULL);
}

void ReadTreeError(char *err, char *token) {
  fprintf(stderr, "Tree parse error: unexpected token '%s' -- %s\n",
	  token == NULL ? "(End of file)" : token,
	  err);
  exit(1);
}

void ReadTreeAddChild(int parent, int child, /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children) {
  assert(parent >= 0);
  assert(child >= 0);
  assert(parents[child] < 0);
  assert(children[parent].nChild < 3);
  parents[child] = parent;
  children[parent].child[children[parent].nChild++] = child;
}

void ReadTreeMaybeAddLeaf(int parent, char *name,
			  hashstrings_t *hashnames, uniquify_t *unique,
			  /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children) {
  hashiterator_t hi = FindMatch(hashnames,name);
  if (HashCount(hashnames,hi) != 1)
    ReadTreeError("not recognized as a sequence name", name);

  int iSeqNonunique = HashFirst(hashnames,hi);
  assert(iSeqNonunique >= 0 && iSeqNonunique < unique->nSeq);
  int iSeqUnique = unique->alnToUniq[iSeqNonunique];
  assert(iSeqUnique >= 0 && iSeqUnique < unique->nUnique);
  /* Either record this leaves' parent (if it is -1) or ignore this leaf (if already seen) */
  if (parents[iSeqUnique] < 0) {
    ReadTreeAddChild(parent, iSeqUnique, /*IN/OUT*/parents, /*IN/OUT*/children);
    if(verbose > 5)
      fprintf(stderr, "Found leaf uniq%d name %s child of %d\n", iSeqUnique, name, parent);
  } else {
    if (verbose > 5)
      fprintf(stderr, "Skipped redundant leaf uniq%d name %s\n", iSeqUnique, name);
  }
}

void ReadTreeRemove(/*IN/OUT*/int *parents, /*IN/OUT*/children_t *children, int node) {
  if(verbose > 5)
    fprintf(stderr,"Removing node %d parent %d\n", node, parents[node]);
  assert(parents[node] >= 0);
  int parent = parents[node];
  parents[node] = -1;
  children_t *pc = &children[parent];
  int oldn;
  for (oldn = 0; oldn < pc->nChild; oldn++) {
    if (pc->child[oldn] == node)
      break;
  }
  assert(oldn < pc->nChild);

  /* move successor nodes back in child list and shorten list */
  int i;
  for (i = oldn; i < pc->nChild-1; i++)
    pc->child[i] = pc->child[i+1];
  pc->nChild--;

  /* add its children to parent's child list */
  children_t *nc = &children[node];
  if (nc->nChild > 0) {
    assert(nc->nChild<=2);
    assert(pc->nChild < 3);
    assert(pc->nChild + nc->nChild <= 3);
    int j;
    for (j = 0; j < nc->nChild; j++) {
      if(verbose > 5)
	fprintf(stderr,"Repointing parent %d to child %d\n", parent, nc->child[j]);
      pc->child[pc->nChild++] = nc->child[j];
      parents[nc->child[j]] = parent;
    }
    nc->nChild = 0;
  }
}  

void ReadTree(/*IN/OUT*/NJ_t *NJ,
	      /*IN*/uniquify_t *unique,
	      /*IN*/hashstrings_t *hashnames,
	      /*READ*/FILE *fpInTree) {
  assert(NJ->nSeq == unique->nUnique);
  /* First, do a preliminary parse of the tree to with non-unique leaves ignored
     We need to store this separately from NJ because it may have too many internal nodes
     (matching sequences show up once in the NJ but could be in multiple places in the tree)
     Will use iUnique as the index of nodes, as in the NJ structure
  */
  int maxnodes = unique->nSeq*2;
  int maxnode = unique->nSeq;
  int *parent = (int*)mymalloc(sizeof(int)*maxnodes);
  children_t *children = (children_t *)mymalloc(sizeof(children_t)*maxnodes);
  int root = maxnode++;
  int i;
  for (i = 0; i < maxnodes; i++) {
    parent[i] = -1;
    children[i].nChild = 0;
  }

  /* The stack is the current path to the root, with the root at the first (top) position */
  int stack_size = 1;
  int *stack = (int*)mymalloc(sizeof(int)*maxnodes);
  stack[0] = root;
  int nDown = 0;
  int nUp = 0;

  char *token;
  token = ReadTreeToken(fpInTree);
  if (token == NULL || *token != '(')
    ReadTreeError("No '(' at start", token);
  /* nDown is still 0 because we have created the root */

  while ((token = ReadTreeToken(fpInTree)) != NULL) {
    if (nDown > 0) {		/* In a stream of parentheses */
      if (*token == '(')
	nDown++;
      else if (*token == ',' || *token == ';' || *token == ':' || *token == ')')
	ReadTreeError("while reading parentheses", token);
      else {
	/* Add intermediate nodes if nDown was > 1 (for nDown=1, the only new node is the leaf) */
	while (nDown-- > 0) {
	  int new = maxnode++;
	  assert(new < maxnodes);
	  ReadTreeAddChild(stack[stack_size-1], new, /*IN/OUT*/parent, /*IN/OUT*/children);
	  if(verbose > 5)
	    fprintf(stderr, "Added internal child %d of %d, stack size increase to %d\n",
		    new, stack[stack_size-1],stack_size+1);
	  stack[stack_size++] = new;
	  assert(stack_size < maxnodes);
	}
	ReadTreeMaybeAddLeaf(stack[stack_size-1], token,
			     hashnames, unique,
			     /*IN/OUT*/parent, /*IN/OUT*/children);
      }
    } else if (nUp > 0) {
      if (*token == ';') {	/* end the tree? */
	if (nUp != stack_size)
	  ReadTreeError("unbalanced parentheses", token);
	else
	  break;
      } else if (*token == ')')
	nUp++;
      else if (*token == '(')
	ReadTreeError("unexpected '(' after ')'", token);
      else if (*token == ':') {
	token = ReadTreeToken(fpInTree);
	/* Read the branch length and ignore it */
	if (token == NULL || (*token != '-' && !isdigit(*token)))
	  ReadTreeError("not recognized as a branch length", token);
      } else if (*token == ',') {
	/* Go back up the stack the correct #times */
	while (nUp-- > 0) {
	  stack_size--;
	  if(verbose > 5)
	    fprintf(stderr, "Up to nUp=%d stack size %d at %d\n",
		    nUp, stack_size, stack[stack_size-1]);
	  if (stack_size <= 0)
	    ReadTreeError("too many ')'", token);
	}
	nUp = 0;
      } else if (*token == '-' || isdigit(*token))
	; 			/* ignore bootstrap value */
      else
	fprintf(stderr, "Warning while parsing tree: non-numeric label %s for internal node\n",
		token);
    } else if (*token == '(') {
      nDown = 1;
    } else if (*token == ')') {
      nUp = 1;
    } else if (*token == ':') {
      token = ReadTreeToken(fpInTree);
      if (token == NULL || (*token != '-' && !isdigit(*token)))
	ReadTreeError("not recognized as a branch length", token);
    } else if (*token == ',') {
      ;				/* do nothing */
    } else if (*token == ';')
      ReadTreeError("unexpected token", token);
    else
      ReadTreeMaybeAddLeaf(stack[stack_size-1], token,
			   hashnames, unique,
			   /*IN/OUT*/parent, /*IN/OUT*/children);
  }

  /* Verify that all sequences were seen */
  for (i = 0; i < unique->nUnique; i++) {
    if (parent[i] < 0) {
      fprintf(stderr, "Alignment sequence %d (unique %d) absent from input tree\n"
	      "The starting tree (the argument to -intree) must include all sequences in the alignment!\n",
	      unique->uniqueFirst[i], i);
      exit(1);
    }
  }

  /* Simplify the tree -- remove all internal nodes with < 2 children
     Keep trying until no nodes get removed
  */
  int nRemoved;
  do {
    nRemoved = 0;
    /* Here stack is the list of nodes we haven't visited yet while doing
       a tree traversal */
    stack_size = 1;
    stack[0] = root;
    while (stack_size > 0) {
      int node = stack[--stack_size];
      if (node >= unique->nUnique) { /* internal node */
	if (children[node].nChild <= 1) {
	  if (node != root) {
	    ReadTreeRemove(/*IN/OUT*/parent,/*IN/OUT*/children,node);
	    nRemoved++;
	  } else if (node == root && children[node].nChild == 1) {
	    int newroot = children[node].child[0];
	    parent[newroot] = -1;
	    children[root].nChild = 0;
	    nRemoved++;
	    if(verbose > 5)
	      fprintf(stderr,"Changed root from %d to %d\n",root,newroot);
	    root = newroot;
	    stack[stack_size++] = newroot;
	  }
	} else {
	  int j;
	  for (j = 0; j < children[node].nChild; j++) {
	    assert(stack_size < maxnodes);
	    stack[stack_size++] = children[node].child[j];
	    if(verbose > 5)
	      fprintf(stderr,"Added %d to stack\n", stack[stack_size-1]);
	  }
	}
      }
    }
  } while (nRemoved > 0);

  /* Simplify the root node to 3 children if it has 2 */
  if (children[root].nChild == 2) {
    for (i = 0; i < 2; i++) {
      int child = children[root].child[i];
      assert(child >= 0 && child < maxnodes);
      if (children[child].nChild == 2) {
	ReadTreeRemove(parent,children,child); /* replace root -> child -> A,B with root->A,B */
	break;
      }
    }
  }

  for (i = 0; i < maxnodes; i++)
    if(verbose > 5)
      fprintf(stderr,"Simplfied node %d has parent %d nchild %d\n",
	      i, parent[i], children[i].nChild);

  /* Map the remaining internal nodes to NJ nodes */
  int *map = (int*)mymalloc(sizeof(int)*maxnodes);
  for (i = 0; i < unique->nUnique; i++)
    map[i] = i;
  for (i = unique->nUnique; i < maxnodes; i++)
    map[i] = -1;
  stack_size = 1;
  stack[0] = root;
  while (stack_size > 0) {
    int node = stack[--stack_size];
    if (node >= unique->nUnique) { /* internal node */
      assert(node == root || children[node].nChild > 1);
      map[node] =  NJ->maxnode++;
      for (i = 0; i < children[node].nChild; i++) {
	assert(stack_size < maxnodes);
	stack[stack_size++] = children[node].child[i];
      }
    }
  }
  for (i = 0; i < maxnodes; i++)
    if(verbose > 5)
      fprintf(stderr,"Map %d to %d (parent %d nchild %d)\n",
	      i, map[i], parent[i], children[i].nChild);

  /* Set NJ->parent, NJ->children, NJ->root */
  NJ->root = map[root];
  int node;
  for (node = 0; node < maxnodes; node++) {
    int njnode = map[node];
    if (njnode >= 0) {
      NJ->child[njnode].nChild = children[node].nChild;
      for (i = 0; i < children[node].nChild; i++) {
	assert(children[node].child[i] >= 0 && children[node].child[i] < maxnodes);
	NJ->child[njnode].child[i] = map[children[node].child[i]];
      }
      if (parent[node] >= 0)
	NJ->parent[njnode] = map[parent[node]];
    }
  }

  /* Make sure that parent/child relationships match */
  for (i = 0; i < NJ->maxnode; i++) {
    children_t *c = &NJ->child[i];
    int j;
    for (j = 0; j < c->nChild;j++)
      assert(c->child[j] >= 0 && c->child[j] < NJ->maxnode && NJ->parent[c->child[j]] == i);
  }
  assert(NJ->parent[NJ->root] < 0);

  map = myfree(map,sizeof(int)*maxnodes);
  stack = myfree(stack,sizeof(int)*maxnodes);
  children = myfree(children,sizeof(children_t)*maxnodes);
  parent = myfree(parent,sizeof(int)*maxnodes);

  /* Compute profiles as balanced -- the NNI stage will recompute these
     profiles anyway
  */
  traversal_t traversal = InitTraversal(NJ);
  node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (node >= NJ->nSeq && node != NJ->root)
      SetProfile(/*IN/OUT*/NJ, node, /*noweight*/-1.0);
  }
  traversal = FreeTraversal(traversal,NJ);
}

/* Print topology using node indices as node names */
void PrintNJInternal(FILE *fp, NJ_t *NJ, bool useLen) {
  if (NJ->nSeq < 4) {
    return;
  }
  typedef struct { int node; int end; } stack_t;
  stack_t *stack = (stack_t *)mymalloc(sizeof(stack_t)*NJ->maxnodes);
  int stackSize = 1;
  stack[0].node = NJ->root;
  stack[0].end = 0;

  while(stackSize>0) {
    stack_t *last = &stack[stackSize-1];
    stackSize--;
    /* Save last, as we are about to overwrite it */
    int node = last->node;
    int end = last->end;

    if (node < NJ->nSeq) {
      if (NJ->child[NJ->parent[node]].child[0] != node) fputs(",",fp);
      fprintf(fp, "%d", node);
      if (useLen)
	fprintf(fp, ":%.4f", NJ->branchlength[node]);
    } else if (end) {
      fprintf(fp, ")%d", node);
      if (useLen)
	fprintf(fp, ":%.4f", NJ->branchlength[node]);
    } else {
            if (node != NJ->root && NJ->child[NJ->parent[node]].child[0] != node) fprintf(fp, ",");
      fprintf(fp, "(");
      stackSize++;
      stack[stackSize-1].node = node;
      stack[stackSize-1].end = 1;
      children_t *c = &NJ->child[node];
      /* put children on in reverse order because we use the last one first */
      int i;
      for (i = c->nChild-1; i >=0; i--) {
	stackSize++;
	stack[stackSize-1].node = c->child[i];
	stack[stackSize-1].end = 0;
      }
    }
  }
  fprintf(fp, ";\n");
  stack = myfree(stack, sizeof(stack_t)*NJ->maxnodes);
}

void PrintNJ(FILE *fp, NJ_t *NJ, char **names, uniquify_t *unique, bool bShowSupport, bool bQuote) {
  /* And print the tree: depth first search
   * The stack contains
   * list of remaining children with their depth
   * parent node, with a flag of -1 so I know to print right-paren
   */
  if (NJ->nSeq==1 && unique->alnNext[unique->uniqueFirst[0]] >= 0) {
    /* Special case -- otherwise we end up with double parens */
    int first = unique->uniqueFirst[0];
    assert(first >= 0 && first < unique->nSeq);
    fprintf(fp, bQuote ? "('%s':0.0" : "(%s:0.0", names[first]);
    int iName = unique->alnNext[first];
    while (iName >= 0) {
      assert(iName < unique->nSeq);
      fprintf(fp, bQuote ? ",'%s':0.0" : ",%s:0.0", names[iName]);
      iName = unique->alnNext[iName];
    }
    fprintf(fp,");\n");
    return;
  }

  typedef struct { int node; int end; } stack_t;
  stack_t *stack = (stack_t *)mymalloc(sizeof(stack_t)*NJ->maxnodes);
  int stackSize = 1;
  stack[0].node = NJ->root;
  stack[0].end = 0;

  while(stackSize>0) {
    stack_t *last = &stack[stackSize-1];
    stackSize--;
    /* Save last, as we are about to overwrite it */
    int node = last->node;
    int end = last->end;

    if (node < NJ->nSeq) {
      if (NJ->child[NJ->parent[node]].child[0] != node) fputs(",",fp);
      int first = unique->uniqueFirst[node];
      assert(first >= 0 && first < unique->nSeq);
      /* Print the name, or the subtree of duplicate names */
      if (unique->alnNext[first] == -1) {
	fprintf(fp, bQuote ? "'%s'" : "%s", names[first]);
      } else {
	fprintf(fp, bQuote ? "('%s':0.0" : "(%s:0.0", names[first]);
	int iName = unique->alnNext[first];
	while (iName >= 0) {
	  assert(iName < unique->nSeq);
	  fprintf(fp, bQuote ? ",'%s':0.0" : ",%s:0.0", names[iName]);
	  iName = unique->alnNext[iName];
	}
	fprintf(fp,")");
      }
      /* Print the branch length */
#ifdef USE_DOUBLE
#define FP_FORMAT "%.9f"
#else
#define FP_FORMAT "%.5f"
#endif
      fprintf(fp, ":" FP_FORMAT, NJ->branchlength[node]);
    } else if (end) {
      if (node == NJ->root)
	fprintf(fp, ")");
      else if (bShowSupport)
	fprintf(fp, ")%.3f:" FP_FORMAT, NJ->support[node], NJ->branchlength[node]);
      else
	fprintf(fp, "):" FP_FORMAT, NJ->branchlength[node]);
    } else {
      if (node != NJ->root && NJ->child[NJ->parent[node]].child[0] != node) fprintf(fp, ",");
      fprintf(fp, "(");
      stackSize++;
      stack[stackSize-1].node = node;
      stack[stackSize-1].end = 1;
      children_t *c = &NJ->child[node];
      /* put children on in reverse order because we use the last one first */
      int i;
      for (i = c->nChild-1; i >=0; i--) {
	stackSize++;
	stack[stackSize-1].node = c->child[i];
	stack[stackSize-1].end = 0;
      }
    }
  }
  fprintf(fp, ";\n");
  stack = myfree(stack, sizeof(stack_t)*NJ->maxnodes);
}

alignment_t *ReadAlignment(/*IN*/FILE *fp, bool bQuote) {
  /* bQuote supports the -quote option */
  int nSeq = 0;
  int nPos = 0;
  char **names = NULL;
  char **seqs = NULL;
  char buf[BUFFER_SIZE] = "";
  if (fgets(buf,sizeof(buf),fp) == NULL) {
    fprintf(stderr, "Error reading header line\n");
    exit(1);
  }
  int nSaved = 100;
  if (buf[0] == '>') {
    /* FASTA, truncate names at any of these */
    char *nameStop = bQuote ? "'\t\r\n" : "(),: \t\r\n";
    char *seqSkip = " \t\r\n";	/* skip these characters in the sequence */
    seqs = (char**)mymalloc(sizeof(char*) * nSaved);
    names = (char**)mymalloc(sizeof(char*) * nSaved);

    do {
      /* loop over lines */
      if (buf[0] == '>') {
	/* truncate the name */
	char *p, *q;
	for (p = buf+1; *p != '\0'; p++) {
	  for (q = nameStop; *q != '\0'; q++) {
	    if (*p == *q) {
	      *p = '\0';
	      break;
	    }
	  }
	  if (*p == '\0') break;
	}

	/* allocate space for another sequence */
	nSeq++;
	if (nSeq > nSaved) {
	  int nNewSaved = nSaved*2;
	  seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved, /*copy*/false);
	  names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved, /*copy*/false);
	  nSaved = nNewSaved;
	}
	names[nSeq-1] = (char*)mymemdup(buf+1,strlen(buf));
	seqs[nSeq-1] = NULL;
      } else {
	/* count non-space characters and append to sequence */
	int nKeep = 0;
	char *p, *q;
	for (p=buf; *p != '\0'; p++) {
	  for (q=seqSkip; *q != '\0'; q++) {
	    if (*p == *q)
	      break;
	  }
	  if (*p != *q)
	    nKeep++;
	}
	int nOld = (seqs[nSeq-1] == NULL) ? 0 : strlen(seqs[nSeq-1]);
	seqs[nSeq-1] = (char*)myrealloc(seqs[nSeq-1], nOld, nOld+nKeep+1, /*copy*/false);
	if (nOld+nKeep > nPos)
	  nPos = nOld + nKeep;
	char *out = seqs[nSeq-1] + nOld;
	for (p=buf; *p != '\0'; p++) {
	  for (q=seqSkip; *q != '\0'; q++) {
	    if (*p == *q)
	      break;
	  }
	  if (*p != *q) {
	    *out = *p;
	    out++;
	  }
	}
	assert(out-seqs[nSeq-1] == nKeep + nOld);
	*out = '\0';
      }
    } while(fgets(buf,sizeof(buf),fp) != NULL);

    if (seqs[nSeq-1] == NULL) {
      fprintf(stderr, "No sequence data for last entry %s\n",names[nSeq-1]);
      exit(1);
    }
    names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nSeq, /*copy*/false);
    seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nSeq, /*copy*/false);
  } else {
    /* PHYLIP interleaved-like format
       Allow arbitrary length names, require spaces between names and sequences
       Allow multiple alignments, either separated by a single empty line (e.g. seqboot output)
       or not.
     */
    if (buf[0] == '\n' || buf[0] == '\r') {
      if (fgets(buf,sizeof(buf),fp) == NULL) {
	fprintf(stderr, "Empty header line followed by EOF\n");
	exit(1);
      }
    }
    if (sscanf(buf, "%d%d", &nSeq, &nPos) != 2
      || nSeq < 1 || nPos < 1) {
      fprintf(stderr, "Error parsing header line:%s\n", buf);
      exit(1);
    }
    names = (char **)mymalloc(sizeof(char*) * nSeq);
    seqs = (char **)mymalloc(sizeof(char*) * nSeq);
    nSaved = nSeq;

    int i;
    for (i = 0; i < nSeq; i++) {
      names[i] = NULL;
      seqs[i] = (char *)mymalloc(nPos+1);	/* null-terminate */
      seqs[i][0] = '\0';
    }
    int iSeq = 0;
    
    while(fgets(buf,sizeof(buf),fp)) {
      if ((buf[0] == '\n' || buf[0] == '\r') && (iSeq == nSeq || iSeq == 0)) {
	iSeq = 0;
      } else {
	int j = 0; /* character just past end of name */
	if (buf[0] == ' ') {
	  if (names[iSeq] == NULL) {
	    fprintf(stderr, "No name in phylip line %s", buf);
	    exit(1);
	  }
	} else {
	  while (buf[j] != '\n' && buf[j] != '\0' && buf[j] != ' ')
	    j++;
	  if (buf[j] != ' ' || j == 0) {
	    fprintf(stderr, "No sequence in phylip line %s", buf);
	    exit(1);
	  }
	  if (iSeq >= nSeq) {
	    fprintf(stderr, "No empty line between sequence blocks (is the sequence count wrong?)\n");
	    exit(1);
	  }
	  if (names[iSeq] == NULL) {
	    /* save the name */
	    names[iSeq] = (char *)mymalloc(j+1);
	    int k;
	    for (k = 0; k < j; k++) names[iSeq][k] = buf[k];
	    names[iSeq][j] = '\0';
	  } else {
	    /* check the name */
	    int k;
	    int match = 1;
	    for (k = 0; k < j; k++) {
	      if (names[iSeq][k] != buf[k]) {
		match = 0;
		break;
	      }
	    }
	    if (!match || names[iSeq][j] != '\0') {
	      fprintf(stderr, "Wrong name in phylip line %s\nExpected %s\n", buf, names[iSeq]);
	      exit(1);
	    }
	  }
	}
	int seqlen = strlen(seqs[iSeq]);
	for (; buf[j] != '\n' && buf[j] != '\0'; j++) {
	  if (buf[j] != ' ') {
	    if (seqlen >= nPos) {
	      fprintf(stderr, "Too many characters (expected %d) for sequence named %s\nSo far have:\n%s\n",
		      nPos, names[iSeq], seqs[iSeq]);
	      exit(1);
	    }
	    seqs[iSeq][seqlen++] = toupper(buf[j]);
	  }
	}
	seqs[iSeq][seqlen] = '\0'; /* null-terminate */
	if(verbose>10) fprintf(stderr,"Read iSeq %d name %s seqsofar %s\n", iSeq, names[iSeq], seqs[iSeq]);
	iSeq++;
	if (iSeq == nSeq && strlen(seqs[0]) == nPos)
	  break; /* finished alignment */
      } /* end else non-empty phylip line */
    }
    if (iSeq != nSeq && iSeq != 0) {
      fprintf(stderr, "Wrong number of sequences: expected %d\n", nSeq);
      exit(1);
    }
  }
  /* Check lengths of sequences */
  int i;
  for (i = 0; i < nSeq; i++) {
    int seqlen = strlen(seqs[i]);
    if (seqlen != nPos) {
      fprintf(stderr, "Wrong number of characters for %s: expected %d but have %d instead.\n"
	      "This sequence may be truncated, or another sequence may be too long.\n",
	      names[i], nPos, seqlen);
      exit(1);
    }
  }
  /* Replace "." with "-" and warn if we find any */
  /* If nucleotide sequences, replace U with T and N with X */
  bool findDot = false;
  for (i = 0; i < nSeq; i++) {
    char *p;
    for (p = seqs[i]; *p != '\0'; p++) {
      if (*p == '.') {
	findDot = true;
	*p = '-';
      }
      if (nCodes == 4 && *p == 'U')
	*p = 'T';
      if (nCodes == 4 && *p == 'N')
	*p = 'X';
    }
  }
  if (findDot)
    fprintf(stderr, "Warning! Found \".\" character(s). These are treated as gaps\n");

  if (ferror(fp)) {
    fprintf(stderr, "Error reading input file\n");
    exit(1);
  }

  alignment_t *align = (alignment_t*)mymalloc(sizeof(alignment_t));
  align->nSeq = nSeq;
  align->nPos = nPos;
  align->names = names;
  align->seqs = seqs;
  align->nSaved = nSaved;
  return(align);
}

void FreeAlignmentSeqs(/*IN/OUT*/alignment_t *aln) {
  assert(aln != NULL);
  int i;
  for (i = 0; i < aln->nSeq; i++)
    aln->seqs[i] = myfree(aln->seqs[i], aln->nPos+1);
}

alignment_t *FreeAlignment(alignment_t *aln) {
  if(aln==NULL)
    return(NULL);
  int i;
  for (i = 0; i < aln->nSeq; i++) {
    aln->names[i] = myfree(aln->names[i],strlen(aln->names[i])+1);
    aln->seqs[i] = myfree(aln->seqs[i], aln->nPos+1);
  }
  aln->names = myfree(aln->names, sizeof(char*)*aln->nSaved);
  aln->seqs = myfree(aln->seqs, sizeof(char*)*aln->nSaved);
  myfree(aln, sizeof(alignment_t));
  return(NULL);
}

char **AlnToConstraints(alignment_t *constraints, uniquify_t *unique, hashstrings_t *hashnames) {
  /* look up constraints as names and map to unique-space */
  char **  uniqConstraints = (char**)mymalloc(sizeof(char*) * unique->nUnique);	
  int i;
  for (i = 0; i < unique->nUnique; i++)
    uniqConstraints[i] = NULL;
  for (i = 0; i < constraints->nSeq; i++) {
    char *name = constraints->names[i];
    char *constraintSeq = constraints->seqs[i];
    hashiterator_t hi = FindMatch(hashnames,name);
    if (HashCount(hashnames,hi) != 1) {
      fprintf(stderr, "Sequence %s from constraints file is not in the alignment\n", name);
      exit(1);
    }
    int iSeqNonunique = HashFirst(hashnames,hi);
    assert(iSeqNonunique >= 0 && iSeqNonunique < unique->nSeq);
    int iSeqUnique = unique->alnToUniq[iSeqNonunique];
    assert(iSeqUnique >= 0 && iSeqUnique < unique->nUnique);
    if (uniqConstraints[iSeqUnique] != NULL) {
      /* Already set a constraint for this group of sequences!
	 Warn that we are ignoring this one unless the constraints match */
      if (strcmp(uniqConstraints[iSeqUnique],constraintSeq) != 0) {
	fprintf(stderr,
		"Warning: ignoring constraints for %s:\n%s\n"
		"Another sequence has the same sequence but different constraints\n",
		name, constraintSeq);
      }
    } else {
      uniqConstraints[iSeqUnique] = constraintSeq;
    }
  }
  return(uniqConstraints);
}


profile_t *SeqToProfile(/*IN/OUT*/NJ_t *NJ,
			char *seq, int nPos,
			/*OPTIONAL*/char *constraintSeq, int nConstraints,
			int iNode,
			unsigned long counts[256]) {
  static unsigned char charToCode[256];
  static int codeSet = 0;
  int c, i;

  if (!codeSet) {
    for (c = 0; c < 256; c++) {
      charToCode[c] = nCodes;
    }
    for (i = 0; codesString[i]; i++) {
      charToCode[codesString[i]] = i;
      charToCode[tolower(codesString[i])] = i;
    }
    charToCode['-'] = NOCODE;
    codeSet=1;
  }

  assert(strlen(seq) == nPos);
  profile_t *profile = NewProfile(nPos,nConstraints);

  for (i = 0; i < nPos; i++) {
    unsigned int character = (unsigned int) seq[i];
    counts[character]++;
    c = charToCode[character];
    if(verbose>10 && i < 2) fprintf(stderr,"pos %d char %c code %d\n", i, seq[i], c);
    /* treat unknowns as gaps */
    if (c == nCodes || c == NOCODE) {
      profile->codes[i] = NOCODE;
      profile->weights[i] = 0.0;
    } else {
      profile->codes[i] = c;
      profile->weights[i] = 1.0;
    }
  }
  if (nConstraints > 0) {
    for (i = 0; i < nConstraints; i++) {
      profile->nOn[i] = 0;
      profile->nOff[i] = 0;
    }
    bool bWarn = false;
    if (constraintSeq != NULL) {
      assert(strlen(constraintSeq) == nConstraints);
      for (i = 0; i < nConstraints; i++) {
	if (constraintSeq[i] == '1') {
	  profile->nOn[i] = 1;
	} else if (constraintSeq[i] == '0') {
	  profile->nOff[i] = 1;
	} else if (constraintSeq[i] != '-') {
	  if (!bWarn) {
	    fprintf(stderr, "Constraint characters in unique sequence %d replaced with gap:", iNode+1);
	    bWarn = true;
	  }
	  fprintf(stderr, " %c%d", constraintSeq[i], i+1);
	  /* For the benefit of ConstraintSequencePenalty -- this is a bit of a hack, as
	     this modifies the value read from the alignment
	  */
	  constraintSeq[i] = '-';
	}
      }
      if (bWarn)
	fprintf(stderr, "\n");
    }
  }
  return profile;
}

void SeqDist(unsigned char *codes1, unsigned char *codes2, int nPos,
	     distance_matrix_t *dmat, 
	     /*OUT*/besthit_t *hit) {
  double top = 0;		/* summed over positions */
  int nUse = 0;
  int i;
  if (dmat==NULL) {
    int nDiff = 0;
    for (i = 0; i < nPos; i++) {
      if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
	nUse++;
	if (codes1[i] != codes2[i]) nDiff++;
      }
    }
    top = (double)nDiff;
  } else {
    for (i = 0; i < nPos; i++) {
      if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
	nUse++;
	top += dmat->distances[(unsigned int)codes1[i]][(unsigned int)codes2[i]];
      }
    }
  }
  hit->weight = (double)nUse;
  hit->dist = nUse > 0 ? top/(double)nUse : 1.0;
  seqOps++;
}

void CorrectedPairDistances(profile_t **profiles, int nProfiles,
			    /*OPTIONAL*/distance_matrix_t *distance_matrix,
			    int nPos,
			    /*OUT*/double *distances) {
  assert(distances != NULL);
  assert(profiles != NULL);
  assert(nProfiles>1 && nProfiles <= 4);
  besthit_t hit[6];
  int iHit,i,j;

  for (iHit=0, i=0; i < nProfiles; i++) {
    for (j=i+1; j < nProfiles; j++, iHit++) {
      ProfileDist(profiles[i],profiles[j],nPos,distance_matrix,/*OUT*/&hit[iHit]);
      distances[iHit] = hit[iHit].dist;
    }
  }
  if (pseudoWeight > 0) {
    /* Estimate the prior distance */
    double dTop = 0;
    double dBottom = 0;
    for (iHit=0; iHit < (nProfiles*(nProfiles-1))/2; iHit++) {
      dTop += hit[iHit].dist * hit[iHit].weight;
      dBottom += hit[iHit].weight;
    }
    double prior = (dBottom > 0.01) ? dTop/dBottom : 3.0;
    for (iHit=0; iHit < (nProfiles*(nProfiles-1))/2; iHit++)
      distances[iHit] = (distances[iHit] * hit[iHit].weight + prior * pseudoWeight)
	/ (hit[iHit].weight + pseudoWeight);
  }
  if (logdist) {
    for (iHit=0; iHit < (nProfiles*(nProfiles-1))/2; iHit++)
      distances[iHit] = LogCorrect(distances[iHit]);
  }
}

/* During the neighbor-joining phase, a join only violates our constraints if
   node1, node2, and other are all represented in the constraint
   and if one of the 3 is split and the other two do not agree
 */
int JoinConstraintPenalty(/*IN*/NJ_t *NJ, int node1, int node2) {
  if (NJ->nConstraints == 0)
    return(0.0);
  int penalty = 0;
  int iC;
  for (iC = 0; iC < NJ->nConstraints; iC++)
    penalty += JoinConstraintPenaltyPiece(NJ, node1, node2, iC);
  return(penalty);
}

int JoinConstraintPenaltyPiece(NJ_t *NJ, int node1, int node2, int iC) {
  profile_t *pOut = NJ->outprofile;
  profile_t *p1 = NJ->profiles[node1];
  profile_t *p2 = NJ->profiles[node2];
  int nOn1 = p1->nOn[iC];
  int nOff1 = p1->nOff[iC];
  int nOn2 = p2->nOn[iC];
  int nOff2 = p2->nOff[iC];
  int nOnOut = pOut->nOn[iC] - nOn1 - nOn2;
  int nOffOut = pOut->nOff[iC] - nOff1 - nOff2;

  if ((nOn1+nOff1) > 0 && (nOn2+nOff2) > 0 && (nOnOut+nOffOut) > 0) {
    /* code is -1 for split, 0 for off, 1 for on */
    int code1 = (nOn1 > 0 && nOff1 > 0) ? -1 : (nOn1 > 0 ? 1 : 0);
    int code2 = (nOn2 > 0 && nOff2 > 0) ? -1 : (nOn2 > 0 ? 1 : 0);
    int code3 = (nOnOut > 0 && nOffOut) > 0 ? -1 : (nOnOut > 0 ? 1 : 0);
    int nSplit = (code1 == -1 ? 1 : 0) + (code2 == -1 ? 1 : 0) + (code3 == -1 ? 1 : 0);
    int nOn = (code1 == 1 ? 1 : 0) + (code2 == 1 ? 1 : 0) + (code3 == 1 ? 1 : 0);
    if (nSplit == 1 && nOn == 1)
      return(SplitConstraintPenalty(nOn1+nOn2, nOff1+nOff2, nOnOut, nOffOut));
  }
  /* else */
  return(0);
}

void QuartetConstraintPenalties(profile_t *profiles[4], int nConstraints, /*OUT*/double penalty[3]) {
  int i;
  for (i=0; i < 3; i++)
    penalty[i] = 0.0;
  if(nConstraints == 0)
    return;
  int iC;
  for (iC = 0; iC < nConstraints; iC++) {
    double part[3];
    if (QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/part)) {
      for (i=0;i<3;i++)
	penalty[i] += part[i];

      if (verbose>2
	  && (fabs(part[ABvsCD]-part[ACvsBD]) > 0.001 || fabs(part[ABvsCD]-part[ADvsBC]) > 0.001))
	fprintf(stderr, "Constraint Penalties at %d: ABvsCD %.3f ACvsBD %.3f ADvsBC %.3f %d/%d %d/%d %d/%d %d/%d\n",
		iC, part[ABvsCD], part[ACvsBD], part[ADvsBC],
		profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
    }
  }
  if (verbose>2)
    fprintf(stderr, "Total Constraint Penalties: ABvsCD %.3f ACvsBD %.3f ADvsBC %.3f\n",
	    penalty[ABvsCD], penalty[ACvsBD], penalty[ADvsBC]);
}

double PairConstraintDistance(int nOn1, int nOff1, int nOn2, int nOff2) {
  double f1 = nOn1/(double)(nOn1+nOff1);
  double f2 = nOn2/(double)(nOn2+nOff2);
  /* 1 - f1 * f2 - (1-f1)*(1-f2) = 1 - f1 * f2 - 1 + f1 + f2 - f1 * f2 */
  return(f1 + f2 - 2.0 * f1 * f2);
}

bool QuartetConstraintPenaltiesPiece(profile_t *profiles[4], int iC, /*OUT*/double piece[3]) {
  int nOn[4];
  int nOff[4];
  int i;
  int nSplit = 0;
  int nPlus = 0;
  int nMinus = 0;
  
  for (i=0; i < 4; i++) {
    nOn[i] = profiles[i]->nOn[iC];
    nOff[i] = profiles[i]->nOff[iC];
    if (nOn[i] + nOff[i] == 0)
      return(false);		/* ignore */
    else if (nOn[i] > 0 && nOff[i] > 0)
      nSplit++;
    else if (nOn[i] > 0)
      nPlus++;
    else
      nMinus++;
  }
  /* If just one of them is split or on the other side and the others all agree, also ignore */
  if (nPlus >= 3 || nMinus >= 3)
    return(false);
  piece[ABvsCD] = constraintWeight
    * (PairConstraintDistance(nOn[0],nOff[0],nOn[1],nOff[1])
       + PairConstraintDistance(nOn[2],nOff[2],nOn[3],nOff[3]));
  piece[ACvsBD] = constraintWeight
    * (PairConstraintDistance(nOn[0],nOff[0],nOn[2],nOff[2])
       + PairConstraintDistance(nOn[1],nOff[1],nOn[3],nOff[3]));
  piece[ADvsBC] = constraintWeight
    * (PairConstraintDistance(nOn[0],nOff[0],nOn[3],nOff[3])
       + PairConstraintDistance(nOn[2],nOff[2],nOn[1],nOff[1]));
  return(true);
}

/* Minimum number of constrained leaves that need to be moved
   to satisfy the constraint (or 0 if constraint is satisfied)
   Defining it this way should ensure that SPR moves that break
   constraints get a penalty
*/
int SplitConstraintPenalty(int nOn1, int nOff1, int nOn2, int nOff2) {
  return(nOn1 + nOff2 < nOn2 + nOff1 ?
	 (nOn1 < nOff2 ? nOn1 : nOff2)
	 : (nOn2 < nOff1 ? nOn2 : nOff1));
}

bool SplitViolatesConstraint(profile_t *profiles[4], int iConstraint) {
  int i;
  int codes[4]; /* 0 for off, 1 for on, -1 for split (quit if not constrained at all) */
  for (i = 0; i < 4; i++) {
    if (profiles[i]->nOn[iConstraint] + profiles[i]->nOff[iConstraint] == 0)
      return(false);
    else if (profiles[i]->nOn[iConstraint] > 0 && profiles[i]->nOff[iConstraint] == 0)
      codes[i] = 1;
    else if (profiles[i]->nOn[iConstraint] == 0 && profiles[i]->nOff[iConstraint] > 0)
      codes[i] = 0;
    else
      codes[i] = -1;
  }
  int n0 = 0;
  int n1 = 0;
  for (i = 0; i < 4; i++) {
    if (codes[i] == 0)
      n0++;
    else if (codes[i] == 1)
      n1++;
  }
  /* 3 on one side means no violation, even if other is code -1
     otherwise must have code != -1 and agreement on the split
   */
  if (n0 >= 3 || n1 >= 3)
    return(false);
  if (n0==2 && n1==2 && codes[0] == codes[1] && codes[2] == codes[3])
    return(false);
  return(true);
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

/* A helper function -- f1 and f2 can be NULL if the corresponding code != NOCODE
*/
double ProfileDistPiece(unsigned int code1, unsigned int code2,
			numeric_t *f1, numeric_t *f2, 
			/*OPTIONAL*/distance_matrix_t *dmat,
			/*OPTIONAL*/numeric_t *codeDist2) {
  if (dmat) {
    if (code1 != NOCODE && code2 != NOCODE) { /* code1 vs code2 */
      return(dmat->distances[code1][code2]);
    } else if (codeDist2 != NULL && code1 != NOCODE) { /* code1 vs. codeDist2 */
      return(codeDist2[code1]);
    } else { /* f1 vs f2 */
      if (f1 == NULL) {
	if(code1 == NOCODE) return(10.0);
	f1 = &dmat->codeFreq[code1][0];
      }
      if (f2 == NULL) {
	if(code2 == NOCODE) return(10.0);
	f2 = &dmat->codeFreq[code2][0];
      }
      return(vector_multiply3_sum(f1,f2,dmat->eigenval,nCodes));
    }
  } else {
    /* no matrix */
    if (code1 != NOCODE) {
      if (code2 != NOCODE) {
	return(code1 == code2 ? 0.0 : 1.0); /* code1 vs code2 */
      } else {
	if(f2 == NULL) return(10.0);
	return(1.0 - f2[code1]); /* code1 vs. f2 */
      }
    } else {
      if (code2 != NOCODE) {
	if(f1 == NULL) return(10.0);
	return(1.0 - f1[code2]); /* f1 vs code2 */
      } else { /* f1 vs. f2 */
	if (f1 == NULL || f2 == NULL) return(10.0);
	double piece = 1.0;
	int k;
	for (k = 0; k < nCodes; k++) {
	  piece -= f1[k] * f2[k];
	}
	return(piece);
      }
    }
  }
  assert(0);
}

/* E.g. GET_FREQ(profile,iPos,iVector)
   Gets the next element of the vectors (and updates iVector), or
   returns NULL if we didn't store a vector
*/
#define GET_FREQ(P,I,IVECTOR) \
(P->weights[I] > 0 && P->codes[I] == NOCODE ? &P->vectors[nCodes*(IVECTOR++)] : NULL)

void ProfileDist(profile_t *profile1, profile_t *profile2, int nPos,
		 /*OPTIONAL*/distance_matrix_t *dmat,
		 /*OUT*/besthit_t *hit) {
  double top = 0;
  double denom = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  int i = 0;
  for (i = 0; i < nPos; i++) {
      numeric_t *f1 = GET_FREQ(profile1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(profile2,i,/*IN/OUT*/iFreq2);
      if (profile1->weights[i] > 0 && profile2->weights[i] > 0) {
	double weight = profile1->weights[i] * profile2->weights[i];
	denom += weight;
	double piece = ProfileDistPiece(profile1->codes[i],profile2->codes[i],f1,f2,dmat,
					profile2->codeDist ? &profile2->codeDist[i*nCodes] : NULL);
	top += weight * piece;
      }
  }
  assert(iFreq1 == profile1->nVectors);
  assert(iFreq2 == profile2->nVectors);
  hit->weight = denom > 0 ? denom : 0.01; /* 0.01 is an arbitrarily low value of weight (normally >>1) */
  hit->dist = denom > 0 ? top/denom : 1;
  profileOps++;
}

/* This should not be called if the update weight is 0, as
   in that case code==NOCODE and in=NULL is possible, and then
   it will fail.
*/
void AddToFreq(/*IN/OUT*/numeric_t *fOut,
	       double weight,
	       unsigned int codeIn, /*OPTIONAL*/numeric_t *fIn,
	       /*OPTIONAL*/distance_matrix_t *dmat) {
  assert(fOut != NULL);
  if (fIn != NULL) {
    vector_add_mult(fOut, fIn, weight, nCodes);
  } else if (dmat) {
    assert(codeIn != NOCODE);
    vector_add_mult(fOut, dmat->codeFreq[codeIn], weight, nCodes);
  } else {
    assert(codeIn != NOCODE);
    fOut[codeIn] += weight;
  }
}

void SetProfile(/*IN/OUT*/NJ_t *NJ, int node, double weight1) {
    children_t *c = &NJ->child[node];
    assert(c->nChild == 2);
    assert(NJ->profiles[c->child[0]] != NULL);
    assert(NJ->profiles[c->child[1]] != NULL);
    if (NJ->profiles[node] != NULL)
      FreeProfile(NJ->profiles[node], NJ->nPos, NJ->nConstraints);
    NJ->profiles[node] = AverageProfile(NJ->profiles[c->child[0]],
					NJ->profiles[c->child[1]],
					NJ->nPos, NJ->nConstraints,
					NJ->distance_matrix,
					weight1);
}

/* bionjWeight is the weight of the first sequence (between 0 and 1),
   or -1 to do the average.
   */
profile_t *AverageProfile(profile_t *profile1, profile_t *profile2,
			  int nPos, int nConstraints,
			  distance_matrix_t *dmat,
			  double bionjWeight) {
  int i;
  if (bionjWeight < 0) {
    bionjWeight = 0.5;
  }

  /* First, set codes and weights and see how big vectors will be */
  profile_t *out = NewProfile(nPos, nConstraints);

  for (i = 0; i < nPos; i++) {
    out->weights[i] = bionjWeight * profile1->weights[i]
      + (1-bionjWeight) * profile2->weights[i];
    out->codes[i] = NOCODE;
    if (out->weights[i] > 0) {
      if (profile1->weights[i] > 0 && profile1->codes[i] != NOCODE
	  && (profile2->weights[i] <= 0 || profile1->codes[i] == profile2->codes[i])) {
	out->codes[i] = profile1->codes[i];
      } else if (profile1->weights[i] <= 0
		 && profile2->weights[i] > 0
		 && profile2->codes[i] != NOCODE) {
	out->codes[i] = profile2->codes[i];
      }
      if (out->codes[i] == NOCODE) out->nVectors++;
    }
  }

  /* Allocate and set the vectors */
  out->vectors = (numeric_t*)mymalloc(sizeof(numeric_t)*nCodes*out->nVectors);
  for (i = 0; i < nCodes * out->nVectors; i++) out->vectors[i] = 0;
  nProfileFreqAlloc += out->nVectors;
  nProfileFreqAvoid += nPos - out->nVectors;
  int iFreqOut = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  for (i=0; i < nPos; i++) {
    numeric_t *f = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
    numeric_t *f1 = GET_FREQ(profile1,i,/*IN/OUT*/iFreq1);
    numeric_t *f2 = GET_FREQ(profile2,i,/*IN/OUT*/iFreq2);
    if (f != NULL) {
      if (profile1->weights[i] > 0)
	AddToFreq(/*IN/OUT*/f, profile1->weights[i] * bionjWeight,
		  profile1->codes[i], f1, dmat);
      if (profile2->weights[i] > 0)
	AddToFreq(/*IN/OUT*/f, profile2->weights[i] * (1.0-bionjWeight),
		  profile2->codes[i], f2, dmat);
      NormalizeFreq(/*IN/OUT*/f, dmat);
    } /* end if computing f */
    if (verbose > 10 && i < 5) {
      fprintf(stderr,"Average profiles: pos %d in-w1 %f in-w2 %f bionjWeight %f to weight %f code %d\n",
	      i, profile1->weights[i], profile2->weights[i], bionjWeight,
	      out->weights[i], out->codes[i]);
      if (f!= NULL) {
	int k;
	for (k = 0; k < nCodes; k++)
	  fprintf(stderr, "\t%c:%f", codesString[k], f ? f[k] : -1.0);
	fprintf(stderr,"\n");
      }
    }
  } /* end loop over positions */
  assert(iFreq1 == profile1->nVectors);
  assert(iFreq2 == profile2->nVectors);
  assert(iFreqOut == out->nVectors);

  /* compute total constraints */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] = profile1->nOn[i] + profile2->nOn[i];
    out->nOff[i] = profile1->nOff[i] + profile2->nOff[i];
  }
  profileAvgOps++;
  return(out);
}

/* Make the (unrotated) frequencies sum to 1
   Simply dividing by total_weight is not ideal because of roundoff error
   So compute total_freq instead
*/
void NormalizeFreq(/*IN/OUT*/numeric_t *freq, distance_matrix_t *dmat) {
  double total_freq = 0;
  int k;
  if (dmat != NULL) {
    /* The total frequency is dot_product(true_frequencies, 1)
       So we rotate the 1 vector by eigeninv (stored in eigentot)
    */
    total_freq = vector_multiply_sum(freq, dmat->eigentot, nCodes);
  } else {
    for (k = 0; k < nCodes; k++)
      total_freq += freq[k];
  }
  if (total_freq > fPostTotalTolerance) {
    numeric_t inverse_weight = 1.0/total_freq;
    vector_multiply_by(/*IN/OUT*/freq, inverse_weight, nCodes);
  } else {
    /* This can happen if we are in a very low-weight region, e.g. if a mostly-gap position gets weighted down
       repeatedly; just set them all to arbitrary but legal values */
    if (dmat == NULL) {
      for (k = 0; k < nCodes; k++)
	freq[k] = 1.0/nCodes;
    } else {
      for (k = 0; k < nCodes; k++)
	freq[k] = dmat->codeFreq[0][k];
    }
  }
}

/* OutProfile() computes the out-profile */
profile_t *OutProfile(profile_t **profiles, int nProfiles,
		      int nPos, int nConstraints,
		      distance_matrix_t *dmat) {
  int i;			/* position */
  int in;			/* profile */
  profile_t *out = NewProfile(nPos, nConstraints);

  double inweight = 1.0/(double)nProfiles;   /* The maximal output weight is 1.0 */

  /* First, set weights -- code is always NOCODE, prevent weight=0 */
  for (i = 0; i < nPos; i++) {
    out->weights[i] = 0;
    for (in = 0; in < nProfiles; in++)
      out->weights[i] += profiles[in]->weights[i] * inweight;
    if (out->weights[i] <= 0) out->weights[i] = 1e-20; /* always store a vector */
    out->nVectors++;
    out->codes[i] = NOCODE;		/* outprofile is normally complicated */
  }

  /* Initialize the frequencies to 0 */
  out->vectors = (numeric_t*)mymalloc(sizeof(numeric_t)*nCodes*out->nVectors);
  for (i = 0; i < nCodes*out->nVectors; i++)
    out->vectors[i] = 0;

  /* Add up the weights, going through each sequence in turn */
  for (in = 0; in < nProfiles; in++) {
    int iFreqOut = 0;
    int iFreqIn = 0;
    for (i = 0; i < nPos; i++) {
      numeric_t *fIn = GET_FREQ(profiles[in],i,/*IN/OUT*/iFreqIn);
      numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      if (profiles[in]->weights[i] > 0)
	AddToFreq(/*IN/OUT*/fOut, profiles[in]->weights[i],
		  profiles[in]->codes[i], fIn, dmat);
    }
    assert(iFreqOut == out->nVectors);
    assert(iFreqIn == profiles[in]->nVectors);
  }

  /* And normalize the frequencies to sum to 1 */
  int iFreqOut = 0;
  for (i = 0; i < nPos; i++) {
    numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
    if (fOut)
      NormalizeFreq(/*IN/OUT*/fOut, dmat);
  }
  assert(iFreqOut == out->nVectors);
  if (verbose > 10) fprintf(stderr,"Average %d profiles\n", nProfiles);
  if(dmat)
    SetCodeDist(/*IN/OUT*/out, nPos, dmat);

  /* Compute constraints */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] = 0;
    out->nOff[i] = 0;
    for (in = 0; in < nProfiles; in++) {
      out->nOn[i] += profiles[in]->nOn[i];
      out->nOff[i] += profiles[in]->nOff[i];
    }
  }
  return(out);
}

void UpdateOutProfile(/*IN/OUT*/profile_t *out, profile_t *old1, profile_t *old2,
		      profile_t *new, int nActiveOld,
		      int nPos, int nConstraints,
		      distance_matrix_t *dmat) {
  int i, k;
  int iFreqOut = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  int iFreqNew = 0;
  assert(nActiveOld > 0);

  for (i = 0; i < nPos; i++) {
    numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
    numeric_t *fOld1 = GET_FREQ(old1,i,/*IN/OUT*/iFreq1);
    numeric_t *fOld2 = GET_FREQ(old2,i,/*IN/OUT*/iFreq2);
    numeric_t *fNew = GET_FREQ(new,i,/*IN/OUT*/iFreqNew);

    assert(out->codes[i] == NOCODE && fOut != NULL); /* No no-vector optimization for outprofiles */
    if (verbose > 3 && i < 3) {
      fprintf(stderr,"Updating out-profile position %d weight %f (mult %f)\n",
	      i, out->weights[i], out->weights[i]*nActiveOld);
    }
    double originalMult = out->weights[i]*nActiveOld;
    double newMult = originalMult + new->weights[i] - old1->weights[i] - old2->weights[i];
    out->weights[i] = newMult/(nActiveOld-1);
    if (out->weights[i] <= 0) out->weights[i] = 1e-20; /* always use the vector */

    for (k = 0; k < nCodes; k++) fOut[k] *= originalMult;
    
    if (old1->weights[i] > 0)
      AddToFreq(/*IN/OUT*/fOut, -old1->weights[i], old1->codes[i], fOld1, dmat);
    if (old2->weights[i] > 0)
      AddToFreq(/*IN/OUT*/fOut, -old2->weights[i], old2->codes[i], fOld2, dmat);
    if (new->weights[i] > 0)
      AddToFreq(/*IN/OUT*/fOut, new->weights[i], new->codes[i], fNew, dmat);

    /* And renormalize */
    NormalizeFreq(/*IN/OUT*/fOut, dmat);

    if (verbose > 2 && i < 3) {
      fprintf(stderr,"Updated out-profile position %d weight %f (mult %f)",
	      i, out->weights[i], out->weights[i]*nActiveOld);
      if(out->weights[i] > 0)
	for (k=0;k<nCodes;k++)
	  fprintf(stderr, " %c:%f", dmat?'?':codesString[k], fOut[k]);
      fprintf(stderr,"\n");
    }
  }
  assert(iFreqOut == out->nVectors);
  assert(iFreq1 == old1->nVectors);
  assert(iFreq2 == old2->nVectors);
  assert(iFreqNew == new->nVectors);
  if(dmat)
    SetCodeDist(/*IN/OUT*/out,nPos,dmat);

  /* update constraints -- note in practice this should be a no-op */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] += new->nOn[i] - old1->nOn[i] - old2->nOn[i];
    out->nOff[i] += new->nOff[i] - old1->nOff[i] - old2->nOff[i];
  }
}

void SetCodeDist(/*IN/OUT*/profile_t *profile, int nPos,
			   distance_matrix_t *dmat) {
  if (profile->codeDist == NULL)
    profile->codeDist = (numeric_t*)mymalloc(sizeof(numeric_t)*nPos*nCodes);
  int i;
  int iFreq = 0;
  for (i = 0; i < nPos; i++) {
    numeric_t *f = GET_FREQ(profile,i,/*IN/OUT*/iFreq);

    int k;
    for (k = 0; k < nCodes; k++)
      profile->codeDist[i*nCodes+k] = ProfileDistPiece(/*code1*/profile->codes[i], /*code2*/k,
						       /*f1*/f, /*f2*/NULL,
						       dmat, NULL);
  }
  assert(iFreq==profile->nVectors);
}


void SetBestHit(int node, NJ_t *NJ, int nActive,
		/*OUT*/besthit_t *bestjoin, /*OUT OPTIONAL*/besthit_t *allhits) {
  assert(NJ->parent[node] <  0);

  bestjoin->i = node;
  bestjoin->j = -1;
  bestjoin->dist = 1e20;
  bestjoin->criterion = 1e20;

  int j;
  besthit_t tmp;

#ifdef OPENMP
  /* Note -- if we are already in a parallel region, this will be ignored */
  #pragma omp parallel for schedule(dynamic, 50)
#endif
  for (j = 0; j < NJ->maxnode; j++) {
    besthit_t *sv = allhits != NULL ? &allhits[j] : &tmp;
    sv->i = node;
    sv->j = j;
    if (NJ->parent[j] >= 0) {
      sv->i = -1;		/* illegal/empty join */
      sv->weight = 0.0;
      sv->criterion = sv->dist = 1e20;
      continue;
    }
    /* Note that we compute self-distances (allow j==node) because the top-hit heuristic
       expects self to be within its top hits, but we exclude those from the bestjoin
       that we return...
    */
    SetDistCriterion(NJ, nActive, /*IN/OUT*/sv);
    if (sv->criterion < bestjoin->criterion && node != j)
      *bestjoin = *sv;
  }
  if (verbose>5) {
    fprintf(stderr, "SetBestHit %d %d %f %f\n", bestjoin->i, bestjoin->j, bestjoin->dist, bestjoin->criterion);
  }
}

void ReadMatrix(char *filename, /*OUT*/numeric_t codes[MAXCODES][MAXCODES], bool checkCodes) {
  char buf[BUFFER_SIZE] = "";
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot read %s\n",filename);
    exit(1);
  }
  if (fgets(buf,sizeof(buf),fp) == NULL) {
    fprintf(stderr, "Error reading header line for %s:\n%s\n", filename, buf);
    exit(1);
  }
  if (checkCodes) {
    int i;
    int iBufPos;
    for (iBufPos=0,i=0;i<nCodes;i++,iBufPos++) {
      if(buf[iBufPos] != codesString[i]) {
	fprintf(stderr,"Header line\n%s\nin file %s does not have expected code %c # %d in %s\n",
		buf, filename, codesString[i], i, codesString);
	exit(1);
      }
      iBufPos++;
      if(buf[iBufPos] != '\n' && buf[iBufPos] != '\r' && buf[iBufPos] != '\0' && buf[iBufPos] != '\t') {
	fprintf(stderr, "Header line in %s should be tab-delimited\n", filename);
	exit(1);
      }
      if (buf[iBufPos] == '\0' && i < nCodes-1) {
	fprintf(stderr, "Header line in %s ends prematurely\n",filename);
	exit(1);
      }
    } /* end loop over codes */
    /* Should be at end, but allow \n because of potential DOS \r\n */
    if(buf[iBufPos] != '\0' && buf[iBufPos] != '\n' && buf[iBufPos] != '\r') {
      fprintf(stderr, "Header line in %s has too many entries\n", filename);
      exit(1);
    }
  }
  int iLine;
  for (iLine = 0; iLine < nCodes; iLine++) {
    buf[0] = '\0';
    if (fgets(buf,sizeof(buf),fp) == NULL) {
      fprintf(stderr, "Cannot read line %d from file %s\n", iLine+2, filename);
      exit(1);
    }
    char *field = strtok(buf,"\t\r\n");
    field = strtok(NULL, "\t");	/* ignore first column */
    int iColumn;
    for (iColumn = 0; iColumn < nCodes && field != NULL; iColumn++, field = strtok(NULL,"\t")) {
      if(sscanf(field,ScanNumericSpec,&codes[iLine][iColumn]) != 1) {
	fprintf(stderr,"Cannot parse field %s in file %s\n", field, filename);
	exit(1);
      }
    }
  }
}

void ReadVector(char *filename, /*OUT*/numeric_t codes[MAXCODES]) {
  FILE *fp = fopen(filename,"r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot read %s\n",filename);
    exit(1);
  }
  int i;
  for (i = 0; i < nCodes; i++) {
    if (fscanf(fp,ScanNumericSpec,&codes[i]) != 1) {
      fprintf(stderr,"Cannot read %d entry of %s\n",i+1,filename);
      exit(1);
    }
  }
  if (fclose(fp) != 0) {
    fprintf(stderr, "Error reading %s\n",filename);
    exit(1);
  }
}

distance_matrix_t *ReadDistanceMatrix(char *prefix) {
  char buffer[BUFFER_SIZE];
  distance_matrix_t *dmat = (distance_matrix_t*)mymalloc(sizeof(distance_matrix_t));

  if(strlen(prefix) > BUFFER_SIZE-20) {
    fprintf(stderr,"Filename %s too long\n", prefix);
    exit(1);
  }

  strcpy(buffer, prefix);
  strcat(buffer, ".distances");
  ReadMatrix(buffer, /*OUT*/dmat->distances, /*checkCodes*/true);

  strcpy(buffer, prefix);
  strcat(buffer, ".inverses");
  ReadMatrix(buffer, /*OUT*/dmat->eigeninv, /*checkCodes*/false);

  strcpy(buffer, prefix);
  strcat(buffer, ".eigenvalues");
  ReadVector(buffer, /*OUT*/dmat->eigenval);

  if(verbose>1) fprintf(stderr, "Read distance matrix from %s\n",prefix);
  SetupDistanceMatrix(/*IN/OUT*/dmat);
  return(dmat);
}

void SetupDistanceMatrix(/*IN/OUT*/distance_matrix_t *dmat) {
  /* Check that the eigenvalues and eigen-inverse are consistent with the
     distance matrix and that the matrix is symmetric */
  int i,j,k;
  for (i = 0; i < nCodes; i++) {
    for (j = 0; j < nCodes; j++) {
      if(fabs(dmat->distances[i][j]-dmat->distances[j][i]) > 1e-6) {
	fprintf(stderr,"Distance matrix not symmetric for %d,%d: %f vs %f\n",
		i+1,j+1,
		dmat->distances[i][j],
		dmat->distances[j][i]);
	exit(1);
      }
      double total = 0.0;
      for (k = 0; k < nCodes; k++)
	total += dmat->eigenval[k] * dmat->eigeninv[k][i] * dmat->eigeninv[k][j];
      if(fabs(total - dmat->distances[i][j]) > 1e-6) {
	fprintf(stderr,"Distance matrix entry %d,%d should be %f but eigen-representation gives %f\n",
		i+1,j+1,dmat->distances[i][j],total);
	exit(1);
      }
    }
  }
  
  /* And compute eigentot */
  for (k = 0; k < nCodes; k++) {
    dmat->eigentot[k] = 0.;
    int j;
    for (j = 0; j < nCodes; j++)
      dmat->eigentot[k] += dmat->eigeninv[k][j];
  }
  
  /* And compute codeFreq */
  int code;
  for(code = 0; code < nCodes; code++) {
    for (k = 0; k < nCodes; k++) {
      dmat->codeFreq[code][k] = dmat->eigeninv[k][code];
    }
  }
  /* And gapFreq */
  for(code = 0; code < nCodes; code++) {
    double gapFreq = 0.0;
    for (k = 0; k < nCodes; k++)
      gapFreq += dmat->codeFreq[k][code];
    dmat->gapFreq[code] = gapFreq / nCodes;
  }

  if(verbose>10) fprintf(stderr, "Made codeFreq\n");
}

nni_t ChooseNNI(profile_t *profiles[4],
		/*OPTIONAL*/distance_matrix_t *dmat,
		int nPos, int nConstraints,
		/*OUT*/double criteria[3]) {
  double d[6];
  CorrectedPairDistances(profiles, 4, dmat, nPos, /*OUT*/d);
  double penalty[3]; 		/* indexed as nni_t */
  QuartetConstraintPenalties(profiles, nConstraints, /*OUT*/penalty);
  criteria[ABvsCD] = d[qAB] + d[qCD] + penalty[ABvsCD];
  criteria[ACvsBD] = d[qAC] + d[qBD] + penalty[ACvsBD];
  criteria[ADvsBC] = d[qAD] + d[qBC] + penalty[ADvsBC];

  nni_t choice = ABvsCD;
  if (criteria[ACvsBD] < criteria[ABvsCD] && criteria[ACvsBD] <= criteria[ADvsBC]) {
    choice = ACvsBD;
  } else if (criteria[ADvsBC] < criteria[ABvsCD] && criteria[ADvsBC] <= criteria[ACvsBD]) {
    choice = ADvsBC;
  }
  if (verbose > 1 && penalty[choice] > penalty[ABvsCD] + 1e-6) {
    fprintf(stderr, "Worsen constraint: from %.3f to %.3f distance %.3f to %.3f: ",
	    penalty[ABvsCD], penalty[choice],
	    criteria[ABvsCD], choice == ACvsBD ? criteria[ACvsBD] : criteria[ADvsBC]);
    int iC;
    for (iC = 0; iC < nConstraints; iC++) {
      double ppart[3];
      if (QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/ppart)) {
	double old_penalty = ppart[ABvsCD];
	double new_penalty = ppart[choice];
	if (new_penalty > old_penalty + 1e-6)
	  fprintf(stderr, " %d (%d/%d %d/%d %d/%d %d/%d)", iC,
		  profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		  profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		  profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		  profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
      }
    }
    fprintf(stderr,"\n");
  }
  if (verbose > 3)
    fprintf(stderr, "NNI scores ABvsCD %.5f ACvsBD %.5f ADvsBC %.5f choice %s\n",
	    criteria[ABvsCD], criteria[ACvsBD], criteria[ADvsBC],
	    choice == ABvsCD ? "AB|CD" : (choice == ACvsBD ? "AC|BD" : "AD|BC"));
  return(choice);
}

profile_t *PosteriorProfile(profile_t *p1, profile_t *p2,
			    double len1, double len2,
			    /*OPTIONAL*/transition_matrix_t *transmat,
			    rates_t *rates,
			    int nPos, int nConstraints) {
  if (len1 < MLMinBranchLength)
    len1 = MLMinBranchLength;
  if (len2 < MLMinBranchLength)
    len2 = MLMinBranchLength;

  int i,j,k;
  profile_t *out = NewProfile(nPos, nConstraints);
  for (i = 0; i < nPos; i++) {
    out->codes[i] = NOCODE;
    out->weights[i] = 1.0;
  }
  out->nVectors = nPos;
  out->vectors = (numeric_t*)mymalloc(sizeof(numeric_t)*nCodes*out->nVectors);
  for (i = 0; i < nCodes * out->nVectors; i++) out->vectors[i] = 0;
  int iFreqOut = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  numeric_t *expeigenRates1 = NULL, *expeigenRates2 = NULL;

  if (transmat != NULL) {
    expeigenRates1 = ExpEigenRates(len1, transmat, rates);
    expeigenRates2 = ExpEigenRates(len2, transmat, rates);
  }

  if (transmat == NULL) {	/* Jukes-Cantor */
    assert(nCodes == 4);

    double *PSame1 = PSameVector(len1, rates);
    double *PDiff1 = PDiffVector(PSame1, rates);
    double *PSame2 = PSameVector(len2, rates);
    double *PDiff2 = PDiffVector(PSame2, rates);

    numeric_t mix1[4], mix2[4];

    for (i=0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      double w1 = p1->weights[i];
      double w2 = p2->weights[i];
      int code1 = p1->codes[i];
      int code2 = p2->codes[i];
      numeric_t *f1 = GET_FREQ(p1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(p2,i,/*IN/OUT*/iFreq2);

      /* First try to store a simple profile */
      if (f1 == NULL && f2 == NULL) {
	if (code1 == NOCODE && code2 == NOCODE) {
	  out->codes[i] = NOCODE;
	  out->weights[i] = 0.0;
	  continue;
	} else if (code1 == NOCODE) {
	  /* Posterior(parent | character & gap, len1, len2) = Posterior(parent | character, len1)
	     = PSame() for matching characters and 1-PSame() for the rest
	     = (pSame - pDiff) * character + (1-(pSame-pDiff)) * gap
	  */
	  out->codes[i] = code2;
	  out->weights[i] = w2 * (PSame2[iRate] - PDiff2[iRate]);
	  continue;
	} else if (code2 == NOCODE) {
	  out->codes[i] = code1;
	  out->weights[i] = w1 * (PSame1[iRate] - PDiff1[iRate]);
	  continue;
	} else if (code1 == code2) {
	  out->codes[i] = code1;
	  double f12code = (w1*PSame1[iRate] + (1-w1)*0.25) * (w2*PSame2[iRate] + (1-w2)*0.25);
	  double f12other = (w1*PDiff1[iRate] + (1-w1)*0.25) * (w2*PDiff2[iRate] + (1-w2)*0.25);
	  /* posterior probability of code1/code2 after scaling */
	  double pcode = f12code/(f12code+3*f12other);
	  /* Now f = w * (code ? 1 : 0) + (1-w) * 0.25, so to get pcode we need
	     fcode = 1/4 + w1*3/4 or w = (f-1/4)*4/3
	   */
	  out->weights[i] = (pcode - 0.25) * 4.0/3.0;
	  /* This can be zero because of numerical problems, I think */
	  if (out->weights[i] < 1e-6) {
	    if (verbose > 1)
	      fprintf(stderr, "Replaced weight %f with %f from w1 %f w2 %f PSame %f %f f12code %f f12other %f\n",
		      out->weights[i], 1e-6,
		      w1, w2,
		      PSame1[iRate], PSame2[iRate],
		      f12code, f12other);
	    out->weights[i] = 1e-6;
	  }
	  continue;
	}
      }
      /* if we did not compute a simple profile, then do the full computation and
         store the full vector
      */
      if (f1 == NULL) {
	for (j = 0; j < 4; j++)
	  mix1[j] = (1-w1)*0.25;
	if(code1 != NOCODE)
	  mix1[code1] += w1;
	f1 = mix1;
      }
      if (f2 == NULL) {
	for (j = 0; j < 4; j++)
	  mix2[j] = (1-w2)*0.25;
	if(code2 != NOCODE)
	  mix2[code2] += w2;
	f2 = mix2;
      }
      out->codes[i] = NOCODE;
      out->weights[i] = 1.0;
      numeric_t *f = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      double lkAB = 0;
      for (j = 0; j < 4; j++) {
	f[j] = (f1[j] * PSame1[iRate] + (1.0-f1[j]) * PDiff1[iRate])
	  * (f2[j] * PSame2[iRate] + (1.0-f2[j]) * PDiff2[iRate]);
	lkAB += f[j];
      }
      double lkABInv = 1.0/lkAB;
      for (j = 0; j < 4; j++)
	f[j] *= lkABInv;
    }
    PSame1 = myfree(PSame1, sizeof(double) * rates->nRateCategories);
    PSame2 = myfree(PSame2, sizeof(double) * rates->nRateCategories);
    PDiff1 = myfree(PDiff1, sizeof(double) * rates->nRateCategories);
    PDiff2 = myfree(PDiff2, sizeof(double) * rates->nRateCategories);
  } else if (nCodes == 4) {	/* matrix model on nucleotides */
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];
    numeric_t f1mix[4], f2mix[4];
    
    for (i=0; i < nPos; i++) {
      if (p1->codes[i] == NOCODE && p2->codes[i] == NOCODE
	  && p1->weights[i] == 0 && p2->weights[i] == 0) {
	/* aligning gap with gap -- just output a gap
	   out->codes[i] is already set to NOCODE so need not set that */
	out->weights[i] = 0;
	continue;
      }
      int iRate = rates->ratecat[i];
      numeric_t *expeigen1 = &expeigenRates1[iRate*4];
      numeric_t *expeigen2 = &expeigenRates2[iRate*4];
      numeric_t *f1 = GET_FREQ(p1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(p2,i,/*IN/OUT*/iFreq2);
      numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      assert(fOut != NULL);

      if (f1 == NULL) {
	f1 = &transmat->codeFreq[p1->codes[i]][0]; /* codeFreq includes an entry for NOCODE */
	double w = p1->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 4; j++)
	    f1mix[j] = w * f1[j] + (1.0-w) * fGap[j];
	  f1 = f1mix;
	}
      }
      if (f2 == NULL) {
	f2 = &transmat->codeFreq[p2->codes[i]][0];
	double w = p2->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 4; j++)
	    f2mix[j] = w * f2[j] + (1.0-w) * fGap[j];
	  f2 = f2mix;
	}
      }
      numeric_t fMult1[4] ALIGNED;	/* rotated1 * expeigen1 */
      numeric_t fMult2[4] ALIGNED;	/* rotated2 * expeigen2 */
#if 0 /* SSE3 is slower */
      vector_multiply(f1, expeigen1, 4, /*OUT*/fMult1);
      vector_multiply(f2, expeigen2, 4, /*OUT*/fMult2);
#else
      for (j = 0; j < 4; j++) {
	fMult1[j] = f1[j]*expeigen1[j];
	fMult2[j] = f2[j]*expeigen2[j];
      }
#endif
      numeric_t fPost[4] ALIGNED;		/* in  unrotated space */
      for (j = 0; j < 4; j++) {
#if 0 /* SSE3 is slower */
	fPost[j] = vector_dot_product_rot(fMult1, fMult2, &transmat->codeFreq[j][0], 4)
	  * transmat->statinv[j]; */
#else
	double out1 = 0;
	double out2 = 0;
	for (k = 0; k < 4; k++) {
	  out1 += fMult1[k] * transmat->codeFreq[j][k];
	  out2 += fMult2[k] * transmat->codeFreq[j][k];
	}
	fPost[j] = out1*out2*transmat->statinv[j];
#endif
      }
      double fPostTot = 0;
      for (j = 0; j < 4; j++)
	fPostTot += fPost[j];
      assert(fPostTot > fPostTotalTolerance);
      double fPostInv = 1.0/fPostTot;
#if 0 /* SSE3 is slower */
      vector_multiply_by(fPost, fPostInv, 4);
#else
      for (j = 0; j < 4; j++)
	fPost[j] *= fPostInv;
#endif

      /* and finally, divide by stat again & rotate to give the new frequencies */
      matrixt_by_vector4(transmat->eigeninvT, fPost, /*OUT*/fOut);
    }  /* end loop over position i */
  } else if (nCodes == 20) {	/* matrix model on amino acids */
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];
    numeric_t f1mix[20] ALIGNED;
    numeric_t f2mix[20] ALIGNED;
    
    for (i=0; i < nPos; i++) {
      if (p1->codes[i] == NOCODE && p2->codes[i] == NOCODE
	  && p1->weights[i] == 0 && p2->weights[i] == 0) {
	/* aligning gap with gap -- just output a gap
	   out->codes[i] is already set to NOCODE so need not set that */
	out->weights[i] = 0;
	continue;
      }
      int iRate = rates->ratecat[i];
      numeric_t *expeigen1 = &expeigenRates1[iRate*20];
      numeric_t *expeigen2 = &expeigenRates2[iRate*20];
      numeric_t *f1 = GET_FREQ(p1,i,/*IN/OUT*/iFreq1);
      numeric_t *f2 = GET_FREQ(p2,i,/*IN/OUT*/iFreq2);
      numeric_t *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      assert(fOut != NULL);

      if (f1 == NULL) {
	f1 = &transmat->codeFreq[p1->codes[i]][0]; /* codeFreq includes an entry for NOCODE */
	double w = p1->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 20; j++)
	    f1mix[j] = w * f1[j] + (1.0-w) * fGap[j];
	  f1 = f1mix;
	}
      }
      if (f2 == NULL) {
	f2 = &transmat->codeFreq[p2->codes[i]][0];
	double w = p2->weights[i];
	if (w > 0.0 && w < 1.0) {
	  for (j = 0; j < 20; j++)
	    f2mix[j] = w * f2[j] + (1.0-w) * fGap[j];
	  f2 = f2mix;
	}
      }
      numeric_t fMult1[20] ALIGNED;	/* rotated1 * expeigen1 */
      numeric_t fMult2[20] ALIGNED;	/* rotated2 * expeigen2 */
      vector_multiply(f1, expeigen1, 20, /*OUT*/fMult1);
      vector_multiply(f2, expeigen2, 20, /*OUT*/fMult2);
      numeric_t fPost[20] ALIGNED;		/* in  unrotated space */
      for (j = 0; j < 20; j++) {
	numeric_t value = vector_dot_product_rot(fMult1, fMult2, &transmat->codeFreq[j][0], 20)
	  * transmat->statinv[j];
	/* Added this logic try to avoid rare numerical problems */
	fPost[j] = value >= 0 ? value : 0;
      }
      double fPostTot = vector_sum(fPost, 20);
      assert(fPostTot > fPostTotalTolerance);
      double fPostInv = 1.0/fPostTot;
      vector_multiply_by(/*IN/OUT*/fPost, fPostInv, 20);
      int ch = -1;		/* the dominant character, if any */
      if (!exactML) {
	for (j = 0; j < 20; j++) {
	  if (fPost[j] >= approxMLminf) {
	    ch = j;
	    break;
	  }
	}
      }

      /* now, see if we can use the approximation 
	 fPost ~= (1 or 0) * w + nearP * (1-w)
	 to avoid rotating */
      double w = 0;
      if (ch >= 0) {
	w = (fPost[ch] - transmat->nearP[ch][ch]) / (1.0 - transmat->nearP[ch][ch]);
	for (j = 0; j < 20; j++) {
	  if (j != ch) {
	    double fRough = (1.0-w) * transmat->nearP[ch][j];
	    if (fRough < fPost[j]  * approxMLminratio) {
	      ch = -1;		/* give up on the approximation */
	      break;
	    }
	  }
	}
      }
      if (ch >= 0) {
	nAAPosteriorRough++;
	double wInvStat = w * transmat->statinv[ch];
	for (j = 0; j < 20; j++)
	  fOut[j] = wInvStat * transmat->codeFreq[ch][j] + (1.0-w) * transmat->nearFreq[ch][j];
      } else {
	/* and finally, divide by stat again & rotate to give the new frequencies */
	nAAPosteriorExact++;
	for (j = 0; j < 20; j++)
	  fOut[j] = vector_multiply_sum(fPost, &transmat->eigeninv[j][0], 20);
      }
    } /* end loop over position i */
  } else {
    assert(0);			/* illegal nCodes */
  }

  if (transmat != NULL) {
    expeigenRates1 = myfree(expeigenRates1, sizeof(numeric_t) * rates->nRateCategories * nCodes);
    expeigenRates2 = myfree(expeigenRates2, sizeof(numeric_t) * rates->nRateCategories * nCodes);
  }

  /* Reallocate out->vectors to be the right size */
  out->nVectors = iFreqOut;
  if (out->nVectors == 0)
    out->vectors = (numeric_t*)myfree(out->vectors, sizeof(numeric_t)*nCodes*nPos);
  else
    out->vectors = (numeric_t*)myrealloc(out->vectors,
				     /*OLDSIZE*/sizeof(numeric_t)*nCodes*nPos,
				     /*NEWSIZE*/sizeof(numeric_t)*nCodes*out->nVectors,
				     /*copy*/true); /* try to save space */
  nProfileFreqAlloc += out->nVectors;
  nProfileFreqAvoid += nPos - out->nVectors;

  /* compute total constraints */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] = p1->nOn[i] + p2->nOn[i];
    out->nOff[i] = p1->nOff[i] + p2->nOff[i];
  }
  nPosteriorCompute++;
  return(out);
}

double *PSameVector(double length, rates_t *rates) {
  double *pSame = mymalloc(sizeof(double) * rates->nRateCategories);
  int iRate;
  for (iRate = 0; iRate < rates->nRateCategories; iRate++)
    pSame[iRate] = 0.25 + 0.75 * exp((-4.0/3.0) * fabs(length*rates->rates[iRate]));
  return(pSame);
}

double *PDiffVector(double *pSame, rates_t *rates) {
  double *pDiff = mymalloc(sizeof(double) * rates->nRateCategories);
  int iRate;
  for (iRate = 0; iRate < rates->nRateCategories; iRate++)
    pDiff[iRate] = (1.0 - pSame[iRate])/3.0;
  return(pDiff);
}

numeric_t *ExpEigenRates(double length, transition_matrix_t *transmat, rates_t *rates) {
  numeric_t *expeigen = mymalloc(sizeof(numeric_t) * nCodes * rates->nRateCategories);
  int iRate, j;
  for (iRate = 0; iRate < rates->nRateCategories; iRate++) {
    for (j = 0; j < nCodes; j++) {
      double relLen = length * rates->rates[iRate];
      /* very short branch lengths lead to numerical problems so prevent them */
      if (relLen < MLMinRelBranchLength)
	relLen  = MLMinRelBranchLength;
      expeigen[iRate*nCodes + j] = exp(relLen * transmat->eigenval[j]);
    }
  }
  return(expeigen);
}

double PairLogLk(profile_t *pA, profile_t *pB, double length, int nPos,
		 /*OPTIONAL*/transition_matrix_t *transmat,
		 rates_t *rates,
		 /*OPTIONAL IN/OUT*/double *site_likelihoods) {
  double lk = 1.0;
  double loglk = 0.0;		/* stores underflow of lk during the loop over positions */
  int i,j;
  assert(rates != NULL && rates->nRateCategories > 0);
  numeric_t *expeigenRates = NULL;
  if (transmat != NULL)
    expeigenRates = ExpEigenRates(length, transmat, rates);

  if (transmat == NULL) {	/* Jukes-Cantor */
    assert (nCodes == 4);
    double *pSame = PSameVector(length, rates);
    double *pDiff = PDiffVector(pSame, rates);
    
    int iFreqA = 0;
    int iFreqB = 0;
    for (i = 0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      double wA = pA->weights[i];
      double wB = pB->weights[i];
      int codeA = pA->codes[i];
      int codeB = pB->codes[i];
      numeric_t *fA = GET_FREQ(pA,i,/*IN/OUT*/iFreqA);
      numeric_t *fB = GET_FREQ(pB,i,/*IN/OUT*/iFreqB);
      double lkAB = 0;

      if (fA == NULL && fB == NULL) {
	if (codeA == NOCODE) {	/* A is all gaps */
	  /* gap to gap is sum(j) 0.25 * (0.25 * pSame + 0.75 * pDiff) = sum(i) 0.25*0.25 = 0.25
	     gap to any character gives the same result
	  */
	  lkAB = 0.25;
	} else if (codeB == NOCODE) { /* B is all gaps */
	  lkAB = 0.25;
	} else if (codeA == codeB) { /* A and B match */
	  lkAB = pSame[iRate] * wA*wB + 0.25 * (1-wA*wB);
	} else {		/* codeA != codeB */
	  lkAB = pDiff[iRate] * wA*wB + 0.25 * (1-wA*wB);
	}
      } else if (fA == NULL) {
	/* Compare codeA to profile of B */
	if (codeA == NOCODE)
	  lkAB = 0.25;
	else
	  lkAB = wA * (pDiff[iRate] + fB[codeA] * (pSame[iRate]-pDiff[iRate])) + (1.0-wA) * 0.25;
	/* because lkAB = wA * P(codeA->B) + (1-wA) * 0.25 
	   P(codeA -> B) = sum(j) P(B==j) * (j==codeA ? pSame : pDiff)
	   = sum(j) P(B==j) * pDiff + 
	   = pDiff + P(B==codeA) * (pSame-pDiff)
	*/
      } else if (fB == NULL) { /* Compare codeB to profile of A */
	if (codeB == NOCODE)
	  lkAB = 0.25;
	else
	  lkAB = wB * (pDiff[iRate] + fA[codeB] * (pSame[iRate]-pDiff[iRate])) + (1.0-wB) * 0.25;
      } else { /* both are full profiles */
	for (j = 0; j < 4; j++)
	  lkAB += fB[j] * (fA[j] * pSame[iRate] + (1-fA[j])* pDiff[iRate]); /* P(A|B) */
      }
      assert(lkAB > 0);
      lk *= lkAB;
      while (lk < LkUnderflow) {
	lk *= LkUnderflowInv;
	loglk -= LogLkUnderflow;
      }
      if (site_likelihoods != NULL)
	site_likelihoods[i] *= lkAB;
    }
    pSame = myfree(pSame, sizeof(double) * rates->nRateCategories);
    pDiff = myfree(pDiff, sizeof(double) * rates->nRateCategories);
  } else if (nCodes == 4) {	/* matrix model on nucleotides */
    int iFreqA = 0;
    int iFreqB = 0;
    numeric_t fAmix[4], fBmix[4];
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];

    for (i = 0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      numeric_t *expeigen = &expeigenRates[iRate*4];
      double wA = pA->weights[i];
      double wB = pB->weights[i];
      if (wA == 0 && wB == 0 && pA->codes[i] == NOCODE && pB->codes[i] == NOCODE) {
	/* Likelihood of A vs B is 1, so nothing changes
	   Do not need to advance iFreqA or iFreqB */
	continue;		
      }
      numeric_t *fA = GET_FREQ(pA,i,/*IN/OUT*/iFreqA);
      numeric_t *fB = GET_FREQ(pB,i,/*IN/OUT*/iFreqB);
      if (fA == NULL)
	fA = &transmat->codeFreq[pA->codes[i]][0];
      if (wA > 0.0 && wA < 1.0) {
	for (j  = 0; j < 4; j++)
	  fAmix[j] = wA*fA[j] + (1.0-wA)*fGap[j];
	fA = fAmix;
      }
      if (fB == NULL)
	fB = &transmat->codeFreq[pB->codes[i]][0];
      if (wB > 0.0 && wB < 1.0) {
	for (j  = 0; j < 4; j++)
	  fBmix[j] = wB*fB[j] + (1.0-wB)*fGap[j];
	fB = fBmix;
      }
      /* SSE3 instructions do not speed this step up:
	 numeric_t lkAB = vector_multiply3_sum(expeigen, fA, fB); */
		// dsp this is where check for <=0 was added in 2.1.1.LG
      double lkAB = 0;
      for (j = 0; j < 4; j++)
	lkAB += expeigen[j]*fA[j]*fB[j];
      assert(lkAB > 0);
      if (site_likelihoods != NULL)
	site_likelihoods[i] *= lkAB;
      lk *= lkAB;
      while (lk < LkUnderflow) {
	lk *= LkUnderflowInv;
	loglk -= LogLkUnderflow;
      }
      while (lk > LkUnderflowInv) {
	lk *= LkUnderflow;
	loglk += LogLkUnderflow;
      }
    }
  } else if (nCodes == 20) {	/* matrix model on amino acids */
    int iFreqA = 0;
    int iFreqB = 0;
    numeric_t fAmix[20], fBmix[20];
    numeric_t *fGap = &transmat->codeFreq[NOCODE][0];

    for (i = 0; i < nPos; i++) {
      int iRate = rates->ratecat[i];
      numeric_t *expeigen = &expeigenRates[iRate*20];
      double wA = pA->weights[i];
      double wB = pB->weights[i];
      if (wA == 0 && wB == 0 && pA->codes[i] == NOCODE && pB->codes[i] == NOCODE) {
	/* Likelihood of A vs B is 1, so nothing changes
	   Do not need to advance iFreqA or iFreqB */
	continue;		
      }
      numeric_t *fA = GET_FREQ(pA,i,/*IN/OUT*/iFreqA);
      numeric_t *fB = GET_FREQ(pB,i,/*IN/OUT*/iFreqB);
      if (fA == NULL)
	fA = &transmat->codeFreq[pA->codes[i]][0];
      if (wA > 0.0 && wA < 1.0) {
	for (j  = 0; j < 20; j++)
	  fAmix[j] = wA*fA[j] + (1.0-wA)*fGap[j];
	fA = fAmix;
      }
      if (fB == NULL)
	fB = &transmat->codeFreq[pB->codes[i]][0];
      if (wB > 0.0 && wB < 1.0) {
	for (j  = 0; j < 20; j++)
	  fBmix[j] = wB*fB[j] + (1.0-wB)*fGap[j];
	fB = fBmix;
      }
      numeric_t lkAB = vector_multiply3_sum(expeigen, fA, fB, 20);
      if (!(lkAB > 0)) {
	/* If this happens, it indicates a numerical problem that needs to be addressed elsewhere,
	   so report all the details */
	fprintf(stderr, "# FastTree.c::PairLogLk -- numerical problem!\n");
	fprintf(stderr, "# This block is intended for loading into R\n");

	fprintf(stderr, "lkAB = %.8g\n", lkAB);
	fprintf(stderr, "Branch_length= %.8g\nalignment_position=%d\nnCodes=%d\nrate_category=%d\nrate=%.8g\n",
		length, i, nCodes, iRate, rates->rates[iRate]);
	fprintf(stderr, "wA=%.8g\nwB=%.8g\n", wA, wB);
	fprintf(stderr, "codeA = %d\ncodeB = %d\n", pA->codes[i], pB->codes[i]);

	fprintf(stderr, "fA = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", fA[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "fB = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", fB[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "stat = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", transmat->stat[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "eigenval = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", transmat->eigenval[j]);
	fprintf(stderr,")\n");

	fprintf(stderr, "expeigen = c(");
	for (j = 0; j < nCodes; j++) fprintf(stderr, "%s %.8g", j==0?"":",", expeigen[j]);
	fprintf(stderr,")\n");

	int k;
	fprintf(stderr, "codeFreq = c(");
	for (j = 0; j < nCodes; j++) for(k = 0; k < nCodes; k++) fprintf(stderr, "%s %.8g", j==0 && k==0?"":",",
									     transmat->codeFreq[j][k]);
	fprintf(stderr,")\n");

	fprintf(stderr, "eigeninv = c(");
	for (j = 0; j < nCodes; j++) for(k = 0; k < nCodes; k++) fprintf(stderr, "%s %.8g", j==0 && k==0?"":",",
									     transmat->eigeninv[j][k]);
	fprintf(stderr,")\n");

	fprintf(stderr, "# Transform into matrices and compute un-rotated vectors for profiles A and B\n");
	fprintf(stderr, "codeFreq = matrix(codeFreq,nrow=20);\n");
	fprintf(stderr, "eigeninv = matrix(eigeninv,nrow=20);\n");
	fputs("unrotA = stat * (eigeninv %*% fA)\n", stderr);
	fputs("unrotB = stat * (eigeninv %*% fB)\n", stderr);
	fprintf(stderr,"# End of R block\n");
      }
      assert(lkAB > 0);
      if (site_likelihoods != NULL)
	site_likelihoods[i] *= lkAB;
      lk *= lkAB;
      while (lk < LkUnderflow) {
	lk *= LkUnderflowInv;
	loglk -= LogLkUnderflow;
      }
      while (lk > LkUnderflowInv) {
	lk *= LkUnderflow;
	loglk += LogLkUnderflow;
      }
    }
  } else {
    assert(0);			/* illegal nCodes */
  }
  if (transmat != NULL)
    expeigenRates = myfree(expeigenRates, sizeof(numeric_t) * rates->nRateCategories * 20);
  loglk += log(lk);
  nLkCompute++;
  return(loglk);
}

double MLQuartetLogLk(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
		      int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		      /*IN*/double branch_lengths[5],
		      /*OPTIONAL OUT*/double *site_likelihoods) {
  profile_t *pAB = PosteriorProfile(pA, pB,
				    branch_lengths[0], branch_lengths[1],
				    transmat,
				    rates,
				    nPos, /*nConstraints*/0);
  profile_t *pCD = PosteriorProfile(pC, pD,
				    branch_lengths[2], branch_lengths[3],
				    transmat,
				    rates,
				    nPos, /*nConstraints*/0);
  if (site_likelihoods != NULL) {
    int i;
    for (i = 0; i < nPos; i++)
      site_likelihoods[i] = 1.0;
  }
  /* Roughly, P(A,B,C,D) = P(A) P(B|A) P(D|C) P(AB | CD) */
  double loglk = PairLogLk(pA, pB, branch_lengths[0]+branch_lengths[1],
			   nPos, transmat, rates, /*OPTIONAL IN/OUT*/site_likelihoods)
    + PairLogLk(pC, pD, branch_lengths[2]+branch_lengths[3],
		nPos, transmat, rates, /*OPTIONAL IN/OUT*/site_likelihoods)
    + PairLogLk(pAB, pCD, branch_lengths[4],
		nPos, transmat, rates, /*OPTIONAL IN/OUT*/site_likelihoods);
  pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);
  pCD = FreeProfile(pCD, nPos, /*nConstraints*/0);
  return(loglk);
}

double PairNegLogLk(double x, void *data) {
  quartet_opt_t *qo = (quartet_opt_t *)data;
  assert(qo != NULL);
  assert(qo->pair1 != NULL && qo->pair2 != NULL);
  qo->nEval++;
  double loglk = PairLogLk(qo->pair1, qo->pair2, x, qo->nPos, qo->transmat, qo->rates, /*site_lk*/NULL);
  assert(loglk < 1e100);
  if (verbose > 5)
    fprintf(stderr, "PairLogLk(%.4f) =  %.4f\n", x, loglk);
  return(-loglk);
}

double MLQuartetOptimize(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
			 int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
			 /*IN/OUT*/double branch_lengths[5],
			 /*OPTIONAL OUT*/bool *pStarTest,
			 /*OPTIONAL OUT*/double *site_likelihoods) {
  int j;
  double start_length[5];
  for (j = 0; j < 5; j++) {
    start_length[j] = branch_lengths[j];
    if (branch_lengths[j] < MLMinBranchLength)
      branch_lengths[j] = MLMinBranchLength;
  }
  quartet_opt_t qopt = { nPos, transmat, rates, /*nEval*/0,
			 /*pair1*/NULL, /*pair2*/NULL };
  double f2x, negloglk;

  if (pStarTest != NULL)
    *pStarTest = false;

  /* First optimize internal branch, then branch to A, B, C, D, in turn
     May use star test to quit after internal branch
   */
  profile_t *pAB = PosteriorProfile(pA, pB,
				    branch_lengths[LEN_A], branch_lengths[LEN_B],
				    transmat, rates, nPos, /*nConstraints*/0);
  profile_t *pCD = PosteriorProfile(pC, pD,
				    branch_lengths[LEN_C], branch_lengths[LEN_D],
				    transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pAB;
  qopt.pair2 = pCD;
  branch_lengths[LEN_I] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_I],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);

  if (pStarTest != NULL) {
    assert(site_likelihoods == NULL);
    double loglkStar = -PairNegLogLk(MLMinBranchLength, &qopt);
    if (loglkStar < -negloglk - closeLogLkLimit) {
      *pStarTest = true;
      double off = PairLogLk(pA, pB,
			     branch_lengths[LEN_A] + branch_lengths[LEN_B],
			     qopt.nPos, qopt.transmat, qopt.rates, /*site_lk*/NULL)
	+ PairLogLk(pC, pD,
		    branch_lengths[LEN_C] + branch_lengths[LEN_D],
		    qopt.nPos, qopt.transmat, qopt.rates, /*site_lk*/NULL);
      pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);
      pCD = FreeProfile(pCD, nPos, /*nConstraints*/0);
      return (-negloglk + off);
    }
  }
  pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);
  profile_t *pBCD = PosteriorProfile(pB, pCD,
				     branch_lengths[LEN_B], branch_lengths[LEN_I],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pA;
  qopt.pair2 = pBCD;
  branch_lengths[LEN_A] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_A],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);
  pBCD = FreeProfile(pBCD, nPos, /*nConstraints*/0);
  profile_t *pACD = PosteriorProfile(pA, pCD,
				     branch_lengths[LEN_A], branch_lengths[LEN_I],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pB;
  qopt.pair2 = pACD;
  branch_lengths[LEN_B] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_B],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);
  pACD = FreeProfile(pACD, nPos, /*nConstraints*/0);
  pCD = FreeProfile(pCD, nPos, /*nConstraints*/0);
  pAB = PosteriorProfile(pA, pB,
			 branch_lengths[LEN_A], branch_lengths[LEN_B],
			 transmat, rates, nPos, /*nConstraints*/0);
  profile_t *pABD = PosteriorProfile(pAB, pD,
				     branch_lengths[LEN_I], branch_lengths[LEN_D],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pC;
  qopt.pair2 = pABD;
  branch_lengths[LEN_C] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_C],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);
  pABD = FreeProfile(pABD, nPos, /*nConstraints*/0);
  profile_t *pABC = PosteriorProfile(pAB, pC,
				     branch_lengths[LEN_I], branch_lengths[LEN_C],
				     transmat, rates, nPos, /*nConstraints*/0);
  qopt.pair1 = pD;
  qopt.pair2 = pABC;
  branch_lengths[LEN_D] = onedimenmin(/*xmin*/MLMinBranchLength,
				      /*xguess*/branch_lengths[LEN_D],
				      /*xmax*/6.0,
				      PairNegLogLk,
				      /*data*/&qopt,
				      /*ftol*/MLFTolBranchLength,
				      /*atol*/MLMinBranchLengthTolerance,
				      /*OUT*/&negloglk,
				      /*OUT*/&f2x);

  /* Compute the total quartet likelihood
     PairLogLk(ABC,D) + PairLogLk(AB,C) + PairLogLk(A,B)
   */
  double loglkABCvsD = -negloglk;
  if (site_likelihoods) {
    for (j = 0; j < nPos; j++)
      site_likelihoods[j] = 1.0;
    PairLogLk(pABC, pD, branch_lengths[LEN_D],
	      qopt.nPos, qopt.transmat, qopt.rates, /*IN/OUT*/site_likelihoods);
  }
  double quartetloglk = loglkABCvsD
    + PairLogLk(pAB, pC, branch_lengths[LEN_I] + branch_lengths[LEN_C],
		qopt.nPos, qopt.transmat, qopt.rates,
		/*IN/OUT*/site_likelihoods)
    + PairLogLk(pA, pB, branch_lengths[LEN_A] + branch_lengths[LEN_B],
		qopt.nPos, qopt.transmat, qopt.rates,
		/*IN/OUT*/site_likelihoods);

  pABC = FreeProfile(pABC, nPos, /*nConstraints*/0);
  pAB = FreeProfile(pAB, nPos, /*nConstraints*/0);

  if (verbose > 3) {
    double loglkStart = MLQuartetLogLk(pA, pB, pC, pD, nPos, transmat, rates, start_length, /*site_lk*/NULL);
    fprintf(stderr, "Optimize loglk from %.5f to %.5f eval %d lengths from\n"
	    "   %.5f %.5f %.5f %.5f %.5f to\n"
	    "   %.5f %.5f %.5f %.5f %.5f\n",
	    loglkStart, quartetloglk, qopt.nEval,
	    start_length[0], start_length[1], start_length[2], start_length[3], start_length[4],
	    branch_lengths[0], branch_lengths[1], branch_lengths[2], branch_lengths[3], branch_lengths[4]);
  }
  return(quartetloglk);
}

nni_t MLQuartetNNI(profile_t *profiles[4],
		   /*OPTIONAL*/transition_matrix_t *transmat,
		   rates_t *rates,
		   int nPos, int nConstraints,
		   /*OUT*/double criteria[3], /* The three potential quartet log-likelihoods */
		   /*IN/OUT*/numeric_t len[5],
		   bool bFast)
{
  int i;
  double lenABvsCD[5] = {len[LEN_A], len[LEN_B], len[LEN_C], len[LEN_D], len[LEN_I]};
  double lenACvsBD[5] = {len[LEN_A], len[LEN_C], len[LEN_B], len[LEN_D], len[LEN_I]};   /* Swap B & C */
  double lenADvsBC[5] = {len[LEN_A], len[LEN_D], len[LEN_C], len[LEN_B], len[LEN_I]};   /* Swap B & D */
  bool bConsiderAC = true;
  bool bConsiderAD = true;
  int iRound;
  int nRounds = mlAccuracy < 2 ? 2 : mlAccuracy;
  double penalty[3];
  QuartetConstraintPenalties(profiles, nConstraints, /*OUT*/penalty);
  if (penalty[ABvsCD] > penalty[ACvsBD] || penalty[ABvsCD] > penalty[ADvsBC])
    bFast = false;
#ifdef OPENMP
      bFast = false;		/* turn off star topology test */
#endif

  for (iRound = 0; iRound < nRounds; iRound++) {
    bool bStarTest = false;
    {
#ifdef OPENMP
      #pragma omp parallel
      #pragma omp sections
#endif
      {
#ifdef OPENMP
        #pragma omp section
#endif
	{
	  criteria[ABvsCD] = MLQuartetOptimize(profiles[0], profiles[1], profiles[2], profiles[3],
					       nPos, transmat, rates,
					       /*IN/OUT*/lenABvsCD,
					       bFast ? &bStarTest : NULL,
					       /*site_likelihoods*/NULL)
	    - penalty[ABvsCD];	/* subtract penalty b/c we are trying to maximize log lk */
	}

#ifdef OPENMP
        #pragma omp section
#else
	if (bStarTest) {
	  nStarTests++;
	  criteria[ACvsBD] = -1e20;
	  criteria[ADvsBC] = -1e20;
	  len[LEN_I] = lenABvsCD[LEN_I];
	  return(ABvsCD);
	}
#endif
	{
	  if (bConsiderAC)
	    criteria[ACvsBD] = MLQuartetOptimize(profiles[0], profiles[2], profiles[1], profiles[3],
						 nPos, transmat, rates,
						 /*IN/OUT*/lenACvsBD, NULL, /*site_likelihoods*/NULL)
	      - penalty[ACvsBD];
	}
	
#ifdef OPENMP
        #pragma omp section
#endif
	{
	  if (bConsiderAD)
	    criteria[ADvsBC] = MLQuartetOptimize(profiles[0], profiles[3], profiles[2], profiles[1],
						 nPos, transmat, rates,
						 /*IN/OUT*/lenADvsBC, NULL, /*site_likelihoods*/NULL)
	      - penalty[ADvsBC];
	}
      }
    } /* end parallel sections */
    if (mlAccuracy < 2) {
      /* If clearly worse then ABvsCD, or have short internal branch length and worse, then
         give up */
      if (criteria[ACvsBD] < criteria[ABvsCD] - closeLogLkLimit
	  || (lenACvsBD[LEN_I] <= 2.0*MLMinBranchLength && criteria[ACvsBD] < criteria[ABvsCD]))
	bConsiderAC = false;
      if (criteria[ADvsBC] < criteria[ABvsCD] - closeLogLkLimit
	  || (lenADvsBC[LEN_I] <= 2.0*MLMinBranchLength && criteria[ADvsBC] < criteria[ABvsCD]))
	bConsiderAD = false;
      if (!bConsiderAC && !bConsiderAD)
	break;
      /* If clearly better than either alternative, then give up
         (Comparison is probably biased in favor of ABvsCD anyway) */
      if (criteria[ACvsBD] > criteria[ABvsCD] + closeLogLkLimit
	  && criteria[ACvsBD] > criteria[ADvsBC] + closeLogLkLimit)
	break;
      if (criteria[ADvsBC] > criteria[ABvsCD] + closeLogLkLimit
	  && criteria[ADvsBC] > criteria[ACvsBD] + closeLogLkLimit)
	break;
    }
  } /* end loop over rounds */

  if (verbose > 2) {
    fprintf(stderr, "Optimized quartet for %d rounds: ABvsCD %.5f ACvsBD %.5f ADvsBC %.5f\n",
	    iRound, criteria[ABvsCD], criteria[ACvsBD], criteria[ADvsBC]);
  }
  if (criteria[ACvsBD] > criteria[ABvsCD] && criteria[ACvsBD] > criteria[ADvsBC]) {
    for (i = 0; i < 5; i++) len[i] = lenACvsBD[i];
    return(ACvsBD);
  } else if (criteria[ADvsBC] > criteria[ABvsCD] && criteria[ADvsBC] > criteria[ACvsBD]) {
    for (i = 0; i < 5; i++) len[i] = lenADvsBC[i];
    return(ADvsBC);
  } else {
    for (i = 0; i < 5; i++) len[i] = lenABvsCD[i];
    return(ABvsCD);
  }
}

double TreeLength(/*IN/OUT*/NJ_t *NJ, bool recomputeProfiles) {
  if (recomputeProfiles) {
    traversal_t traversal2 = InitTraversal(NJ);
    int j = NJ->root;
    while((j = TraversePostorder(j, NJ, /*IN/OUT*/traversal2, /*pUp*/NULL)) >= 0) {
      /* nothing to do for leaves or root */
      if (j >= NJ->nSeq && j != NJ->root)
	SetProfile(/*IN/OUT*/NJ, j, /*noweight*/-1.0);
    }
    traversal2 = FreeTraversal(traversal2,NJ);
  }
  UpdateBranchLengths(/*IN/OUT*/NJ);
  double total_len = 0;
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++)
    total_len += NJ->branchlength[iNode];
  return(total_len);
}

double TreeLogLk(/*IN*/NJ_t *NJ, /*OPTIONAL OUT*/double *site_loglk) {
  int i;
  if (NJ->nSeq < 2)
    return(0.0);
  double loglk = 0.0;
  double *site_likelihood = NULL;
  if (site_loglk != NULL) {
    site_likelihood = mymalloc(sizeof(double)*NJ->nPos);
    for (i = 0; i < NJ->nPos; i++) {
      site_likelihood[i] = 1.0;
      site_loglk[i] = 0.0;
    }
  }
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    int nChild = NJ->child[node].nChild;
    if (nChild == 0)
      continue;
    assert(nChild >= 2);
    int *children = NJ->child[node].child;
    double loglkchild = PairLogLk(NJ->profiles[children[0]], NJ->profiles[children[1]],
				  NJ->branchlength[children[0]]+NJ->branchlength[children[1]],
				  NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/site_likelihood);
    loglk += loglkchild;
    if (site_likelihood != NULL) {
      /* prevent underflows */
      for (i = 0; i < NJ->nPos; i++) {
	while(site_likelihood[i] < LkUnderflow) {
	  site_likelihood[i] *= LkUnderflowInv;
	  site_loglk[i] -= LogLkUnderflow;
	}
      }
    }
    if (verbose > 2)
      fprintf(stderr, "At %d: LogLk(%d:%.4f,%d:%.4f) = %.3f\n",
	      node,
	      children[0], NJ->branchlength[children[0]],
	      children[1], NJ->branchlength[children[1]],
	      loglkchild);
    if (NJ->child[node].nChild == 3) {
      assert(node == NJ->root);
      /* Infer the common parent of the 1st two to define the third... */
      profile_t *pAB = PosteriorProfile(NJ->profiles[children[0]],
					NJ->profiles[children[1]],
					NJ->branchlength[children[0]],
					NJ->branchlength[children[1]],
					NJ->transmat, &NJ->rates,
					NJ->nPos, /*nConstraints*/0);
      double loglkup = PairLogLk(pAB, NJ->profiles[children[2]],
				 NJ->branchlength[children[2]],
				 NJ->nPos, NJ->transmat, &NJ->rates,
				 /*IN/OUT*/site_likelihood);
      loglk += loglkup;
      if (verbose > 2)
	fprintf(stderr, "At root %d: LogLk((%d/%d),%d:%.3f) = %.3f\n",
		node, children[0], children[1], children[2],
		NJ->branchlength[children[2]],
		loglkup);
      pAB = FreeProfile(pAB, NJ->nPos, NJ->nConstraints);
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  if (site_likelihood != NULL) {
    for (i = 0; i < NJ->nPos; i++) {
      site_loglk[i] += log(site_likelihood[i]);
    }
    site_likelihood = myfree(site_likelihood, sizeof(double)*NJ->nPos);
  }

  /* For Jukes-Cantor, with a tree of size 4, if the children of the root are
     (A,B), C, and D, then
     P(ABCD) = P(A) P(B|A) P(C|AB) P(D|ABC)
     
     Above we compute P(B|A) P(C|AB) P(D|ABC) -- note P(B|A) is at the child of root
     and P(C|AB) P(D|ABC) is at root.

     Similarly if the children of the root are C, D, and (A,B), then
     P(ABCD) = P(C|D) P(A|B) P(AB|CD) P(D), and above we compute that except for P(D)

     So we need to multiply by P(A) = 0.25, so we pay log(4) at each position
     (if ungapped). Each gapped position in any sequence reduces the payment by log(4)

     For JTT or GTR, we are computing P(A & B) and the posterior profiles are scaled to take
     the prior into account, so we do not need any correction.
     codeFreq[NOCODE] is scaled x higher so that P(-) = 1 not P(-)=1/nCodes, so gaps
     do not need to be corrected either.
   */

  if (nCodes == 4 && NJ->transmat == NULL) {
    int nGaps = 0;
    double logNCodes = log((double)nCodes);
    for (i = 0; i < NJ->nPos; i++) {
      int nGapsThisPos = 0;
      for (node = 0; node < NJ->nSeq; node++) {
	unsigned char *codes = NJ->profiles[node]->codes;
	if (codes[i] == NOCODE)
	  nGapsThisPos++;
      }
      nGaps += nGapsThisPos;
      if (site_loglk != NULL) {
	site_loglk[i] += nGapsThisPos * logNCodes;
	if (nCodes == 4 && NJ->transmat == NULL)
	  site_loglk[i] -= logNCodes;
      }
    }
    loglk -= NJ->nPos * logNCodes;
    loglk += nGaps * logNCodes;	/* do not pay for gaps -- only Jukes-Cantor */
  }
  return(loglk);
}

void SetMLGtr(/*IN/OUT*/NJ_t *NJ, /*OPTIONAL IN*/double *freq_in, /*OPTIONAL WRITE*/FILE *fpLog) {
  int i;
  assert(nCodes==4);
  gtr_opt_t gtr;
  gtr.NJ = NJ;
  gtr.fpLog = fpLog;
  if (freq_in != NULL) {
    for (i=0; i<4; i++)
      gtr.freq[i]=freq_in[i];
  } else {
    /* n[] and sum were int in FastTree 2.1.9 and earlier -- this
       caused gtr analyses to fail on analyses with >2e9 positions */
    long n[4] = {1,1,1,1};	/* pseudocounts */
    for (i=0; i<NJ->nSeq; i++) {
      unsigned char *codes = NJ->profiles[i]->codes;
      int iPos;
      for (iPos=0; iPos<NJ->nPos; iPos++)
	if (codes[iPos] < 4)
	  n[codes[iPos]]++;
    }
    long sum = n[0]+n[1]+n[2]+n[3];
    for (i=0; i<4; i++)
      gtr.freq[i] = n[i]/(double)sum;
  }
  for (i=0; i<6; i++)
    gtr.rates[i] = 1.0;
  int nRounds = mlAccuracy < 2 ? 2 : mlAccuracy;
  for (i = 0; i < nRounds; i++) {
    for (gtr.iRate = 0; gtr.iRate < 6; gtr.iRate++) {
      ProgressReport("Optimizing GTR model, step %d of %d", i*6+gtr.iRate+1, 12, 0, 0);
      double negloglk, f2x;
      gtr.rates[gtr.iRate] = onedimenmin(/*xmin*/0.05,
					 /*xguess*/gtr.rates[gtr.iRate],
					 /*xmax*/20.0,
					 GTRNegLogLk,
					 /*data*/&gtr,
					 /*ftol*/0.001,
					 /*atol*/0.0001,
					 /*OUT*/&negloglk,
					 /*OUT*/&f2x);
    }
  }
  /* normalize gtr so last rate is 1 -- specifying that rate separately is useful for optimization only */
  for (i = 0; i < 5; i++)
    gtr.rates[i] /= gtr.rates[5];
  gtr.rates[5] = 1.0;
  if (verbose) {
    fprintf(stderr, "GTR Frequencies: %.4f %.4f %.4f %.4f\n", gtr.freq[0], gtr.freq[1], gtr.freq[2], gtr.freq[3]);
    fprintf(stderr, "GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f\n",
	    gtr.rates[0],gtr.rates[1],gtr.rates[2],gtr.rates[3],gtr.rates[4],gtr.rates[5]);
  }
  if (fpLog != NULL) {
    fprintf(fpLog, "GTRFreq\t%.4f\t%.4f\t%.4f\t%.4f\n", gtr.freq[0], gtr.freq[1], gtr.freq[2], gtr.freq[3]);
    fprintf(fpLog, "GTRRates\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	    gtr.rates[0],gtr.rates[1],gtr.rates[2],gtr.rates[3],gtr.rates[4],gtr.rates[5]);
  }
  myfree(NJ->transmat, sizeof(transition_matrix_t));
  NJ->transmat = CreateGTR(gtr.rates, gtr.freq);
  RecomputeMLProfiles(/*IN/OUT*/NJ);
  OptimizeAllBranchLengths(/*IN/OUT*/NJ);
}

double GTRNegLogLk(double x, void *data) {
  
  gtr_opt_t *gtr = (gtr_opt_t*)data;
  assert(nCodes == 4);
  assert(gtr->NJ != NULL);
  assert(gtr->iRate >= 0 && gtr->iRate < 6);
  assert(x > 0);
  transition_matrix_t *old = gtr->NJ->transmat;
  double rates[6];
  int i;
  for (i = 0; i < 6; i++)
    rates[i] = gtr->rates[i];
  rates[gtr->iRate] = x;

  FILE *fpLog = gtr->fpLog;
  if (fpLog)
    fprintf(fpLog, "GTR_Opt\tfreq %.5f %.5f %.5f %.5f rates %.5f %.5f %.5f %.5f %.5f %.5f\n",
          gtr->freq[0], gtr->freq[1], gtr->freq[2], gtr->freq[3],
          rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]);

  gtr->NJ->transmat = CreateGTR(rates, gtr->freq);
  RecomputeMLProfiles(/*IN/OUT*/gtr->NJ);
  double loglk = TreeLogLk(gtr->NJ, /*site_loglk*/NULL);
  myfree(gtr->NJ->transmat, sizeof(transition_matrix_t));
  gtr->NJ->transmat = old;
  /* Do not recompute profiles -- assume the caller will do that */
  if (verbose > 2)
    fprintf(stderr, "GTR LogLk(%.5f %.5f %.5f %.5f %.5f %.5f) = %f\n",
	    rates[0], rates[1], rates[2], rates[3], rates[4], rates[5], loglk);
  if (fpLog)
    fprintf(fpLog, "GTR_Opt\tGTR LogLk(%.5f %.5f %.5f %.5f %.5f %.5f) = %f\n",
	    rates[0], rates[1], rates[2], rates[3], rates[4], rates[5], loglk);
  return(-loglk);
}

/* Caller must free the resulting vector of n rates */
numeric_t *MLSiteRates(int nRateCategories) {
  /* Even spacing from 1/nRate to nRate */
  double logNCat = log((double)nRateCategories);
  double logMinRate = -logNCat;
  double logMaxRate = logNCat;
  double logd = (logMaxRate-logMinRate)/(double)(nRateCategories-1);

  numeric_t *rates = mymalloc(sizeof(numeric_t)*nRateCategories);
  int i;
  for (i = 0; i < nRateCategories; i++)
    rates[i] = exp(logMinRate + logd*(double)i);
  return(rates);
}

double *MLSiteLikelihoodsByRate(/*IN*/NJ_t *NJ, /*IN*/numeric_t *rates, int nRateCategories) {
  double *site_loglk = mymalloc(sizeof(double)*NJ->nPos*nRateCategories);

  /* save the original rates */
  assert(NJ->rates.nRateCategories > 0);
  numeric_t *oldRates = NJ->rates.rates;
  NJ->rates.rates = mymalloc(sizeof(numeric_t) * NJ->rates.nRateCategories);

  /* Compute site likelihood for each rate */
  int iPos;
  int iRate;
  for (iRate = 0; iRate  < nRateCategories; iRate++) {
    int i;
    for (i = 0; i < NJ->rates.nRateCategories; i++)
      NJ->rates.rates[i] = rates[iRate];
    RecomputeMLProfiles(/*IN/OUT*/NJ);
    double loglk = TreeLogLk(NJ, /*OUT*/&site_loglk[NJ->nPos*iRate]);
    ProgressReport("Site likelihoods with rate category %d of %d", iRate+1, nRateCategories, 0, 0);
    if(verbose > 2) {
      fprintf(stderr, "Rate %.3f Loglk %.3f SiteLogLk", rates[iRate], loglk);
      for (iPos = 0; iPos < NJ->nPos; iPos++)
	fprintf(stderr,"\t%.3f", site_loglk[NJ->nPos*iRate + iPos]);
      fprintf(stderr,"\n");
    }
  }

  /* restore original rates and profiles */
  myfree(NJ->rates.rates, sizeof(numeric_t) * NJ->rates.nRateCategories);
  NJ->rates.rates = oldRates;
  RecomputeMLProfiles(/*IN/OUT*/NJ);

  return(site_loglk);
}

void SetMLRates(/*IN/OUT*/NJ_t *NJ, int nRateCategories) {
  assert(nRateCategories > 0);
  AllocRateCategories(/*IN/OUT*/&NJ->rates, 1, NJ->nPos); /* set to 1 category of rate 1 */
  if (nRateCategories == 1) {
    RecomputeMLProfiles(/*IN/OUT*/NJ);
    return;
  }
  numeric_t *rates = MLSiteRates(nRateCategories);
  double *site_loglk = MLSiteLikelihoodsByRate(/*IN*/NJ, /*IN*/rates, nRateCategories);

  /* Select best rate for each site, correcting for the prior
     For a prior, use a gamma distribution with shape parameter 3, scale 1/3, so
     Prior(rate) ~ rate**2 * exp(-3*rate)
     log Prior(rate) = C + 2 * log(rate) - 3 * rate
  */
  double sumRates = 0;
  int iPos;
  int iRate;
  for (iPos = 0; iPos < NJ->nPos; iPos++) {
    int iBest = -1;
    double dBest = -1e20;
    for (iRate = 0; iRate < nRateCategories; iRate++) {
      double site_loglk_with_prior = site_loglk[NJ->nPos*iRate + iPos]
	+ 2.0 * log(rates[iRate]) - 3.0 * rates[iRate];
      if (site_loglk_with_prior > dBest) {
	iBest = iRate;
	dBest = site_loglk_with_prior;
      }
    }
    if (verbose > 2)
      fprintf(stderr, "Selected rate category %d rate %.3f for position %d\n",
	      iBest, rates[iBest], iPos+1);
    NJ->rates.ratecat[iPos] = iBest;
    sumRates += rates[iBest];
  }
  site_loglk = myfree(site_loglk, sizeof(double)*NJ->nPos*nRateCategories);

  /* Force the rates to average to 1 */
  double avgRate = sumRates/NJ->nPos;
  for (iRate = 0; iRate < nRateCategories; iRate++)
    rates[iRate] /= avgRate;
  
  /* Save the rates */
  NJ->rates.rates = myfree(NJ->rates.rates, sizeof(numeric_t) * NJ->rates.nRateCategories);
  NJ->rates.rates = rates;
  NJ->rates.nRateCategories = nRateCategories;

  /* Update profiles based on rates */
  RecomputeMLProfiles(/*IN/OUT*/NJ);

  if (verbose) {
    fprintf(stderr, "Switched to using %d rate categories (CAT approximation)\n", nRateCategories);
    fprintf(stderr, "Rate categories were divided by %.3f so that average rate = 1.0\n", avgRate);
    fprintf(stderr, "CAT-based log-likelihoods may not be comparable across runs\n");
    if (!gammaLogLk)
      fprintf(stderr, "Use -gamma for approximate but comparable Gamma(20) log-likelihoods\n");
  }
}

double GammaLogLk(/*IN*/siteratelk_t *s, /*OPTIONAL OUT*/double *gamma_loglk_sites) {
  int iRate, iPos;
  double *dRate = mymalloc(sizeof(double) * s->nRateCats);
  for (iRate = 0; iRate < s->nRateCats; iRate++) {
    /* The probability density for each rate is approximated by the total
       density between the midpoints */
    double pMin = iRate == 0 ? 0.0 :
      PGamma(s->mult * (s->rates[iRate-1] + s->rates[iRate])/2.0, s->alpha);
    double pMax = iRate == s->nRateCats-1 ? 1.0 :
      PGamma(s->mult * (s->rates[iRate]+s->rates[iRate+1])/2.0, s->alpha);
    dRate[iRate] = pMax-pMin;
  }

  double loglk = 0.0;
  for (iPos = 0; iPos < s->nPos; iPos++) {
    /* Prevent underflow on large trees by comparing to maximum loglk */
    double maxloglk = -1e20;
    for (iRate = 0; iRate < s->nRateCats; iRate++) {
      double site_loglk = s->site_loglk[s->nPos*iRate + iPos];
      if (site_loglk > maxloglk)
	maxloglk = site_loglk;
    }
    double rellk = 0; /* likelihood scaled by exp(maxloglk) */
    for (iRate = 0; iRate < s->nRateCats; iRate++) {
      double lk = exp(s->site_loglk[s->nPos*iRate + iPos] - maxloglk);
      rellk += lk * dRate[iRate];
    }
    double loglk_site = maxloglk + log(rellk);
    loglk += loglk_site;
    if (gamma_loglk_sites != NULL)
      gamma_loglk_sites[iPos] = loglk_site;
  }
  dRate = myfree(dRate, sizeof(double)*s->nRateCats);
  return(loglk);
}

double OptAlpha(double alpha, void *data) {
  siteratelk_t *s = (siteratelk_t *)data;
  s->alpha = alpha;
  return(-GammaLogLk(s, NULL));
}

double OptMult(double mult, void *data) {
  siteratelk_t *s = (siteratelk_t *)data;
  s->mult = mult;
  return(-GammaLogLk(s, NULL));
}

/* Input site_loglk must be for each rate */
double RescaleGammaLogLk(int nPos, int nRateCats, /*IN*/numeric_t *rates, /*IN*/double *site_loglk,
			 /*OPTIONAL*/FILE *fpLog) {
  siteratelk_t s = { /*mult*/1.0, /*alpha*/1.0, nPos, nRateCats, rates, site_loglk };
  double fx, f2x;
  int i;
  fx = -GammaLogLk(&s, NULL);
  if (verbose>2)
    fprintf(stderr, "Optimizing alpha, starting at loglk %.3f\n", -fx);
  for (i = 0; i < 10; i++) {
    ProgressReport("Optimizing alpha round %d", i+1, 0, 0, 0);
    double start = fx;
    s.alpha = onedimenmin(0.01, s.alpha, 10.0, OptAlpha, &s, 0.001, 0.001, &fx, &f2x);
    if (verbose>2)
      fprintf(stderr, "Optimize alpha round %d to %.3f lk %.3f\n", i+1, s.alpha, -fx);
    s.mult = onedimenmin(0.01, s.mult, 10.0, OptMult, &s, 0.001, 0.001, &fx, &f2x);
    if (verbose>2)
      fprintf(stderr, "Optimize mult round %d to %.3f lk %.3f\n", i+1, s.mult, -fx);
    if (fx > start - 0.001) {
      if (verbose>2)
	fprintf(stderr, "Optimizing alpha & mult converged\n");
      break;
    }
  }

  double *gamma_loglk_sites = mymalloc(sizeof(double) * nPos);
  double gammaLogLk = GammaLogLk(&s, /*OUT*/gamma_loglk_sites);
  if (verbose > 0)
    fprintf(stderr, "Gamma(%d) LogLk = %.3f alpha = %.3f rescaling lengths by %.3f\n",
	    nRateCats, gammaLogLk, s.alpha, 1/s.mult);
  if (fpLog) {
    int iPos;
    int iRate;
    fprintf(fpLog, "Gamma%dLogLk\t%.3f\tApproximate\tAlpha\t%.3f\tRescale\t%.3f\n",
	    nRateCats, gammaLogLk, s.alpha, 1/s.mult);
    fprintf(fpLog, "Gamma%d\tSite\tLogLk", nRateCats);
    for (iRate = 0; iRate < nRateCats; iRate++)
      fprintf(fpLog, "\tr=%.3f", rates[iRate]/s.mult);
    fprintf(fpLog,"\n");
    for (iPos = 0; iPos < nPos; iPos++) {
      fprintf(fpLog, "Gamma%d\t%d\t%.3f", nRateCats, iPos, gamma_loglk_sites[iPos]);
      for (iRate = 0; iRate < nRateCats; iRate++)
	fprintf(fpLog, "\t%.3f", site_loglk[nPos*iRate + iPos]);
      fprintf(fpLog,"\n");
    }
  }
  gamma_loglk_sites = myfree(gamma_loglk_sites, sizeof(double) * nPos);
  return(1.0/s.mult);
}

double MLPairOptimize(profile_t *pA, profile_t *pB,
		      int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		      /*IN/OUT*/double *branch_length) {
  quartet_opt_t qopt = { nPos, transmat, rates,
			 /*nEval*/0, /*pair1*/pA, /*pair2*/pB };
  double f2x,negloglk;
  *branch_length = onedimenmin(/*xmin*/MLMinBranchLength,
			       /*xguess*/*branch_length,
			       /*xmax*/6.0,
			       PairNegLogLk,
			       /*data*/&qopt,
			       /*ftol*/MLFTolBranchLength,
			       /*atol*/MLMinBranchLengthTolerance,
			       /*OUT*/&negloglk,
			       /*OUT*/&f2x);
  return(-negloglk);		/* the log likelihood */
}

void OptimizeAllBranchLengths(/*IN/OUT*/NJ_t *NJ) {
  if (NJ->nSeq < 2)
    return;
  if (NJ->nSeq == 2) {
    int parent = NJ->root;
    assert(NJ->child[parent].nChild==2);
    int nodes[2] = { NJ->child[parent].child[0], NJ->child[parent].child[1] };
    double length = 1.0;
    (void)MLPairOptimize(NJ->profiles[nodes[0]], NJ->profiles[nodes[1]],
			 NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/&length);
    NJ->branchlength[nodes[0]] = length/2.0;
    NJ->branchlength[nodes[1]] = length/2.0;
    return;
  };

  traversal_t traversal = InitTraversal(NJ);
  profile_t **upProfiles = UpProfiles(NJ);
  int node = NJ->root;
  int iDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    int nChild = NJ->child[node].nChild;
    if (nChild > 0) {
      if ((iDone % 100) == 0)
	ProgressReport("ML Lengths %d of %d splits", iDone+1, NJ->maxnode - NJ->nSeq, 0, 0);
      iDone++;

      /* optimize the branch lengths between self, parent, and children,
         with two iterations
      */
      assert(nChild == 2 || nChild == 3);
      int nodes[3] = { NJ->child[node].child[0],
		       NJ->child[node].child[1],
		       nChild == 3 ? NJ->child[node].child[2] : node };
      profile_t *profiles[3] = { NJ->profiles[nodes[0]],
			   NJ->profiles[nodes[1]], 
			   nChild == 3 ? NJ->profiles[nodes[2]]
			   : GetUpProfile(/*IN/OUT*/upProfiles, NJ, node, /*useML*/true) };
      int iter;
      for (iter = 0; iter < 2; iter++) {
	int i;
	for (i = 0; i < 3; i++) {
	  profile_t *pA = profiles[i];
	  int b1 = (i+1) % 3;
	  int b2 = (i+2) % 3;
	  profile_t *pB = PosteriorProfile(profiles[b1], profiles[b2],
					   NJ->branchlength[nodes[b1]],
					   NJ->branchlength[nodes[b2]],
					   NJ->transmat, &NJ->rates, NJ->nPos, /*nConstraints*/0);
	  double len = NJ->branchlength[nodes[i]];
	  if (len < MLMinBranchLength)
	    len = MLMinBranchLength;
	  (void)MLPairOptimize(pA, pB, NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/&len);
	  NJ->branchlength[nodes[i]] = len;
	  pB = FreeProfile(pB, NJ->nPos, /*nConstraints*/0);
	  if (verbose>3)
	    fprintf(stderr, "Optimize length for %d to %.3f\n",
		    nodes[i], NJ->branchlength[nodes[i]]);
	}
      }
      if (node != NJ->root) {
	RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node, /*useML*/true);
	DeleteUpProfile(upProfiles, NJ, node);
      }
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}

void RecomputeMLProfiles(/*IN/OUT*/NJ_t *NJ) {
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (NJ->child[node].nChild == 2) {
      NJ->profiles[node] = FreeProfile(NJ->profiles[node], NJ->nPos, NJ->nConstraints);
      int *children = NJ->child[node].child;
      NJ->profiles[node] = PosteriorProfile(NJ->profiles[children[0]], NJ->profiles[children[1]],
					    NJ->branchlength[children[0]], NJ->branchlength[children[1]],
					    NJ->transmat, &NJ->rates, NJ->nPos, NJ->nConstraints);
    }
  }
  traversal = FreeTraversal(traversal, NJ);
}

void RecomputeProfiles(/*IN/OUT*/NJ_t *NJ, /*OPTIONAL*/distance_matrix_t *dmat) {
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (NJ->child[node].nChild == 2) {
      int *child = NJ->child[node].child;
      NJ->profiles[node] = FreeProfile(NJ->profiles[node], NJ->nPos, NJ->nConstraints);
      NJ->profiles[node] = AverageProfile(NJ->profiles[child[0]], NJ->profiles[child[1]],
					  NJ->nPos, NJ->nConstraints,
					  dmat, /*unweighted*/-1.0);
    }
  }
  traversal = FreeTraversal(traversal,NJ);
}

int NNI(/*IN/OUT*/NJ_t *NJ, int iRound, int nRounds, bool useML,
	/*IN/OUT*/nni_stats_t *stats,
	/*OUT*/double *dMaxDelta) {
  /* For each non-root node N, with children A,B, sibling C, and uncle D,
     we compare the current topology AB|CD to the alternate topologies
     AC|BD and AD|BC, by using the 4 relevant profiles.

     If useML is true, it uses quartet maximum likelihood, and it
     updates branch lengths as it goes.

     If useML is false, it uses the minimum-evolution criterion with
     log-corrected distances on profiles.  (If logdist is false, then
     the log correction is not done.) If useML is false, then NNI()
     does NOT modify the branch lengths.

     Regardless of whether it changes the topology, it recomputes the
     profile for the node, using the pairwise distances and BIONJ-like
     weightings (if bionj is set). The parent's profile has changed,
     but recomputing it is not necessary because we will visit it
     before we need it (we use postorder, so we may visit the sibling
     and its children before we visit the parent, but we never
     consider an ancestor's profile, so that is OK). When we change
     the parent's profile, this alters the uncle's up-profile, so we
     remove that.  Finally, if the topology has changed, we remove the
     up-profiles of the nodes.

     If we do an NNI during post-order traversal, the result is a bit
     tricky. E.g. if we are at node N, and have visited its children A
     and B but not its uncle C, and we do an NNI that swaps B & C,
     then the post-order traversal will visit C, and its children, but
     then on the way back up, it will skip N, as it has already
     visited it.  So, the profile of N will not be recomputed: any
     changes beneath C will not be reflected in the profile of N, and
     the profile of N will be slightly stale. This will be corrected
     on the next round of NNIs.
  */
  double supportThreshold = useML ? treeLogLkDelta : MEMinDelta;
  int i;
  *dMaxDelta = 0.0;
  int nNNIThisRound = 0;

  if (NJ->nSeq <= 3)
    return(0);			/* nothing to do */
  if (verbose > 2) {
    fprintf(stderr, "Beginning round %d of NNIs with ml? %d\n", iRound, useML?1:0);
    PrintNJInternal(/*WRITE*/stderr, NJ, /*useLen*/useML && iRound > 0 ? 1 : 0);
  }
  /* For each node the upProfile or NULL */
  profile_t **upProfiles = UpProfiles(NJ);

  traversal_t traversal = InitTraversal(NJ);

  /* Identify nodes we can skip traversing into */
  int node;
  if (fastNNI) {
    for (node = 0; node < NJ->maxnode; node++) {
      if (node != NJ->root
	  && node >= NJ->nSeq
	  && stats[node].age >= 2
	  && stats[node].subtreeAge >= 2
	  && stats[node].support > supportThreshold) {
	int nodeABCD[4];
	SetupABCD(NJ, node, NULL, NULL, /*OUT*/nodeABCD, useML);
	for (i = 0; i < 4; i++)
	  if (stats[nodeABCD[i]].age == 0 && stats[nodeABCD[i]].support > supportThreshold)
	    break;
	if (i == 4) {
	  SkipTraversalInto(node, /*IN/OUT*/traversal);
	  if (verbose > 2)
	    fprintf(stderr, "Skipping subtree at %d: child %d %d parent %d age %d subtreeAge %d support %.3f\n",
		    node, nodeABCD[0], nodeABCD[1], NJ->parent[node],
		    stats[node].age, stats[node].subtreeAge, stats[node].support);
	}
      }
    }
  }

  int iDone = 0;
  bool bUp;
  node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, &bUp)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */
    if (bUp) {
      if(verbose > 2)
	fprintf(stderr, "Going up back to node %d\n", node);
      /* No longer needed */
      for (i = 0; i < NJ->child[node].nChild; i++)
	DeleteUpProfile(upProfiles, NJ, NJ->child[node].child[i]);
      DeleteUpProfile(upProfiles, NJ, node);
      RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node, useML);
      continue;
    }
    if ((iDone % 100) == 0) {
      char buf[100];
      sprintf(buf, "%s NNI round %%d of %%d, %%d of %%d splits", useML ? "ML" : "ME");
      if (iDone > 0)
	sprintf(buf+strlen(buf), ", %d changes", nNNIThisRound);
      if (nNNIThisRound > 0)
	sprintf(buf+strlen(buf), " (max delta %.3f)", *dMaxDelta);
      ProgressReport(buf, iRound+1, nRounds, iDone+1, NJ->maxnode - NJ->nSeq);
    }
    iDone++;

    profile_t *profiles[4];
    int nodeABCD[4];
    /* Note -- during the first round of ML NNIs, we use the min-evo-based branch lengths,
       which may be suboptimal */
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, useML);

    /* Given our 4 profiles, consider doing a swap */
    int nodeA = nodeABCD[0];
    int nodeB = nodeABCD[1];
    int nodeC = nodeABCD[2];
    int nodeD = nodeABCD[3];

    nni_t choice = ABvsCD;

    if (verbose > 2)
      fprintf(stderr,"Considering NNI around %d: Swap A=%d B=%d C=%d D=up(%d) or parent %d\n",
	      node, nodeA, nodeB, nodeC, nodeD, NJ->parent[node]);
    if (verbose > 3 && useML) {
      double len[5] = { NJ->branchlength[nodeA], NJ->branchlength[nodeB], NJ->branchlength[nodeC], NJ->branchlength[nodeD],
			NJ->branchlength[node] };
      for (i=0; i < 5; i++)
	if (len[i] < MLMinBranchLength)
	  len[i] = MLMinBranchLength;
      fprintf(stderr, "Starting quartet likelihood %.3f len %.3f %.3f %.3f %.3f %.3f\n",
	      MLQuartetLogLk(profiles[0],profiles[1],profiles[2],profiles[3],NJ->nPos,NJ->transmat,&NJ->rates,len, /*site_lk*/NULL),
	      len[0], len[1], len[2], len[3], len[4]);
    }

    numeric_t newlength[5];
    double criteria[3];
    if (useML) {
      for (i = 0; i < 4; i++)
	newlength[i] = NJ->branchlength[nodeABCD[i]];
      newlength[4] = NJ->branchlength[node];
      bool bFast = mlAccuracy < 2 && stats[node].age > 0;
      choice = MLQuartetNNI(profiles, NJ->transmat, &NJ->rates, NJ->nPos, NJ->nConstraints,
			    /*OUT*/criteria, /*IN/OUT*/newlength, bFast);
    } else {
      choice = ChooseNNI(profiles, NJ->distance_matrix, NJ->nPos, NJ->nConstraints,
			 /*OUT*/criteria);
      /* invert criteria so that higher is better, as in ML case, to simplify code below */
      for (i = 0; i < 3; i++)
	criteria[i] = -criteria[i];
    }
    
    if (choice == ACvsBD) {
      /* swap B and C */
      ReplaceChild(/*IN/OUT*/NJ, node, nodeB, nodeC);
      ReplaceChild(/*IN/OUT*/NJ, NJ->parent[node], nodeC, nodeB);
    } else if (choice == ADvsBC) {
      /* swap A and C */
      ReplaceChild(/*IN/OUT*/NJ, node, nodeA, nodeC);
      ReplaceChild(/*IN/OUT*/NJ, NJ->parent[node], nodeC, nodeA);
    }
    
    if (useML) {
      /* update branch length for the internal branch, and of any
	 branches that lead to leaves, b/c those will not are not
	 the internal branch for NNI and would not otherwise be set.
      */
      if (choice == ADvsBC) {
	/* For ADvsBC, MLQuartetNNI swaps B with D, but we swap A with C */
	double length2[5] = { newlength[LEN_C], newlength[LEN_D],
			      newlength[LEN_A], newlength[LEN_B],
			      newlength[LEN_I] };
	int i;
	for (i = 0; i < 5; i++) newlength[i] = length2[i];
	/* and swap A and C */
	double tmp = newlength[LEN_A];
	newlength[LEN_A] = newlength[LEN_C];
	newlength[LEN_C] = tmp;
      } else if (choice == ACvsBD) {
	/* swap B and C */
	double tmp = newlength[LEN_B];
	newlength[LEN_B] = newlength[LEN_C];
	newlength[LEN_C] = tmp;
      }
      
      NJ->branchlength[node] = newlength[LEN_I];
      NJ->branchlength[nodeA] = newlength[LEN_A];
      NJ->branchlength[nodeB] = newlength[LEN_B];
      NJ->branchlength[nodeC] = newlength[LEN_C];
      NJ->branchlength[nodeD] = newlength[LEN_D];
    }
    
    if (verbose>2 && (choice != ABvsCD || verbose > 2))
      fprintf(stderr,"NNI around %d: Swap A=%d B=%d C=%d D=out(C) -- choose %s %s %.4f\n",
	      node, nodeA, nodeB, nodeC,
	      choice == ACvsBD ? "AC|BD" : (choice == ABvsCD ? "AB|CD" : "AD|BC"),
	      useML ? "delta-loglk" : "-deltaLen",
	      criteria[choice] - criteria[ABvsCD]);
    if(verbose >= 3 && slow && useML)
      fprintf(stderr, "Old tree lk -- %.4f\n", TreeLogLk(NJ, /*site_likelihoods*/NULL));
    
    /* update stats, *dMaxDelta, etc. */
    if (choice == ABvsCD) {
      stats[node].age++;
    } else {
      if (useML)
	nML_NNI++;
      else
	nNNI++;
      nNNIThisRound++;
      stats[node].age = 0;
      stats[nodeA].age = 0;
      stats[nodeB].age = 0;
      stats[nodeC].age = 0;
      stats[nodeD].age = 0;
    }
    stats[node].delta = criteria[choice] - criteria[ABvsCD]; /* 0 if ABvsCD */
    if (stats[node].delta > *dMaxDelta)
      *dMaxDelta = stats[node].delta;
    
    /* support is improvement of score for self over better of alternatives */
    stats[node].support = 1e20;
    for (i = 0; i < 3; i++)
      if (choice != i && criteria[choice]-criteria[i] < stats[node].support)
	stats[node].support = criteria[choice]-criteria[i];
    
    /* subtreeAge is the number of rounds since self or descendent had a significant improvement */
    if (stats[node].delta > supportThreshold)
      stats[node].subtreeAge = 0;
    else {
      stats[node].subtreeAge++;
      for (i = 0; i < 2; i++) {
	int child = NJ->child[node].child[i];
	if (stats[node].subtreeAge > stats[child].subtreeAge)
	  stats[node].subtreeAge = stats[child].subtreeAge;
      }
    }

    /* update profiles and free up unneeded up-profiles */
    if (choice == ABvsCD) {
      /* No longer needed */
      DeleteUpProfile(upProfiles, NJ, nodeA);
      DeleteUpProfile(upProfiles, NJ, nodeB);
      DeleteUpProfile(upProfiles, NJ, nodeC);
      RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node, useML);
      if(slow && useML)
	UpdateForNNI(NJ, node, upProfiles, useML);
    } else {
      UpdateForNNI(NJ, node, upProfiles, useML);
    }
    if(verbose > 2 && slow && useML) {
      /* Note we recomputed profiles back up to root already if slow */
      PrintNJInternal(/*WRITE*/stderr, NJ, /*useLen*/true);
      fprintf(stderr, "New tree lk -- %.4f\n", TreeLogLk(NJ, /*site_likelihoods*/NULL));
    }
  } /* end postorder traversal */
  traversal = FreeTraversal(traversal,NJ);
  if (verbose>=2) {
    int nUp = 0;
    for (i = 0; i < NJ->maxnodes; i++)
      if (upProfiles[i] != NULL)
	nUp++;
    fprintf(stderr, "N up profiles at end of NNI:  %d\n", nUp);
  }
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  return(nNNIThisRound);
}

nni_stats_t *InitNNIStats(NJ_t *NJ) {
  nni_stats_t *stats = mymalloc(sizeof(nni_stats_t)*NJ->maxnode);
  const int LargeAge = 1000000;
  int i;
  for (i = 0; i < NJ->maxnode; i++) {
    stats[i].delta = 0;
    stats[i].support = 0;
    if (i == NJ->root || i < NJ->nSeq) {
      stats[i].age = LargeAge;
      stats[i].subtreeAge = LargeAge;
    } else {
      stats[i].age = 0;
      stats[i].subtreeAge = 0;
    }
  }
  return(stats);
}

nni_stats_t *FreeNNIStats(nni_stats_t *stats, NJ_t *NJ) {
  return(myfree(stats, sizeof(nni_stats_t)*NJ->maxnode));
}

int FindSPRSteps(/*IN/OUT*/NJ_t *NJ, 
		 int nodeMove,	 /* the node to move multiple times */
		 int nodeAround, /* sibling or parent of node to NNI to start the chain */
		 /*IN/OUT*/profile_t **upProfiles,
		 /*OUT*/spr_step_t *steps,
		 int maxSteps,
		 bool bFirstAC) {
  int iStep;
  for (iStep = 0; iStep < maxSteps; iStep++) {
    if (NJ->child[nodeAround].nChild != 2)
      break;			/* no further to go */

    /* Consider the NNIs around nodeAround */
    profile_t *profiles[4];
    int nodeABCD[4];
    SetupABCD(NJ, nodeAround, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);
    double criteria[3];
    (void) ChooseNNI(profiles, NJ->distance_matrix, NJ->nPos, NJ->nConstraints,
		     /*OUT*/criteria);

    /* Do & save the swap */
    spr_step_t *step = &steps[iStep];
    if (iStep == 0 ? bFirstAC : criteria[ACvsBD] < criteria[ADvsBC]) {
      /* swap B & C to put AC together */
      step->deltaLength = criteria[ACvsBD] - criteria[ABvsCD];
      step->nodes[0] = nodeABCD[1];
      step->nodes[1] = nodeABCD[2];
    } else {
      /* swap AC to put AD together */
      step->deltaLength = criteria[ADvsBC] - criteria[ABvsCD];
      step->nodes[0] = nodeABCD[0];
      step->nodes[1] = nodeABCD[2];
    }

    if (verbose>3) {
      fprintf(stderr, "SPR chain step %d for %d around %d swap %d %d deltaLen %.5f\n",
	      iStep+1, nodeAround, nodeMove, step->nodes[0], step->nodes[1], step->deltaLength);
      if (verbose>4)
	PrintNJInternal(stderr, NJ, /*useLen*/false);
    }
    ReplaceChild(/*IN/OUT*/NJ, nodeAround, step->nodes[0], step->nodes[1]);
    ReplaceChild(/*IN/OUT*/NJ, NJ->parent[nodeAround], step->nodes[1], step->nodes[0]);
    UpdateForNNI(/*IN/OUT*/NJ, nodeAround, /*IN/OUT*/upProfiles, /*useML*/false);

    /* set the new nodeAround -- either parent(nodeMove) or sibling(nodeMove) --
       so that it different from current nodeAround
     */
    int newAround[2] = { NJ->parent[nodeMove], Sibling(NJ, nodeMove) };
    if (NJ->parent[nodeMove] == NJ->root)
      RootSiblings(NJ, nodeMove, /*OUT*/newAround);
    assert(newAround[0] == nodeAround || newAround[1] == nodeAround);
    assert(newAround[0] != newAround[1]);
    nodeAround = newAround[newAround[0] == nodeAround ? 1 : 0];
  }
  return(iStep);
}

void UnwindSPRStep(/*IN/OUT*/NJ_t *NJ,
		   /*IN*/spr_step_t *step,
		   /*IN/OUT*/profile_t **upProfiles) {
  int parents[2];
  int i;
  for (i = 0; i < 2; i++) {
    assert(step->nodes[i] >= 0 && step->nodes[i] < NJ->maxnodes);
    parents[i] = NJ->parent[step->nodes[i]];
    assert(parents[i] >= 0);
  }
  assert(parents[0] != parents[1]);
  ReplaceChild(/*IN/OUT*/NJ, parents[0], step->nodes[0], step->nodes[1]);
  ReplaceChild(/*IN/OUT*/NJ, parents[1], step->nodes[1], step->nodes[0]);
  int iYounger = 0;
  if (NJ->parent[parents[0]] == parents[1]) {
    iYounger = 0;
  } else {
    assert(NJ->parent[parents[1]] == parents[0]);
    iYounger = 1;
  }
  UpdateForNNI(/*IN/OUT*/NJ, parents[iYounger], /*IN/OUT*/upProfiles, /*useML*/false);
}

/* Update the profile of node and its ancestor, and delete nearby out-profiles */
void UpdateForNNI(/*IN/OUT*/NJ_t *NJ, int node, /*IN/OUT*/profile_t **upProfiles,
		  bool useML) {
  int i;
  if (slow) {
    /* exhaustive update */
    for (i = 0; i < NJ->maxnodes; i++)
      DeleteUpProfile(upProfiles, NJ, i);

    /* update profiles back to root */
    int ancestor;
    for (ancestor = node; ancestor >= 0; ancestor = NJ->parent[ancestor])
      RecomputeProfile(/*IN/OUT*/NJ, upProfiles, ancestor, useML);

    /* remove any up-profiles made while doing that*/
    for (i = 0; i < NJ->maxnodes; i++)
      DeleteUpProfile(upProfiles, NJ, i);
  } else {
    /* if fast, only update around self
       note that upProfile(parent) is still OK after an NNI, but
       up-profiles of uncles may not be
    */
    DeleteUpProfile(upProfiles, NJ, node);
    for (i = 0; i < NJ->child[node].nChild; i++)
      DeleteUpProfile(upProfiles, NJ, NJ->child[node].child[i]);
    assert(node != NJ->root);
    int parent = NJ->parent[node];
    int neighbors[2] = { parent, Sibling(NJ, node) };
    if (parent == NJ->root)
      RootSiblings(NJ, node, /*OUT*/neighbors);
    DeleteUpProfile(upProfiles, NJ, neighbors[0]);
    DeleteUpProfile(upProfiles, NJ, neighbors[1]);
    int uncle = Sibling(NJ, parent);
    if (uncle >= 0)
      DeleteUpProfile(upProfiles, NJ, uncle);
    RecomputeProfile(/*IN/OUT*/NJ, upProfiles, node, useML);
    RecomputeProfile(/*IN/OUT*/NJ, upProfiles, parent, useML);
  }
}

void SPR(/*IN/OUT*/NJ_t *NJ, int maxSPRLength, int iRound, int nRounds) {
  /* Given a non-root node N with children A,B, sibling C, and uncle D,
     we can try to move A by doing three types of moves (4 choices):
     "down" -- swap A with a child of B (if B is not a leaf) [2 choices]
     "over" -- swap B with C
     "up" -- swap A with D
     We follow down moves with down moves, over moves with down moves, and
     up moves with either up or over moves. (Other choices are just backing
     up and hence useless.)

     As with NNIs, we keep track of up-profiles as we go. However, some of the regular
     profiles may also become "stale" so it is a bit trickier.

     We store the traversal before we do SPRs to avoid any possible infinite loop
  */
  double last_tot_len = 0.0;
  if (NJ->nSeq <= 3 || maxSPRLength < 1)
    return;
  if (slow)
    last_tot_len = TreeLength(NJ, /*recomputeLengths*/true);
  int *nodeList = mymalloc(sizeof(int) * NJ->maxnodes);
  int nodeListLen = 0;
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    nodeList[nodeListLen++] = node;
  }
  assert(nodeListLen == NJ->maxnode);
  traversal = FreeTraversal(traversal,NJ);

  profile_t **upProfiles = UpProfiles(NJ);
  spr_step_t *steps = mymalloc(sizeof(spr_step_t) * maxSPRLength); /* current chain of SPRs */

  int i;
  for (i = 0; i < nodeListLen; i++) {
    node = nodeList[i];
    if ((i % 100) == 0)
      ProgressReport("SPR round %3d of %3d, %d of %d nodes",
		     iRound+1, nRounds, i+1, nodeListLen);
    if (node == NJ->root)
      continue; /* nothing to do for root */
    /* The nodes to NNI around */
    int nodeAround[2] = { NJ->parent[node], Sibling(NJ, node) };
    if (NJ->parent[node] == NJ->root) {
      /* NNI around both siblings instead */
      RootSiblings(NJ, node, /*OUT*/nodeAround);
    }
    bool bChanged = false;
    int iAround;
    for (iAround = 0; iAround < 2 && bChanged == false; iAround++) {
      int ACFirst;
      for (ACFirst = 0; ACFirst < 2 && bChanged == false; ACFirst++) {
	if(verbose > 3)
	  PrintNJInternal(stderr, NJ, /*useLen*/false);
	int chainLength = FindSPRSteps(/*IN/OUT*/NJ, node, nodeAround[iAround],
				       upProfiles, /*OUT*/steps, maxSPRLength, (bool)ACFirst);
	double dMinDelta = 0.0;
	int iCBest = -1;
	double dTotDelta = 0.0;
	int iC;
	for (iC = 0; iC < chainLength; iC++) {
	  dTotDelta += steps[iC].deltaLength;
	  if (dTotDelta < dMinDelta) {
	    dMinDelta = dTotDelta;
	    iCBest = iC;
	  }
	}
      
	if (verbose>3) {
	  fprintf(stderr, "SPR %s %d around %d chainLength %d of %d deltaLength %.5f swaps:",
		  iCBest >= 0 ? "move" : "abandoned",
		  node,nodeAround[iAround],iCBest+1,chainLength,dMinDelta);
	  for (iC = 0; iC < chainLength; iC++)
	    fprintf(stderr, " (%d,%d)%.4f", steps[iC].nodes[0], steps[iC].nodes[1], steps[iC].deltaLength);
	  fprintf(stderr,"\n");
	}
	for (iC = chainLength - 1; iC > iCBest; iC--)
	  UnwindSPRStep(/*IN/OUT*/NJ, /*IN*/&steps[iC], /*IN/OUT*/upProfiles);
	if(verbose > 3)
	  PrintNJInternal(stderr, NJ, /*useLen*/false);
	while (slow && iCBest >= 0) {
	  double expected_tot_len = last_tot_len + dMinDelta;
	  double new_tot_len = TreeLength(NJ, /*recompute*/true);
	  if (verbose > 2)
	    fprintf(stderr, "Total branch-length is now %.4f was %.4f expected %.4f\n",
		    new_tot_len, last_tot_len, expected_tot_len);
	  if (new_tot_len < last_tot_len) {
	    last_tot_len = new_tot_len;
	    break;		/* no rewinding necessary */
	  }
	  if (verbose > 2)
	    fprintf(stderr, "Rewinding SPR to %d\n",iCBest);
	  UnwindSPRStep(/*IN/OUT*/NJ, /*IN*/&steps[iCBest], /*IN/OUT*/upProfiles);
	  dMinDelta -= steps[iCBest].deltaLength;
	  iCBest--;
	}
	if (iCBest >= 0)
	  bChanged = true;
      }	/* loop over which step to take at 1st NNI */
    } /* loop over which node to pivot around */

    if (bChanged) {
      nSPR++;		/* the SPR move is OK */
      /* make sure all the profiles are OK */
      int j;
      for (j = 0; j < NJ->maxnodes; j++)
	DeleteUpProfile(upProfiles, NJ, j);
      int ancestor;
      for (ancestor = NJ->parent[node]; ancestor >= 0; ancestor = NJ->parent[ancestor])
	RecomputeProfile(/*IN/OUT*/NJ, upProfiles, ancestor, /*useML*/false);
    }
  } /* end loop over subtrees to prune & regraft */
  steps = myfree(steps, sizeof(spr_step_t) * maxSPRLength);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  nodeList = myfree(nodeList, sizeof(int) * NJ->maxnodes);
}

void RecomputeProfile(/*IN/OUT*/NJ_t *NJ, /*IN/OUT*/profile_t **upProfiles, int node,
		      bool useML) {
  if (node < NJ->nSeq || node == NJ->root)
    return;			/* no profile to compute */
  assert(NJ->child[node].nChild==2);

  profile_t *profiles[4];
  double weight = 0.5;
  if (useML || !bionj) {
    profiles[0] = NJ->profiles[NJ->child[node].child[0]];
    profiles[1] = NJ->profiles[NJ->child[node].child[1]];
  } else {
    int nodeABCD[4];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, useML);
    weight = QuartetWeight(profiles, NJ->distance_matrix, NJ->nPos);
  }
  if (verbose>3) {
    if (useML) {
      fprintf(stderr, "Recompute %d from %d %d lengths %.4f %.4f\n",
	      node,
	      NJ->child[node].child[0],
	      NJ->child[node].child[1],
	      NJ->branchlength[NJ->child[node].child[0]],
	      NJ->branchlength[NJ->child[node].child[1]]);
    } else {
      fprintf(stderr, "Recompute %d from %d %d weight %.3f\n",
	      node, NJ->child[node].child[0], NJ->child[node].child[1], weight);
    }
  }
  NJ->profiles[node] = FreeProfile(NJ->profiles[node], NJ->nPos, NJ->nConstraints);
  if (useML) {
    NJ->profiles[node] = PosteriorProfile(profiles[0], profiles[1],
					  NJ->branchlength[NJ->child[node].child[0]],
					  NJ->branchlength[NJ->child[node].child[1]],
					  NJ->transmat, &NJ->rates, NJ->nPos, NJ->nConstraints);
  } else {
    NJ->profiles[node] = AverageProfile(profiles[0], profiles[1],
					NJ->nPos, NJ->nConstraints,
					NJ->distance_matrix, weight);
  }
}

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

/* Resets the children entry of parent and also the parent entry of newchild */
void ReplaceChild(/*IN/OUT*/NJ_t *NJ, int parent, int oldchild, int newchild) {
  NJ->parent[newchild] = parent;

  int iChild;
  for (iChild = 0; iChild < NJ->child[parent].nChild; iChild++) {
    if (NJ->child[parent].child[iChild] == oldchild) {
      NJ->child[parent].child[iChild] = newchild;
      return;
    }
  }
  assert(0);
}

/* Recomputes all branch lengths

   For internal branches such as (A,B) vs. (C,D), uses the formula 

   length(AB|CD) = (d(A,C)+d(A,D)+d(B,C)+d(B,D))/4 - d(A,B)/2 - d(C,D)/2

   (where all distances are profile distances - diameters).

   For external branches (e.g. to leaves) A vs. (B,C), use the formula

   length(A|BC) = (d(A,B)+d(A,C)-d(B,C))/2
*/
void UpdateBranchLengths(/*IN/OUT*/NJ_t *NJ) {
  if (NJ->nSeq < 2)
    return;
  else if (NJ->nSeq == 2) {
    int root = NJ->root;
    int nodeA = NJ->child[root].child[0];
    int nodeB = NJ->child[root].child[1];
    besthit_t h;
    ProfileDist(NJ->profiles[nodeA],NJ->profiles[nodeB],
		NJ->nPos, NJ->distance_matrix, /*OUT*/&h);
    if (logdist)
      h.dist = LogCorrect(h.dist);
    NJ->branchlength[nodeA] = h.dist/2.0;
    NJ->branchlength[nodeB] = h.dist/2.0;
    return;
  }

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    /* reset branch length of node (distance to its parent) */
    if (node == NJ->root)
      continue; /* no branch length to set */
    if (node < NJ->nSeq) { /* a leaf */
      profile_t *profileA = NJ->profiles[node];
      profile_t *profileB = NULL;
      profile_t *profileC = NULL;

      int sib = Sibling(NJ,node);
      if (sib == -1) { /* at root, have 2 siblings */
	int sibs[2];
	RootSiblings(NJ, node, /*OUT*/sibs);
	profileB = NJ->profiles[sibs[0]];
	profileC = NJ->profiles[sibs[1]];
      } else {
	profileB = NJ->profiles[sib];
	profileC = GetUpProfile(/*IN/OUT*/upProfiles, NJ, NJ->parent[node], /*useML*/false);
      }
      profile_t *profiles[3] = {profileA,profileB,profileC};
      double d[3]; /*AB,AC,BC*/
      CorrectedPairDistances(profiles, 3, NJ->distance_matrix, NJ->nPos, /*OUT*/d);
      /* d(A,BC) = (dAB+dAC-dBC)/2 */
      NJ->branchlength[node] = (d[0]+d[1]-d[2])/2.0;
    } else {
      profile_t *profiles[4];
      int nodeABCD[4];
      SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);
      double d[6];
      CorrectedPairDistances(profiles, 4, NJ->distance_matrix, NJ->nPos, /*OUT*/d);
      NJ->branchlength[node] = (d[qAC]+d[qAD]+d[qBC]+d[qBD])/4.0 - (d[qAB]+d[qCD])/2.0;
      
      /* no longer needed */
      DeleteUpProfile(upProfiles, NJ, nodeABCD[0]);
      DeleteUpProfile(upProfiles, NJ, nodeABCD[1]);
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}

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

void ReliabilityNJ(/*IN/OUT*/NJ_t *NJ, int nBootstrap) {
  /* For each non-root node N, with children A,B, parent P, sibling C, and grandparent G,
     we test the reliability of the split (A,B) versus rest by comparing the profiles
     of A, B, C, and the "up-profile" of P.

     Each node's upProfile is the average of its sibling's (down)-profile + its parent's up-profile
     (If node's parent is the root, then there are two siblings and we don't need an up-profile)

     To save memory, we do depth-first-search down from the root, and we only keep
     up-profiles for nodes in the active path.
  */
  if (NJ->nSeq <= 3 || nBootstrap <= 0)
    return;			/* nothing to do */
  int *col = ResampleColumns(NJ->nPos, nBootstrap);

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  int iNodesDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */

    if(iNodesDone > 0 && (iNodesDone % 100) == 0)
      ProgressReport("Local bootstrap for %6d of %6d internal splits", iNodesDone, NJ->nSeq-3, 0, 0);
    iNodesDone++;

    profile_t *profiles[4];
    int nodeABCD[4];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);

    NJ->support[node] = SplitSupport(profiles[0], profiles[1], profiles[2], profiles[3],
				     NJ->distance_matrix,
				     NJ->nPos,
				     nBootstrap,
				     col);

    /* no longer needed */
    DeleteUpProfile(upProfiles, NJ, nodeABCD[0]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[1]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[2]);
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  col = myfree(col, sizeof(int)*((size_t)NJ->nPos)*nBootstrap);
}

profile_t *NewProfile(int nPos, int nConstraints) {
  profile_t *profile = (profile_t *)mymalloc(sizeof(profile_t));
  profile->weights = mymalloc(sizeof(numeric_t)*nPos);
  profile->codes = mymalloc(sizeof(unsigned char)*nPos);
  profile->vectors = NULL;
  profile->nVectors = 0;
  profile->codeDist = NULL;
  if (nConstraints == 0) {
    profile->nOn = NULL;
    profile->nOff = NULL;
  } else {
    profile->nOn = mymalloc(sizeof(int)*nConstraints);
    profile->nOff = mymalloc(sizeof(int)*nConstraints);
  }
  return(profile);
}

profile_t *FreeProfile(profile_t *profile, int nPos, int nConstraints) {
    if(profile==NULL) return(NULL);
    myfree(profile->codes, nPos);
    myfree(profile->weights, nPos);
    myfree(profile->vectors, sizeof(numeric_t)*nCodes*profile->nVectors);
    myfree(profile->codeDist, sizeof(numeric_t)*nCodes*nPos);
    if (nConstraints > 0) {
      myfree(profile->nOn, sizeof(int)*nConstraints);
      myfree(profile->nOff,  sizeof(int)*nConstraints);
    }
    return(myfree(profile, sizeof(profile_t)));
}

void SetupABCD(NJ_t *NJ, int node,
	       /* the 4 profiles; the last one is an outprofile */
	       /*OPTIONAL OUT*/profile_t *profiles[4], 
	       /*OPTIONAL IN/OUT*/profile_t **upProfiles,
	       /*OUT*/int nodeABCD[4],
	       bool useML) {
  int parent = NJ->parent[node];
  assert(parent >= 0);
  assert(NJ->child[node].nChild == 2);
  nodeABCD[0] = NJ->child[node].child[0]; /*A*/
  nodeABCD[1] = NJ->child[node].child[1]; /*B*/

  profile_t *profile4 = NULL;
  if (parent == NJ->root) {
    int sibs[2];
    RootSiblings(NJ, node, /*OUT*/sibs);
    nodeABCD[2] = sibs[0];
    nodeABCD[3] = sibs[1];
    if (profiles == NULL)
      return;
    profile4 = NJ->profiles[sibs[1]];
  } else {
    nodeABCD[2] = Sibling(NJ,node);
    assert(nodeABCD[2] >= 0);
    nodeABCD[3] = parent;
    if (profiles == NULL)
      return;
    profile4 = GetUpProfile(upProfiles,NJ,parent,useML);
  }
  assert(upProfiles != NULL);
  int i;
  for (i = 0; i < 3; i++)
    profiles[i] = NJ->profiles[nodeABCD[i]];
  profiles[3] = profile4;
}


int Sibling(NJ_t *NJ, int node) {
  int parent = NJ->parent[node];
  if (parent < 0 || parent == NJ->root)
    return(-1);
  int iChild;
  for(iChild=0;iChild<NJ->child[parent].nChild;iChild++) {
    if(NJ->child[parent].child[iChild] != node)
      return (NJ->child[parent].child[iChild]);
  }
  assert(0);
  return(-1);
}

void RootSiblings(NJ_t *NJ, int node, /*OUT*/int sibs[2]) {
  assert(NJ->parent[node] == NJ->root);
  assert(NJ->child[NJ->root].nChild == 3);

  int nSibs = 0;
  int iChild;
  for(iChild=0; iChild < NJ->child[NJ->root].nChild; iChild++) {
    int child = NJ->child[NJ->root].child[iChild];
    if (child != node) sibs[nSibs++] = child;
  }
  assert(nSibs==2);
}

void TestSplitsML(/*IN/OUT*/NJ_t *NJ, /*OUT*/SplitCount_t *splitcount, int nBootstrap) {
  const double tolerance = 1e-6;
  splitcount->nBadSplits = 0;
  splitcount->nConstraintViolations = 0;
  splitcount->nBadBoth = 0;
  splitcount->nSplits = 0;
  splitcount->dWorstDeltaUnconstrained = 0;
  splitcount->dWorstDeltaConstrained = 0;

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  int *col = nBootstrap > 0 ? ResampleColumns(NJ->nPos, nBootstrap) : NULL;
  double *site_likelihoods[3];
  int choice;
  for (choice = 0; choice < 3; choice++)
    site_likelihoods[choice] = mymalloc(sizeof(double)*NJ->nPos);

  int iNodesDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */
    
    if(iNodesDone > 0 && (iNodesDone % 100) == 0)
      ProgressReport("ML split tests for %6d of %6d internal splits", iNodesDone, NJ->nSeq-3, 0, 0);
    iNodesDone++;

    profile_t *profiles[4];
    int nodeABCD[4];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/true);
    double loglk[3];
    double len[5];
    int i;
    for (i = 0; i < 4; i++)
      len[i] = NJ->branchlength[nodeABCD[i]];
    len[4] = NJ->branchlength[node];
    double lenABvsCD[5] = {len[LEN_A], len[LEN_B], len[LEN_C], len[LEN_D], len[LEN_I]};
    double lenACvsBD[5] = {len[LEN_A], len[LEN_C], len[LEN_B], len[LEN_D], len[LEN_I]};   /* Swap B & C */
    double lenADvsBC[5] = {len[LEN_A], len[LEN_D], len[LEN_C], len[LEN_B], len[LEN_I]};   /* Swap B & D */

    {
#ifdef OPENMP
      #pragma omp parallel
      #pragma omp sections
#endif
      {
#ifdef OPENMP
      #pragma omp section
#endif
	{
	  /* Lengths are already optimized for ABvsCD */
	  loglk[ABvsCD] = MLQuartetLogLk(profiles[0], profiles[1], profiles[2], profiles[3],
					 NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenABvsCD,
					 /*OUT*/site_likelihoods[ABvsCD]);
	}

#ifdef OPENMP
      #pragma omp section
#endif
	{
	  loglk[ACvsBD] = MLQuartetOptimize(profiles[0], profiles[2], profiles[1], profiles[3],
					    NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenACvsBD, /*pStarTest*/NULL,
					    /*OUT*/site_likelihoods[ACvsBD]);
	}

#ifdef OPENMP
      #pragma omp section
#endif
	{
	  loglk[ADvsBC] = MLQuartetOptimize(profiles[0], profiles[3], profiles[2], profiles[1],
					    NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenADvsBC, /*pStarTest*/NULL,
					    /*OUT*/site_likelihoods[ADvsBC]);
	}
      }
    }

    /* do a second pass on the better alternative if it is close */
    if (loglk[ACvsBD] > loglk[ADvsBC]) {
      if (mlAccuracy > 1 || loglk[ACvsBD] > loglk[ABvsCD] - closeLogLkLimit) {
	loglk[ACvsBD] = MLQuartetOptimize(profiles[0], profiles[2], profiles[1], profiles[3],
					  NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenACvsBD, /*pStarTest*/NULL,
					  /*OUT*/site_likelihoods[ACvsBD]);
      }
    } else {
      if (mlAccuracy > 1 || loglk[ADvsBC] > loglk[ABvsCD] - closeLogLkLimit) {
	loglk[ADvsBC] = MLQuartetOptimize(profiles[0], profiles[3], profiles[2], profiles[1],
					  NJ->nPos, NJ->transmat, &NJ->rates, /*IN/OUT*/lenADvsBC, /*pStarTest*/NULL,
					  /*OUT*/site_likelihoods[ADvsBC]);
      }
    }

    if (loglk[ABvsCD] >= loglk[ACvsBD] && loglk[ABvsCD] >= loglk[ADvsBC])
      choice = ABvsCD;
    else if (loglk[ACvsBD] >= loglk[ABvsCD] && loglk[ACvsBD] >= loglk[ADvsBC])
      choice = ACvsBD;
    else
      choice = ADvsBC;
    bool badSplit = loglk[choice] > loglk[ABvsCD] + treeLogLkDelta; /* ignore small changes in likelihood */

    /* constraint penalties, indexed by nni_t (lower is better) */
    double p[3];
    QuartetConstraintPenalties(profiles, NJ->nConstraints, /*OUT*/p);
    bool bBadConstr = p[ABvsCD] > p[ACvsBD] + tolerance || p[ABvsCD] > p[ADvsBC] + tolerance;
    bool violateConstraint = false;
    int iC;
    for (iC=0; iC < NJ->nConstraints; iC++) {
      if (SplitViolatesConstraint(profiles, iC)) {
	violateConstraint = true;
	break;
      }
    }
    splitcount->nSplits++;
    if (violateConstraint)
      splitcount->nConstraintViolations++;
    if (badSplit)
      splitcount->nBadSplits++;
    if (badSplit && bBadConstr)
      splitcount->nBadBoth++;
    if (badSplit) {
      double delta = loglk[choice] - loglk[ABvsCD];
      /* If ABvsCD is favored over the more likely NNI by constraints,
	 then this is probably a bad split because of the constraint */
      if (p[choice] > p[ABvsCD] + tolerance)
	splitcount->dWorstDeltaConstrained = MAX(delta, splitcount->dWorstDeltaConstrained);
      else
	splitcount->dWorstDeltaUnconstrained = MAX(delta, splitcount->dWorstDeltaUnconstrained);
    }
    if (nBootstrap>0)
      NJ->support[node] = badSplit ? 0.0 : SHSupport(NJ->nPos, nBootstrap, col, loglk, site_likelihoods);

    /* No longer needed */
    DeleteUpProfile(upProfiles, NJ, nodeABCD[0]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[1]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[2]);
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  if (nBootstrap>0)
    col = myfree(col, sizeof(int)*((size_t)NJ->nPos)*nBootstrap);
  for (choice = 0; choice < 3; choice++)
    site_likelihoods[choice] = myfree(site_likelihoods[choice], sizeof(double)*NJ->nPos);
}
    

void TestSplitsMinEvo(NJ_t *NJ, /*OUT*/SplitCount_t *splitcount) {
  const double tolerance = 1e-6;
  splitcount->nBadSplits = 0;
  splitcount->nConstraintViolations = 0;
  splitcount->nBadBoth = 0;
  splitcount->nSplits = 0;
  splitcount->dWorstDeltaUnconstrained = 0.0;
  splitcount->dWorstDeltaConstrained = 0.0;

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal, /*pUp*/NULL)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */

    profile_t *profiles[4];
    int nodeABCD[4];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);

    if (verbose>2)
      fprintf(stderr,"Testing Split around %d: A=%d B=%d C=%d D=up(%d) or node parent %d\n",
	      node, nodeABCD[0], nodeABCD[1], nodeABCD[2], nodeABCD[3], NJ->parent[node]);

    double d[6];		/* distances, perhaps log-corrected distances, no constraint penalties */
    CorrectedPairDistances(profiles, 4, NJ->distance_matrix, NJ->nPos, /*OUT*/d);

    /* alignment-based scores for each split (lower is better) */
    double sABvsCD = d[qAB] + d[qCD];
    double sACvsBD = d[qAC] + d[qBD];
    double sADvsBC = d[qAD] + d[qBC];

    /* constraint penalties, indexed by nni_t (lower is better) */
    double p[3];
    QuartetConstraintPenalties(profiles, NJ->nConstraints, /*OUT*/p);

    int nConstraintsViolated = 0;
    int iC;
    for (iC=0; iC < NJ->nConstraints; iC++) {
      if (SplitViolatesConstraint(profiles, iC)) {
	nConstraintsViolated++;
	if (verbose > 2) {
	  double penalty[3] = {0.0,0.0,0.0};
	  (void)QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/penalty);
	  fprintf(stderr, "Violate constraint %d at %d (children %d %d) penalties %.3f %.3f %.3f %d/%d %d/%d %d/%d %d/%d\n",
		  iC, node, NJ->child[node].child[0], NJ->child[node].child[1],
		  penalty[ABvsCD], penalty[ACvsBD], penalty[ADvsBC],
		  profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		  profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		  profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		  profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
	}
      }
    }

    double delta = sABvsCD - MIN(sACvsBD,sADvsBC);
    bool bBadDist = delta > tolerance;
    bool bBadConstr = p[ABvsCD] > p[ACvsBD] + tolerance || p[ABvsCD] > p[ADvsBC] + tolerance;

    splitcount->nSplits++;
    if (bBadDist) {
      nni_t choice = sACvsBD < sADvsBC ? ACvsBD : ADvsBC;
      /* If ABvsCD is favored over the shorter NNI by constraints,
	 then this is probably a bad split because of the constraint */
      if (p[choice] > p[ABvsCD] + tolerance)
	splitcount->dWorstDeltaConstrained = MAX(delta, splitcount->dWorstDeltaConstrained);
      else
	splitcount->dWorstDeltaUnconstrained = MAX(delta, splitcount->dWorstDeltaUnconstrained);
    }
	    
    if (nConstraintsViolated > 0)
      splitcount->nConstraintViolations++; /* count splits with any violations, not #constraints in a splits */
    if (bBadDist)
      splitcount->nBadSplits++;
    if (bBadDist && bBadConstr)
      splitcount->nBadBoth++;
    if (bBadConstr && verbose > 2) {
      /* Which NNI would be better */
      double dist_advantage = 0;
      double constraint_penalty = 0;
      if (p[ACvsBD] < p[ADvsBC]) {
	dist_advantage = sACvsBD - sABvsCD;
	constraint_penalty = p[ABvsCD] - p[ACvsBD];
      } else {
	dist_advantage = sADvsBC - sABvsCD;
	constraint_penalty = p[ABvsCD] - p[ADvsBC];
      }
      fprintf(stderr, "Violate constraints %d distance_advantage %.3f constraint_penalty %.3f (children %d %d):",
	      node, dist_advantage, constraint_penalty,
	      NJ->child[node].child[0], NJ->child[node].child[1]);
      /* list the constraints with a penalty, meaning that ABCD all have non-zero
         values and that AB|CD worse than others */
      for (iC = 0; iC < NJ->nConstraints; iC++) {
	double ppart[6];
	if (QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/ppart)) {
	  if (ppart[qAB] + ppart[qCD] > ppart[qAD] + ppart[qBC] + tolerance
	      || ppart[qAB] + ppart[qCD] > ppart[qAC] + ppart[qBD] + tolerance)
	    fprintf(stderr, " %d (%d/%d %d/%d %d/%d %d/%d)", iC,
		    profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		    profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		    profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		    profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
	}
      }
      fprintf(stderr, "\n");
    }
    
    /* no longer needed */
    DeleteUpProfile(upProfiles, NJ, nodeABCD[0]);
    DeleteUpProfile(upProfiles, NJ, nodeABCD[1]);
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}

/* Computes support for (A,B),(C,D) compared to that for (A,C),(B,D) and (A,D),(B,C) */
double SplitSupport(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
		    /*OPTIONAL*/distance_matrix_t *dmat,
		    int nPos,
		    int nBootstrap,
		    int *col) {
  int i,j;
  long lPos = nPos; 		/* to avoid overflow when multiplying */

  /* Note distpieces are weighted */
  double *distpieces[6];
  double *weights[6];
  for (j = 0; j < 6; j++) {
    distpieces[j] = (double*)mymalloc(sizeof(double)*nPos);
    weights[j] = (double*)mymalloc(sizeof(double)*nPos);
  }

  int iFreqA = 0;
  int iFreqB = 0;
  int iFreqC = 0;
  int iFreqD = 0;
  for (i = 0; i < nPos; i++) {
    numeric_t *fA = GET_FREQ(pA, i, /*IN/OUT*/iFreqA);
    numeric_t *fB = GET_FREQ(pB, i, /*IN/OUT*/iFreqB);
    numeric_t *fC = GET_FREQ(pC, i, /*IN/OUT*/iFreqC);
    numeric_t *fD = GET_FREQ(pD, i, /*IN/OUT*/iFreqD);

    weights[qAB][i] = pA->weights[i] * pB->weights[i];
    weights[qAC][i] = pA->weights[i] * pC->weights[i];
    weights[qAD][i] = pA->weights[i] * pD->weights[i];
    weights[qBC][i] = pB->weights[i] * pC->weights[i];
    weights[qBD][i] = pB->weights[i] * pD->weights[i];
    weights[qCD][i] = pC->weights[i] * pD->weights[i];

    distpieces[qAB][i] = weights[qAB][i] * ProfileDistPiece(pA->codes[i], pB->codes[i], fA, fB, dmat, NULL);
    distpieces[qAC][i] = weights[qAC][i] * ProfileDistPiece(pA->codes[i], pC->codes[i], fA, fC, dmat, NULL);
    distpieces[qAD][i] = weights[qAD][i] * ProfileDistPiece(pA->codes[i], pD->codes[i], fA, fD, dmat, NULL);
    distpieces[qBC][i] = weights[qBC][i] * ProfileDistPiece(pB->codes[i], pC->codes[i], fB, fC, dmat, NULL);
    distpieces[qBD][i] = weights[qBD][i] * ProfileDistPiece(pB->codes[i], pD->codes[i], fB, fD, dmat, NULL);
    distpieces[qCD][i] = weights[qCD][i] * ProfileDistPiece(pC->codes[i], pD->codes[i], fC, fD, dmat, NULL);
  }
  assert(iFreqA == pA->nVectors);
  assert(iFreqB == pB->nVectors);
  assert(iFreqC == pC->nVectors);
  assert(iFreqD == pD->nVectors);

  double totpieces[6];
  double totweights[6];
  double dists[6];
  for (j = 0; j < 6; j++) {
    totpieces[j] = 0.0;
    totweights[j] = 0.0;
    for (i = 0; i < nPos; i++) {
      totpieces[j] += distpieces[j][i];
      totweights[j] += weights[j][i];
    }
    dists[j] = totweights[j] > 0.01 ? totpieces[j]/totweights[j] : 3.0;
    if (logdist)
      dists[j] = LogCorrect(dists[j]);
  }

  /* Support1 = Support(AB|CD over AC|BD) = d(A,C)+d(B,D)-d(A,B)-d(C,D)
     Support2 = Support(AB|CD over AD|BC) = d(A,D)+d(B,C)-d(A,B)-d(C,D)
  */
  double support1 = dists[qAC] + dists[qBD] - dists[qAB] - dists[qCD];
  double support2 = dists[qAD] + dists[qBC] - dists[qAB] - dists[qCD];

  if (support1 < 0 || support2 < 0) {
    nSuboptimalSplits++;	/* Another split seems superior */
  }

  assert(nBootstrap > 0);
  int nSupport = 0;

  int iBoot;
  for (iBoot=0;iBoot<nBootstrap;iBoot++) {
    int *colw = &col[lPos*iBoot];

    for (j = 0; j < 6; j++) {
      double totp = 0;
      double totw = 0;
      double *d = distpieces[j];
      double *w = weights[j];
      for (i=0; i<nPos; i++) {
	int c = colw[i];
	totp += d[c];
	totw += w[c];
      }
      dists[j] = totw > 0.01 ? totp/totw : 3.0;
      if (logdist)
	dists[j] = LogCorrect(dists[j]);
    }
    support1 = dists[qAC] + dists[qBD] - dists[qAB] - dists[qCD];
    support2 = dists[qAD] + dists[qBC] - dists[qAB] - dists[qCD];
    if (support1 > 0 && support2 > 0)
      nSupport++;
  } /* end loop over bootstrap replicates */

  for (j = 0; j < 6; j++) {
    distpieces[j] = myfree(distpieces[j], sizeof(double)*nPos);
    weights[j] = myfree(weights[j], sizeof(double)*nPos);
  }
  return( nSupport/(double)nBootstrap );
}

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


void SetDistCriterion(/*IN/OUT*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *hit) {
  if (hit->i < NJ->nSeq && hit->j < NJ->nSeq) {
    SeqDist(NJ->profiles[hit->i]->codes,
	    NJ->profiles[hit->j]->codes,
	    NJ->nPos, NJ->distance_matrix, /*OUT*/hit);
  } else {
    ProfileDist(NJ->profiles[hit->i],
		NJ->profiles[hit->j],
		NJ->nPos, NJ->distance_matrix, /*OUT*/hit);
    hit->dist -= (NJ->diameter[hit->i] + NJ->diameter[hit->j]);
  }
  hit->dist += constraintWeight
    * (double)JoinConstraintPenalty(NJ, hit->i, hit->j);
  SetCriterion(NJ,nActive,/*IN/OUT*/hit);
}

void SetCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join) {
  if(join->i < 0
     || join->j < 0
     || NJ->parent[join->i] >= 0
     || NJ->parent[join->j] >= 0)
    return;
  assert(NJ->nOutDistActive[join->i] >= nActive);
  assert(NJ->nOutDistActive[join->j] >= nActive);

  int nDiffAllow = tophitsMult > 0 ? (int)(nActive*staleOutLimit) : 0;
  if (NJ->nOutDistActive[join->i] - nActive > nDiffAllow)
    SetOutDistance(NJ, join->i, nActive);
  if (NJ->nOutDistActive[join->j] - nActive > nDiffAllow)
    SetOutDistance(NJ, join->j, nActive);
  double outI = NJ->outDistances[join->i];
  if (NJ->nOutDistActive[join->i] != nActive)
    outI *= (nActive-1)/(double)(NJ->nOutDistActive[join->i]-1);
  double outJ = NJ->outDistances[join->j];
  if (NJ->nOutDistActive[join->j] != nActive)
    outJ *= (nActive-1)/(double)(NJ->nOutDistActive[join->j]-1);
  join->criterion = join->dist - (outI+outJ)/(double)(nActive-2);
  if (verbose > 2 && nActive <= 5) {
    fprintf(stderr, "Set Criterion to join %d %d with nActive=%d dist+penalty %.3f criterion %.3f\n",
	    join->i, join->j, nActive, join->dist, join->criterion);
  }
}

int ActiveAncestor(/*IN*/NJ_t *NJ, int iNode) {
  if (iNode < 0)
    return(iNode);
  while(NJ->parent[iNode] >= 0)
    iNode = NJ->parent[iNode];
  return(iNode);
}

bool GetVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		/*IN/OUT*/top_hits_t *tophits,
		int iNode, /*OUT*/besthit_t *visible) {
  if (iNode < 0 || NJ->parent[iNode] >= 0)
    return(false);
  hit_t *v = &tophits->visible[iNode];
  if (v->j < 0 || NJ->parent[v->j] >= 0)
    return(false);
  *visible = HitToBestHit(iNode, *v);
  SetCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/visible);  
  return(true);
}

void UpdateVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		   /*IN*/besthit_t *tophitsNode,
		   int nTopHits,
		  /*IN/OUT*/top_hits_t *tophits) {
  int iHit;

  for(iHit = 0; iHit < nTopHits; iHit++) {
    besthit_t *hit = &tophitsNode[iHit];
    if (hit->i < 0) continue;	/* possible empty entries */
    assert(NJ->parent[hit->i] < 0);
    assert(hit->j >= 0 && NJ->parent[hit->j] < 0);
    besthit_t visible;
    bool bSuccess = GetVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, hit->j, /*OUT*/&visible);
    if (!bSuccess || hit->criterion < visible.criterion) {
      if (bSuccess)
	nVisibleUpdate++;
      hit_t *v = &tophits->visible[hit->j];
      v->j = hit->i;
      v->dist = hit->dist;
      UpdateTopVisible(NJ, nActive, hit->j, v, /*IN/OUT*/tophits);
      if(verbose>5) fprintf(stderr,"NewVisible %d %d %f\n",
			    hit->j,v->j,v->dist);
    }
  } /* end loop over hits */
}

/* Update the top-visible list to perhaps include visible[iNode] */
void UpdateTopVisible(/*IN*/NJ_t * NJ, int nActive,
		      int iIn, /*IN*/hit_t *hit,
		      /*IN/OUT*/top_hits_t *tophits) {
  assert(tophits != NULL);
  bool bIn = false; 		/* placed in the list */
  int i;

  /* First, if the list is not full, put it in somewhere */
  for (i = 0; i < tophits->nTopVisible && !bIn; i++) {
    int iNode = tophits->topvisible[i];
    if (iNode == iIn) {
      /* this node is already in the top hit list */
      bIn = true;
    } else if (iNode < 0 || NJ->parent[iNode] >= 0) {
      /* found an empty spot */
      bIn = true;
      tophits->topvisible[i] = iIn;
    }
  }

  int iPosWorst = -1;
  double dCriterionWorst = -1e20;
  if (!bIn) {
    /* Search for the worst hit */
    for (i = 0; i < tophits->nTopVisible && !bIn; i++) {
      int iNode = tophits->topvisible[i];
      assert(iNode >= 0 && NJ->parent[iNode] < 0 && iNode != iIn);
      besthit_t visible;
      if (!GetVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, iNode, /*OUT*/&visible)) {
	/* found an empty spot */
	tophits->topvisible[i] = iIn;
	bIn = true;
      } else if (visible.i == hit->j && visible.j == iIn) {
	/* the reverse hit is already in the top hit list */
	bIn = true;
      } else if (visible.criterion >= dCriterionWorst) {
	iPosWorst = i;
	dCriterionWorst = visible.criterion;
      }
    }
  }

  if (!bIn && iPosWorst >= 0) {
    besthit_t visible = HitToBestHit(iIn, *hit);
    SetCriterion(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/&visible);
    if (visible.criterion < dCriterionWorst) {
      if (verbose > 2) {
	int iOld = tophits->topvisible[iPosWorst];
	fprintf(stderr, "TopVisible replace %d=>%d with %d=>%d\n",
		iOld, tophits->visible[iOld].j, visible.i, visible.j);
      }
      tophits->topvisible[iPosWorst] = iIn;
    }
  }

  if (verbose > 2) {
    fprintf(stderr, "Updated TopVisible: ");
    for (i = 0; i < tophits->nTopVisible; i++) {
      int iNode = tophits->topvisible[i];
      if (iNode >= 0 && NJ->parent[iNode] < 0) {
	besthit_t bh = HitToBestHit(iNode, tophits->visible[iNode]);
	SetDistCriterion(NJ, nActive, &bh);
	fprintf(stderr, " %d=>%d:%.4f", bh.i, bh.j, bh.criterion);
      }
    }
    fprintf(stderr,"\n");
  }
}

/* Recompute the topvisible list */
void ResetTopVisible(/*IN/UPDATE*/NJ_t *NJ,
		     int nActive,
		     /*IN/OUT*/top_hits_t *tophits) {
  besthit_t *visibleSorted = mymalloc(sizeof(besthit_t)*nActive);
  int nVisible = 0;		/* #entries in visibleSorted */
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    /* skip joins involving stale nodes */
    if (NJ->parent[iNode] >= 0)
      continue;
    besthit_t v;
    if (GetVisible(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, iNode, /*OUT*/&v)) {
      assert(nVisible < nActive);
      visibleSorted[nVisible++] = v;
    }
  }
  assert(nVisible > 0);
    
  qsort(/*IN/OUT*/visibleSorted,nVisible,sizeof(besthit_t),CompareHitsByCriterion);
    
  /* Only keep the top m items, and try to avoid duplicating i->j with j->i
     Note that visible(i) -> j does not necessarily imply visible(j) -> i,
     so we store what the pairing was (or -1 for not used yet)
   */
  int *inTopVisible = malloc(sizeof(int) * NJ->maxnodes);
  int i;
  for (i = 0; i < NJ->maxnodes; i++)
    inTopVisible[i] = -1;

  if (verbose > 2)
    fprintf(stderr, "top-hit search: nActive %d nVisible %d considering up to %d items\n",
	    nActive, nVisible, tophits->m);

  /* save the sorted indices in topvisible */
  int iSave = 0;
  for (i = 0; i < nVisible && iSave < tophits->nTopVisible; i++) {
    besthit_t *v = &visibleSorted[i];
    if (inTopVisible[v->i] != v->j) { /* not seen already */
      tophits->topvisible[iSave++] = v->i;
      inTopVisible[v->i] = v->j;
      inTopVisible[v->j] = v->i;
    }
  }
  while(iSave < tophits->nTopVisible)
    tophits->topvisible[iSave++] = -1;
  myfree(visibleSorted, sizeof(besthit_t)*nActive);
  myfree(inTopVisible, sizeof(int) * NJ->maxnodes);
  tophits->topvisibleAge = 0;
  if (verbose > 2) {
    fprintf(stderr, "Reset TopVisible: ");
    for (i = 0; i < tophits->nTopVisible; i++) {
      int iNode = tophits->topvisible[i];
      if (iNode < 0)
	break;
      fprintf(stderr, " %d=>%d", iNode, tophits->visible[iNode].j);
    }
    fprintf(stderr,"\n");
  }
}

int NGaps(/*IN*/NJ_t *NJ, int iNode) {
  assert(iNode < NJ->nSeq);
  int nGaps = 0;
  int p;
  for(p=0; p<NJ->nPos; p++) {
    if (NJ->profiles[iNode]->codes[p] == NOCODE)
      nGaps++;
  }
  return(nGaps);
}

char *OpenMPString(void) {
#ifdef OPENMP
  static char buf[100];
  sprintf(buf, ", OpenMP (%d threads)", omp_get_max_threads());
  return(buf);
#else
  return("");
#endif
}

/******************************************************************************/
/* Minimization of a 1-dimensional function by Brent's method (Numerical Recipes)            
 * Borrowed from Tree-Puzzle 5.1 util.c under GPL
 * Modified by M.N.P to pass in the accessory data for the optimization function,
 * to use 2x bounds around the starting guess and expand them if necessary,
 * and to use both a fractional and an absolute tolerance
 */

#define ITMAX 100
#define CGOLD 0.3819660
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Brents method in one dimension */
double brent(double ax, double bx, double cx, double (*f)(double, void *), void *data,
	     double ftol, double atol,
	     double *foptx, double *f2optx, double fax, double fbx, double fcx)
{
	int iter;
	double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if (fax < fcx) {
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else {
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;	
	}
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=ftol*fabs(x);
		tol2=2.0*(tol1+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))
		    || fabs(a-b) < atol) {
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
				(v*v*xw + x*x*wv + w*w*vx);
			return x;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u,data);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
		(v*v*xw + x*x*wv + w*w*vx);
	return x;
} /* brent */
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN

/* one-dimensional minimization - as input a lower and an upper limit and a trial
  value for the minimum is needed: xmin < xguess < xmax
  the function and a fractional tolerance has to be specified
  onedimenmin returns the optimal x value and the value of the function
  and its second derivative at this point
  */
double onedimenmin(double xmin, double xguess, double xmax, double (*f)(double,void*), void *data,
		   double ftol, double atol,
		   /*OUT*/double *fx, /*OUT*/double *f2x)
{
	double optx, ax, bx, cx, fa, fb, fc;
		
	/* first attempt to bracketize minimum */
	if (xguess == xmin) {
	  ax = xmin;
	  bx = 2.0*xguess;
	  cx = 10.0*xguess;
	} else if (xguess <= 2.0 * xmin) {
	  ax = xmin;
	  bx = xguess;
	  cx = 5.0*xguess;
	} else {
	  ax = 0.5*xguess;
	  bx = xguess;
	  cx = 2.0*xguess;
	}
	if (cx > xmax)
	  cx = xmax;
	if (bx >= cx)
	  bx = 0.5*(ax+cx);
	if (verbose > 4)
	  fprintf(stderr, "onedimenmin lo %.4f guess %.4f hi %.4f range %.4f %.4f\n",
		  ax, bx, cx, xmin, xmax);
	/* ideally this range includes the true minimum, i.e.,
	   fb < fa and fb < fc
	   if not, we gradually expand the boundaries until it does,
	   or we near the boundary of the allowed range and use that
	*/
	fa = (*f)(ax,data);
	fb = (*f)(bx,data);
	fc = (*f)(cx,data);
	while(fa < fb && ax > xmin) {
	  ax = (ax+xmin)/2.0;
	  if (ax < 2.0*xmin)	/* give up on shrinking the region */
	    ax = xmin;
	  fa = (*f)(ax,data);
	}
	while(fc < fb && cx < xmax) {
	  cx = (cx+xmax)/2.0;
	  if (cx > xmax * 0.95)
	    cx = xmax;
	  fc = (*f)(cx,data);
	}
	optx = brent(ax, bx, cx, f, data, ftol, atol, fx, f2x, fa, fb, fc);

	if (verbose > 4)
	  fprintf(stderr, "onedimenmin reaches optimum f(%.4f) = %.4f f2x %.4f\n", optx, *fx, *f2x);
	return optx; /* return optimal x */
} /* onedimenmin */

/* Numerical code for the gamma distribution is modified from the PhyML 3 code
   (GNU public license) of Stephane Guindon
*/

double LnGamma (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;
   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-(double)log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*(double)log(x) - x + .918938533204673
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;
}

double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=(double)exp(p*(double)log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}

double PGamma(double x, double alpha)
{
  /* scale = 1/alpha */
  return IncompleteGamma(x*alpha,alpha,LnGamma(alpha));
}

/* helper function to subtract timval structures */
/* Subtract the `struct timeval' values X and Y,
        storing the result in RESULT.
        Return 1 if the difference is negative, otherwise 0.  */
int     timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
  
  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
  
  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

double clockDiff(/*IN*/struct timeval *clock_start) {
  struct timeval time_now, elapsed;
  gettimeofday(/*OUT*/&time_now,NULL);
  timeval_subtract(/*OUT*/&elapsed,/*IN*/&time_now,/*IN*/clock_start);
  return(elapsed.tv_sec + elapsed.tv_usec*1e-6);
}


hashstrings_t *MakeHashtable(char **strings, int nStrings) {
  hashstrings_t *hash = (hashstrings_t*)mymalloc(sizeof(hashstrings_t));
  hash->nBuckets = 8*nStrings;
  hash->buckets = (hashbucket_t*)mymalloc(sizeof(hashbucket_t) * hash->nBuckets);
  int i;
  for (i=0; i < hash->nBuckets; i++) {
    hash->buckets[i].string = NULL;
    hash->buckets[i].nCount = 0;
    hash->buckets[i].first = -1;
  }
  for (i=0; i < nStrings; i++) {
    hashiterator_t hi = FindMatch(hash, strings[i]);
    if (hash->buckets[hi].string == NULL) {
      /* save a unique entry */
      assert(hash->buckets[hi].nCount == 0);
      hash->buckets[hi].string = strings[i];
      hash->buckets[hi].nCount = 1;
      hash->buckets[hi].first = i;
    } else {
      /* record a duplicate entry */
      assert(hash->buckets[hi].string != NULL);
      assert(strcmp(hash->buckets[hi].string, strings[i]) == 0);
      assert(hash->buckets[hi].first >= 0);
      hash->buckets[hi].nCount++;
    }
  }
  return(hash);
}

hashstrings_t *FreeHashtable(hashstrings_t* hash) {
  if (hash != NULL) {
    myfree(hash->buckets, sizeof(hashbucket_t) * hash->nBuckets);
    myfree(hash, sizeof(hashstrings_t));
  }
  return(NULL);
}

#define MAXADLER 65521
hashiterator_t FindMatch(hashstrings_t *hash, char *string) {
  /* Adler-32 checksum */
  unsigned int hashA = 1;
  unsigned int hashB = 0;
  char *p;
  for (p = string; *p != '\0'; p++) {
    hashA = ((unsigned int)*p + hashA);
    hashB = hashA+hashB;
  }
  hashA %= MAXADLER;
  hashB %= MAXADLER;
  hashiterator_t hi = (hashB*65536+hashA) % hash->nBuckets;
  while(hash->buckets[hi].string != NULL
	&& strcmp(hash->buckets[hi].string, string) != 0) {
    hi++;
    if (hi >= hash->nBuckets)
      hi = 0;
  }
  return(hi);
}

char *GetHashString(hashstrings_t *hash, hashiterator_t hi) {
  return(hash->buckets[hi].string);
}

int HashCount(hashstrings_t *hash, hashiterator_t hi) {
  return(hash->buckets[hi].nCount);
}

int HashFirst(hashstrings_t *hash, hashiterator_t hi) {
  return(hash->buckets[hi].first);
}

uniquify_t *UniquifyAln(alignment_t *aln) {
    int nUniqueSeq = 0;
    char **uniqueSeq = (char**)mymalloc(aln->nSeq * sizeof(char*)); /* iUnique -> seq */
    int *uniqueFirst = (int*)mymalloc(aln->nSeq * sizeof(int)); /* iUnique -> iFirst in aln */
    int *alnNext = (int*)mymalloc(aln->nSeq * sizeof(int)); /* i in aln -> next, or -1 */
    int *alnToUniq = (int*)mymalloc(aln->nSeq * sizeof(int)); /* i in aln -> iUnique; many -> -1 */

    int i;
    for (i = 0; i < aln->nSeq; i++) {
      uniqueSeq[i] = NULL;
      uniqueFirst[i] = -1;
      alnNext[i] = -1;
      alnToUniq[i] = -1;
    }
    hashstrings_t *hashseqs = MakeHashtable(aln->seqs, aln->nSeq);
    for (i=0; i<aln->nSeq; i++) {
      hashiterator_t hi = FindMatch(hashseqs,aln->seqs[i]);
      int first = HashFirst(hashseqs,hi);
      if (first == i) {
	uniqueSeq[nUniqueSeq] = aln->seqs[i];
	uniqueFirst[nUniqueSeq] = i;
	alnToUniq[i] = nUniqueSeq;
	nUniqueSeq++;
      } else {
	int last = first;
	while (alnNext[last] != -1)
	  last = alnNext[last];
	assert(last>=0);
	alnNext[last] = i;
	assert(alnToUniq[last] >= 0 && alnToUniq[last] < nUniqueSeq);
	alnToUniq[i] = alnToUniq[last];
      }
    }
    assert(nUniqueSeq>0);
    hashseqs = FreeHashtable(hashseqs);

    uniquify_t *uniquify = (uniquify_t*)mymalloc(sizeof(uniquify_t));
    uniquify->nSeq = aln->nSeq;
    uniquify->nUnique = nUniqueSeq;
    uniquify->uniqueFirst = uniqueFirst;
    uniquify->alnNext = alnNext;
    uniquify->alnToUniq = alnToUniq;
    uniquify->uniqueSeq = uniqueSeq;
    return(uniquify);
}

uniquify_t *FreeUniquify(uniquify_t *unique) {
  if (unique != NULL) {
    myfree(unique->uniqueFirst, sizeof(int)*unique->nSeq);
    myfree(unique->alnNext, sizeof(int)*unique->nSeq);
    myfree(unique->alnToUniq, sizeof(int)*unique->nSeq);
    myfree(unique->uniqueSeq, sizeof(char*)*unique->nSeq);
    myfree(unique,sizeof(uniquify_t));
    unique = NULL;
  }
  return(unique);
}

traversal_t InitTraversal(NJ_t *NJ) {
  traversal_t worked = (bool*)mymalloc(sizeof(bool)*NJ->maxnodes);
  int i;
  for (i=0; i<NJ->maxnodes; i++)
    worked[i] = false;
  return(worked);
}

void SkipTraversalInto(int node, /*IN/OUT*/traversal_t traversal) {
  traversal[node] = true;
}

int TraversePostorder(int node, NJ_t *NJ, /*IN/OUT*/traversal_t traversal,
		      /*OPTIONAL OUT*/bool *pUp) {
  if (pUp)
    *pUp = false;
  while(1) {
    assert(node >= 0);

    /* move to a child if possible */
    bool found = false;
    int iChild;
    for (iChild=0; iChild < NJ->child[node].nChild; iChild++) {
      int child = NJ->child[node].child[iChild];
      if (!traversal[child]) {
	node = child;
	found = true;
	break;
      }
    }
    if (found)
      continue; /* keep moving down */
    if (!traversal[node]) {
      traversal[node] = true;
      return(node);
    }
    /* If we've already done this node, need to move up */
    if (node == NJ->root)
      return(-1); /* nowhere to go -- done traversing */
    node = NJ->parent[node];
    /* If we go up to someplace that was already marked as visited, this is due
       to a change in topology, so return it marked as "up" */
    if (pUp && traversal[node]) {
      *pUp = true;
      return(node);
    }
  }
}

traversal_t FreeTraversal(traversal_t traversal, NJ_t *NJ) {
  myfree(traversal, sizeof(bool)*NJ->maxnodes);
  return(NULL);
}

profile_t **UpProfiles(NJ_t *NJ) {
  profile_t **upProfiles = (profile_t**)mymalloc(sizeof(profile_t*)*NJ->maxnodes);
  int i;
  for (i=0; i<NJ->maxnodes; i++) upProfiles[i] = NULL;
  return(upProfiles);
}

profile_t *GetUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int outnode, bool useML) {
  assert(outnode != NJ->root && outnode >= NJ->nSeq); /* not for root or leaves */
  if (upProfiles[outnode] != NULL)
    return(upProfiles[outnode]);

  int depth;
  int *pathToRoot = PathToRoot(NJ, outnode, /*OUT*/&depth);
  int i;
  /* depth-1 is root */
  for (i = depth-2; i>=0; i--) {
    int node = pathToRoot[i];

    if (upProfiles[node] == NULL) {
      /* Note -- SetupABCD may call GetUpProfile, but it should do it farther
	 up in the path to the root
      */
      profile_t *profiles[4];
      int nodeABCD[4];
      SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, useML);
      if (useML) {
	/* If node is a child of root, then the 4th profile is of the 2nd root-sibling of node
	   Otherwise, the 4th profile is the up-profile of the parent of node, and that
	   is the branch-length we need
	 */
	double lenC = NJ->branchlength[nodeABCD[2]];
	double lenD = NJ->branchlength[nodeABCD[3]];
	if (verbose > 3) {
	  fprintf(stderr, "Computing UpProfile for node %d with lenC %.4f lenD %.4f pair-loglk %.3f\n",
		  node, lenC, lenD,
		  PairLogLk(profiles[2],profiles[3],lenC+lenD,NJ->nPos,NJ->transmat,&NJ->rates, /*site_lk*/NULL));
	  PrintNJInternal(stderr, NJ, /*useLen*/true);
	}
	upProfiles[node] = PosteriorProfile(/*C*/profiles[2], /*D*/profiles[3],
					    lenC, lenD,
					    NJ->transmat, &NJ->rates, NJ->nPos, NJ->nConstraints);
      } else {
	profile_t *profilesCDAB[4] = { profiles[2], profiles[3], profiles[0], profiles[1] };
	double weight = QuartetWeight(profilesCDAB, NJ->distance_matrix, NJ->nPos);
	if (verbose>3)
	  fprintf(stderr, "Compute upprofile of %d from %d and parents (vs. children %d %d) with weight %.3f\n",
		  node, nodeABCD[2], nodeABCD[0], nodeABCD[1], weight);
	upProfiles[node] = AverageProfile(profiles[2], profiles[3],
					  NJ->nPos, NJ->nConstraints,
					  NJ->distance_matrix,
					  weight);
      }
    }
  }
  FreePath(pathToRoot,NJ);
  assert(upProfiles[outnode] != NULL);
  return(upProfiles[outnode]);
}

profile_t *DeleteUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node) {
  assert(node>=0 && node < NJ->maxnodes);
  if (upProfiles[node] != NULL)
    upProfiles[node] = FreeProfile(upProfiles[node], NJ->nPos, NJ->nConstraints); /* returns NULL */
  return(NULL);
}

profile_t **FreeUpProfiles(profile_t **upProfiles, NJ_t *NJ) {
  int i;
  int nUsed = 0;
  for (i=0; i < NJ->maxnodes; i++) {
    if (upProfiles[i] != NULL)
      nUsed++;
    DeleteUpProfile(upProfiles, NJ, i);
  }
  myfree(upProfiles, sizeof(profile_t*)*NJ->maxnodes);
  if (verbose >= 3)
    fprintf(stderr,"FreeUpProfiles -- freed %d\n", nUsed);
  return(NULL);
}

int *PathToRoot(NJ_t *NJ, int node, /*OUT*/int *outDepth) {
  int *pathToRoot = (int*)mymalloc(sizeof(int)*NJ->maxnodes);
  int depth = 0;
  int ancestor = node;
  while(ancestor >= 0) {
    pathToRoot[depth] = ancestor;
    ancestor = NJ->parent[ancestor];
    depth++;
  }
  *outDepth = depth;
  return(pathToRoot);
}

int *FreePath(int *path, NJ_t *NJ) {
  myfree(path, sizeof(int)*NJ->maxnodes);
  return(NULL);
}

transition_matrix_t *CreateGTR(double *r/*ac ag at cg ct gt*/, double *f/*acgt*/) {
  double matrix[MAXCODES][MAXCODES];
  assert(nCodes==4);
  int i, j;
  /* Place rates onto a symmetric matrix, but correct by f(target), so that
     stationary distribution f[] is maintained
     Leave diagonals as 0 (CreateTransitionMatrix will fix them)
  */
  int imat = 0;
  for (i = 0; i < nCodes; i++) {
    matrix[i][i] = 0;
    for (j = i+1; j < nCodes; j++) {
      double rate = r[imat++];
      assert(rate > 0);
      /* Want t(matrix) * f to be 0 */
      matrix[i][j] = rate * f[i];
      matrix[j][i] = rate * f[j];
    }
  }
  /* Compute average mutation rate */
  double total_rate = 0;
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      total_rate += f[i] * matrix[i][j];
  assert(total_rate > 1e-6);
  double inv = 1.0/total_rate;
  for (i = 0; i < nCodes; i++)
    for (j = 0; j < nCodes; j++)
      matrix[i][j] *= inv;
  return(CreateTransitionMatrix(matrix,f));
}
