#include "include_stuff.h"
#include "datastructs.h"
#include "hyper_parameters.h"
#include "support_functions.h"

int FastTree(FastTreeOptions_t opt) {

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
   		for (i=0; i < opt.argc; i++)
    		fprintf(fpLog, " %s", opt.argv[i]);
    	fprintf(fpLog,"\n");
    	fflush(fpLog);
  	}


	FILE *fps[2] = {NULL,NULL};
	int nFPs = 0;
	if (verbose)
		fps[nFPs++] = stderr;
	if (fpLog != NULL)
		fps[nFPs++] = fpLog;
  
	if (!make_matrix) {
		ReportSetting(opt, nFPs, fps);
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
      		MakeMatrix(aln, distance_matrix);
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
			ReadConstraints(&constraints, &uniqConstraints, &hashnames, &constraintsFile, &fpConstraints, bQuote, iAln, &unique);
			

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

				PrintStats(nFPs, fps, clock_start, NJ, aln, splitcount, uniqConstraints, MLnniToDo, MLlen, fpLog, nniToDo, spr);
			
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
