#ifndef HYPER_PARAMETERS_H
#define HYPER_PARAMETERS_H

#include "include_stuff.h"
#include "datastructs.h"

/* Global variables */
/* Options */
extern int verbose;
extern int showProgress;
extern int slow;
extern int fastest;
extern bool useTopHits2nd;	/* use the second-level top hits heuristic? */
extern int bionj;
extern double tophitsMult;	/* 0 means compare nodes to all other nodes */
extern double tophitsClose;	/* Parameter for how close is close; also used as a coverage req. */
extern double topvisibleMult;	/* nTopVisible = m * topvisibleMult; 1 or 2 did not make much difference
				   in either running time or accuracy so I chose a compromise. */

extern double tophitsRefresh;	/* Refresh if fraction of top-hit-length drops to this */
extern double tophits2Mult;	/* Second-level top heuristic -- only with -fastest */
extern int tophits2Safety;		/* Safety factor for second level of top-hits heuristic */
extern double tophits2Refresh;	/* Refresh 2nd-level top hits if drops down to this fraction of length */

extern double staleOutLimit;	/* nActive changes by at most this amount before we recompute 
				   an out-distance. (Only applies if using the top-hits heuristic) */
extern double fResetOutProfile;	/* Recompute out profile from scratch if nActive has changed
				   by more than this proportion, and */
extern int nResetOutProfile;	/* nActive has also changed more than this amount */
extern int nCodes;			/* 20 if protein, 4 if nucleotide */
extern bool useMatrix;		/* If false, use %different as the uncorrected distance */
extern bool logdist;		/* If true, do a log-correction (scoredist-like or Jukes-Cantor)
				   but only during NNIs and support values, not during neighbor-joining */
extern double pseudoWeight;      /* The weight of pseudocounts to avoid artificial long branches when
				   nearby sequences in the tree have little or no overlap
				   (off by default). The prior distance is based on
				   all overlapping positions among the quartet or triplet under
				   consideration. The log correction takes place after the
				   pseudocount is used. */
extern double constraintWeight;/* Cost of violation of a topological constraint in evolutionary distance
				   or likelihood */
extern double MEMinDelta;	/* Changes of less than this in tree-length are discounted for
				   purposes of identifying fixed subtrees */
extern bool fastNNI;
extern bool gammaLogLk;	/* compute gamma likelihood without reoptimizing branch lengths? */

/* Maximum likelihood options and constants */
/* These are used to rescale likelihood values and avoid taking a logarithm at each position */
extern const double LkUnderflow;
extern const double LkUnderflowInv;
extern const double LogLkUnderflow;
extern const double Log2;

extern const double MLMinBranchLengthTolerance;
extern const double MLFTolBranchLength;
extern const double MLMinBranchLength;
extern const double MLMinRelBranchLength;
extern const double fPostTotalTolerance;

extern int mlAccuracy;		/* Rounds of optimization of branch lengths; 1 means do 2nd round only if close */
extern double closeLogLkLimit;	/* If partial optimization of an NNI looks like it would decrease the log likelihood
				   by this much or more then do not optimize it further */
extern double treeLogLkDelta;	/* Give up if tree log-lk changes by less than this; NNIs that change
				   likelihood by less than this also are considered unimportant
				   by some heuristics */
extern bool exactML;		/* Exact or approximate posterior distributions for a.a.s */
extern double approxMLminf;	/* Only try to approximate posterior distributions if max. value is at least this high */
extern double approxMLminratio;/* Ratio of approximated/true posterior values must be at least this high */
extern double approxMLnearT;	/* 2nd component of near-constant posterior distribution uses this time scale */
extern const int nDefaultRateCats;

/* Performance and memory usage */
extern long profileOps;		/* Full profile-based distance operations */
extern long outprofileOps;		/* How many of profileOps are comparisons to outprofile */
extern long seqOps;		/* Faster leaf-based distance operations */
extern long profileAvgOps;		/* Number of profile-average steps */
extern long nHillBetter;		/* Number of hill-climbing steps */
extern long nCloseUsed;		/* Number of "close" neighbors we avoid full search for */
extern long nClose2Used;		/* Number of "close" neighbors we use 2nd-level top hits for */
extern long nRefreshTopHits;	/* Number of full-blown searches (interior nodes) */
extern long nVisibleUpdate;		/* Number of updates of the visible set */
extern long nNNI;			/* Number of NNI changes performed */
extern long nSPR;			/* Number of SPR changes performed */
extern long nML_NNI;		/* Number of max-lik. NNI changes performed */
extern long nSuboptimalSplits;	/* # of splits that are rejected given final tree (during bootstrap) */
extern long nSuboptimalConstrained; /* Bad splits that are due to constraints */
extern long nConstraintViolations;	/* Number of constraint violations */
extern long nProfileFreqAlloc;
extern long nProfileFreqAvoid;
extern long szAllAlloc;
extern long mymallocUsed;		/* useful allocations by mymalloc */
extern long maxmallocHeap;		/* Maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */
extern long nLkCompute;		/* # of likelihood computations for pairs of probability vectors */
extern long nPosteriorCompute;		/* # of computations of posterior probabilities */
extern long nAAPosteriorExact;	/* # of times compute exact AA posterior */
extern long nAAPosteriorRough;	/* # of times use rough approximation */
extern long nStarTests;		/* # of times we use star test to avoid testing an NNI */

/* Protein character set */
extern unsigned char *codesStringAA;
extern unsigned char *codesStringNT;
extern unsigned char *codesString;

/* Matrix and stationary distributions */
extern distance_matrix_t matrixBLOSUM45;
extern double matrixJTT92[MAXCODES][MAXCODES];
extern double statJTT92[MAXCODES];
extern double matrixLG08[MAXCODES][MAXCODES];
extern double statLG08[MAXCODES];
extern double matrixWAG01[MAXCODES][MAXCODES];
extern double statWAG01[MAXCODES];

#endif
