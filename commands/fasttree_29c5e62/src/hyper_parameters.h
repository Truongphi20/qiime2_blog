#ifndef HYPER_PARAMETERS_H
#define HYPER_PARAMETERS_H

#include "include_stuff.h"

/* Global variables */
/* Options */
int verbose = 1;
int showProgress = 1;
int slow = 0;
int fastest = 0;
bool useTopHits2nd = false;	/* use the second-level top hits heuristic? */
int bionj = 0;
double tophitsMult = 1.0;	/* 0 means compare nodes to all other nodes */
double tophitsClose = -1.0;	/* Parameter for how close is close; also used as a coverage req. */
double topvisibleMult = 1.5;	/* nTopVisible = m * topvisibleMult; 1 or 2 did not make much difference
				   in either running time or accuracy so I chose a compromise. */

double tophitsRefresh = 0.8;	/* Refresh if fraction of top-hit-length drops to this */
double tophits2Mult = 1.0;	/* Second-level top heuristic -- only with -fastest */
int tophits2Safety = 3;		/* Safety factor for second level of top-hits heuristic */
double tophits2Refresh = 0.6;	/* Refresh 2nd-level top hits if drops down to this fraction of length */

double staleOutLimit = 0.01;	/* nActive changes by at most this amount before we recompute 
				   an out-distance. (Only applies if using the top-hits heuristic) */
double fResetOutProfile = 0.02;	/* Recompute out profile from scratch if nActive has changed
				   by more than this proportion, and */
int nResetOutProfile = 200;	/* nActive has also changed more than this amount */
int nCodes=20;			/* 20 if protein, 4 if nucleotide */
bool useMatrix=true;		/* If false, use %different as the uncorrected distance */
bool logdist = true;		/* If true, do a log-correction (scoredist-like or Jukes-Cantor)
				   but only during NNIs and support values, not during neighbor-joining */
double pseudoWeight = 0.0;      /* The weight of pseudocounts to avoid artificial long branches when
				   nearby sequences in the tree have little or no overlap
				   (off by default). The prior distance is based on
				   all overlapping positions among the quartet or triplet under
				   consideration. The log correction takes place after the
				   pseudocount is used. */
double constraintWeight = 100.0;/* Cost of violation of a topological constraint in evolutionary distance
				   or likelihood */
double MEMinDelta = 1.0e-4;	/* Changes of less than this in tree-length are discounted for
				   purposes of identifying fixed subtrees */
bool fastNNI = true;
bool gammaLogLk = false;	/* compute gamma likelihood without reoptimizing branch lengths? */

/* Maximum likelihood options and constants */
/* These are used to rescale likelihood values and avoid taking a logarithm at each position */
const double LkUnderflow = 1.0e-4;
const double LkUnderflowInv = 1.0e4;
const double LogLkUnderflow = 9.21034037197618; /* -log(LkUnderflowInv) */
const double Log2 = 0.693147180559945;
/* These are used to limit the optimization of branch lengths.
   Also very short branch lengths can create numerical problems.
   In version 2.1.7, the minimum branch lengths (MLMinBranchLength and MLMinRelBranchLength)
   were increased to prevent numerical problems in rare cases.
   In version 2.1.8, to provide useful branch lengths for genome-wide alignments,
   the minimum branch lengths were dramatically decreased if USE_DOUBLE is defined.
*/
#ifndef USE_DOUBLE
const double MLMinBranchLengthTolerance = 1.0e-4; /* absolute tolerance for optimizing branch lengths */
const double MLFTolBranchLength = 0.001; /* fractional tolerance for optimizing branch lengths */
const double MLMinBranchLength = 5.0e-4; /* minimum value for branch length */
const double MLMinRelBranchLength = 2.5e-4; /* minimum of rate * length */
const double fPostTotalTolerance = 1.0e-10; /* posterior vector must sum to at least this before rescaling */
#else
const double MLMinBranchLengthTolerance = 1.0e-9;
const double MLFTolBranchLength = 0.001;
const double MLMinBranchLength = 5.0e-9;
const double MLMinRelBranchLength = 2.5e-9;
const double fPostTotalTolerance = 1.0e-20;
#endif

int mlAccuracy = 1;		/* Rounds of optimization of branch lengths; 1 means do 2nd round only if close */
double closeLogLkLimit = 5.0;	/* If partial optimization of an NNI looks like it would decrease the log likelihood
				   by this much or more then do not optimize it further */
double treeLogLkDelta = 0.1;	/* Give up if tree log-lk changes by less than this; NNIs that change
				   likelihood by less than this also are considered unimportant
				   by some heuristics */
bool exactML = true;		/* Exact or approximate posterior distributions for a.a.s */
double approxMLminf = 0.95;	/* Only try to approximate posterior distributions if max. value is at least this high */
double approxMLminratio = 2/3.0;/* Ratio of approximated/true posterior values must be at least this high */
double approxMLnearT = 0.2;	/* 2nd component of near-constant posterior distribution uses this time scale */
const int nDefaultRateCats = 20;

/* Performance and memory usage */
long profileOps = 0;		/* Full profile-based distance operations */
long outprofileOps = 0;		/* How many of profileOps are comparisons to outprofile */
long seqOps = 0;		/* Faster leaf-based distance operations */
long profileAvgOps = 0;		/* Number of profile-average steps */
long nHillBetter = 0;		/* Number of hill-climbing steps */
long nCloseUsed = 0;		/* Number of "close" neighbors we avoid full search for */
long nClose2Used = 0;		/* Number of "close" neighbors we use 2nd-level top hits for */
long nRefreshTopHits = 0;	/* Number of full-blown searches (interior nodes) */
long nVisibleUpdate = 0;		/* Number of updates of the visible set */
long nNNI = 0;			/* Number of NNI changes performed */
long nSPR = 0;			/* Number of SPR changes performed */
long nML_NNI = 0;		/* Number of max-lik. NNI changes performed */
long nSuboptimalSplits = 0;	/* # of splits that are rejected given final tree (during bootstrap) */
long nSuboptimalConstrained = 0; /* Bad splits that are due to constraints */
long nConstraintViolations = 0;	/* Number of constraint violations */
long nProfileFreqAlloc = 0;
long nProfileFreqAvoid = 0;
long szAllAlloc = 0;
long mymallocUsed = 0;		/* useful allocations by mymalloc */
long maxmallocHeap = 0;		/* Maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */
long nLkCompute = 0;		/* # of likelihood computations for pairs of probability vectors */
long nPosteriorCompute = 0;	/* # of computations of posterior probabilities */
long nAAPosteriorExact = 0;	/* # of times compute exact AA posterior */
long nAAPosteriorRough = 0;	/* # of times use rough approximation */
long nStarTests = 0;		/* # of times we use star test to avoid testing an NNI */

/* Protein character set */
unsigned char *codesStringAA = (unsigned char*) "ARNDCQEGHILKMFPSTWYV";
unsigned char *codesStringNT = (unsigned char*) "ACGT";
unsigned char *codesString = NULL;

#endif