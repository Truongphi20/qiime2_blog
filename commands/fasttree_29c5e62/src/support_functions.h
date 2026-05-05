#ifndef SUPPORT_FUNC_H
#define SUPPORT_FUNC_H

#include "datastructs.h"
#include "include_stuff.h"
#include "hyper_parameters.h"

int FastTree(FastTreeOptions_t opt);
void ReportSetting(FastTreeOptions_t opt, int nFPs, FILE *fps[2]);
void MakeMatrix(alignment_t *aln, distance_matrix_t *distance_matrix);
void ReadConstraints(
    alignment_t **constraints, 
    char ***uniqConstraints, 
    hashstrings_t **hashnames, 
    char **constraintsFile,
    FILE **fpConstraints,
    bool bQuote,
    int iAln,
    uniquify_t **unique 
);
void PrintStats(
    int nFPs,
    FILE *fps[2],
    struct timeval clock_start,
    NJ_t *NJ,
    alignment_t *aln,
    SplitCount_t splitcount,
    char **uniqConstraints,
    int MLnniToDo,
    bool MLlen,
    FILE *fpLog,
    int nniToDo,
    int spr
);

void InitOptions(FastTreeOptions_t *opt);
void ParseCommandLine(int argc, char **argv, FastTreeOptions_t *opt);

distance_matrix_t *ReadDistanceMatrix(char *prefix);
void SetupDistanceMatrix(/*IN/OUT*/distance_matrix_t *); /* set eigentot, codeFreq, gapFreq */
void ReadMatrix(char *filename, /*OUT*/numeric_t codes[MAXCODES][MAXCODES], bool check_codes);
void ReadVector(char *filename, /*OUT*/numeric_t codes[MAXCODES]);
alignment_t *ReadAlignment(/*READ*/FILE *fp, bool bQuote); /* Returns a list of strings (exits on failure) */
alignment_t *FreeAlignment(alignment_t *); /* returns NULL */
void FreeAlignmentSeqs(/*IN/OUT*/alignment_t *);

/* Takes as input the transpose of the matrix V, with i -> j
   This routine takes care of setting the diagonals
*/
transition_matrix_t *CreateTransitionMatrix(/*IN*/double matrix[MAXCODES][MAXCODES],
					    /*IN*/double stat[MAXCODES]);
transition_matrix_t *CreateGTR(double *gtrrates/*ac,ag,at,cg,ct,gt*/, double *gtrfreq/*ACGT*/);
transition_matrix_t *ReadAATransitionMatrix(/*IN*/char *filename);

/* For converting profiles from 1 rotation to another, or converts NULL to NULL */
distance_matrix_t *TransMatToDistanceMat(transition_matrix_t *transmat);

/* Allocates memory, initializes leaf profiles */
NJ_t *InitNJ(char **sequences, int nSeqs, int nPos,
	     /*IN OPTIONAL*/char **constraintSeqs, int nConstraints,
	     /*IN OPTIONAL*/distance_matrix_t *,
	     /*IN OPTIONAL*/transition_matrix_t *);

NJ_t *FreeNJ(NJ_t *NJ); /* returns NULL */
void FastNJ(/*IN/OUT*/NJ_t *NJ); /* Does the joins */
void ReliabilityNJ(/*IN/OUT*/NJ_t *NJ, int nBootstrap);	  /* Estimates the reliability of the joins */

/* nni_stats_t is meaningless for leaves and root, so all of those entries
   will just be high (for age) or 0 (for delta)
*/
typedef struct {
  int age;	    /* number of rounds since this node was modified by an NNI */
  int subtreeAge;   /* number of rounds since self or descendent had a significant improvement */
  double delta;	    /* improvement in score for this node (or 0 if no change) */
  double support;   /* improvement of score for self over better of alternatives */
} nni_stats_t;

/* One round of nearest-neighbor interchanges according to the
   minimum-evolution or approximate maximum-likelihood criterion.
   If doing maximum likelihood then this modifies the branch lengths.
   age is the # of rounds since a node was NNId
   Returns the # of topological changes performed
*/
int NNI(/*IN/OUT*/NJ_t *NJ, int iRound, int nRounds, bool useML,
	/*IN/OUT*/nni_stats_t *stats,
	/*OUT*/double *maxDeltaCriterion);
nni_stats_t *InitNNIStats(NJ_t *NJ);
nni_stats_t *FreeNNIStats(nni_stats_t *, NJ_t *NJ);	/* returns NULL */

/* One round of subtree-prune-regraft moves (minimum evolution) */
void SPR(/*IN/OUT*/NJ_t *NJ, int maxSPRLength, int iRound, int nRounds);

/* Recomputes all branch lengths by minimum evolution criterion*/
void UpdateBranchLengths(/*IN/OUT*/NJ_t *NJ);

/* Recomputes all branch lengths and, optionally, internal profiles */
double TreeLength(/*IN/OUT*/NJ_t *NJ, bool recomputeProfiles);

void TestSplitsMinEvo(NJ_t *NJ, /*OUT*/SplitCount_t *splitcount);

/* Sets SH-like support values if nBootstrap>0 */
void TestSplitsML(/*IN/OUT*/NJ_t *NJ, /*OUT*/SplitCount_t *splitcount, int nBootstrap);

/* Pick columns for resampling, stored as returned_vector[iBoot*nPos + j] */
int *ResampleColumns(int nPos, int nBootstrap);

/* Use out-profile and NJ->totdiam to recompute out-distance for node iNode
   Only does this computation if the out-distance is "stale" (nOutDistActive[iNode] != nActive)
   Note "IN/UPDATE" for NJ always means that we may update out-distances but otherwise
   make no changes.
 */
void SetOutDistance(/*IN/UPDATE*/NJ_t *NJ, int iNode, int nActive);

/* Always sets join->criterion; may update NJ->outDistance and NJ->nOutDistActive,
   assumes join's weight and distance are already set,
   and that the constraint penalty (if any) is included in the distance
*/
void SetCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join);

/* Computes weight and distance (which includes the constraint penalty)
   and then sets the criterion (maybe update out-distances)
*/
void SetDistCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join);

/* If join->i or join->j are inactive nodes, replaces them with their active ancestors.
   After doing this, if i == j, or either is -1, sets weight to 0 and dist and criterion to 1e20
      and returns false (not a valid join)
   Otherwise, if i or j changed, recomputes the distance and criterion.
   Note that if i and j are unchanged then the criterion could be stale
   If bUpdateDist is false, and i or j change, then it just sets dist to a negative number
*/
bool UpdateBestHit(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join,
		   bool bUpdateDist);

/* This recomputes the criterion, or returns false if the visible node
   is no longer active.
*/
bool GetVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/top_hits_t *tophits,
		int iNode, /*OUT*/besthit_t *visible);

int ActiveAncestor(/*IN*/NJ_t *NJ, int node);

/* Compute the constraint penalty for a join. This is added to the "distance"
   by SetCriterion */
int JoinConstraintPenalty(/*IN*/NJ_t *NJ, int node1, int node2);
int JoinConstraintPenaltyPiece(NJ_t *NJ, int node1, int node2, int iConstraint);

/* Helper function for computing the number of constraints violated by
   a split, represented as counts of on and off on each side */
int SplitConstraintPenalty(int nOn1, int nOff1, int nOn2, int nOff2);

/* Reports the (min. evo.) support for the (1,2) vs. (3,4) split
   col[iBoot*nPos+j] is column j for bootstrap iBoot
*/
double SplitSupport(profile_t *p1, profile_t *p2, profile_t *p3, profile_t *p4,
		    /*OPTIONAL*/distance_matrix_t *dmat,
		    int nPos,
		    int nBootstrap,
		    int *col);

/* Returns SH-like support given resampling spec. (in col) and site likelihods
   for the three quartets
*/
double SHSupport(int nPos, int nBoostrap, int *col, double loglk[3], double *site_likelihoods[3]);

profile_t *SeqToProfile(/*IN/OUT*/NJ_t *NJ,
			char *seq, int nPos,
			/*OPTIONAL*/char *constraintSeqs, int nConstraints,
			int iNode,
			unsigned long counts[256]);

/* ProfileDist and SeqDist only set the dist and weight fields
   If using an outprofile, use the second argument of ProfileDist
   for better performance.

   These produce uncorrected distances.
*/
void ProfileDist(profile_t *profile1, profile_t *profile2, int nPos,
		 /*OPTIONAL*/distance_matrix_t *distance_matrix,
		 /*OUT*/besthit_t *hit);
void SeqDist(unsigned char *codes1, unsigned char *codes2, int nPos,
	     /*OPTIONAL*/distance_matrix_t *distance_matrix,
	     /*OUT*/besthit_t *hit);

/* Computes all pairs of profile distances, applies pseudocounts
   if pseudoWeight > 0, and applies log-correction if logdist is true.
   The lower index is compared to the higher index, e.g. for profiles
   A,B,C,D the comparison will be as in quartet_pair_t
*/
typedef enum {qAB,qAC,qAD,qBC,qBD,qCD} quartet_pair_t;
void CorrectedPairDistances(profile_t **profiles, int nProfiles,
			    /*OPTIONAL*/distance_matrix_t *distance_matrix,
			    int nPos,
			    /*OUT*/double *distances);

/* output is indexed by nni_t
   To ensure good behavior while evaluating a subtree-prune-regraft move as a series
   of nearest-neighbor interchanges, this uses a distance-ish model of constraints,
   as given by PairConstraintDistance(), rather than
   counting the number of violated splits (which is what FastTree does
   during neighbor-joining).
   Thus, penalty values may well be >0 even if no constraints are violated, but the
   relative scores for the three NNIs will be correct.
 */
void QuartetConstraintPenalties(profile_t *profiles[4], int nConstraints, /*OUT*/double d[3]);

double PairConstraintDistance(int nOn1, int nOff1, int nOn2, int nOff2);

/* the split is consistent with the constraint if any of the profiles have no data
   or if three of the profiles have the same uniform value (all on or all off)
   or if AB|CD = 00|11 or 11|00 (all uniform)
 */
bool SplitViolatesConstraint(profile_t *profiles[4], int iConstraint);

/* If false, no values were set because this constraint was not relevant.
   output is for the 3 splits
*/
bool QuartetConstraintPenaltiesPiece(profile_t *profiles[4], int iConstraint, /*OUT*/double penalty[3]);

/* Apply Jukes-Cantor or scoredist-like log(1-d) transform
   to correct the distance for multiple substitutions.
*/
double LogCorrect(double distance);

/* AverageProfile is used to do a weighted combination of nodes
   when doing a join. If weight is negative, then the value is ignored and the profiles
   are averaged. The weight is *not* adjusted for the gap content of the nodes.
   Also, the weight does not affect the representation of the constraints
*/
profile_t *AverageProfile(profile_t *profile1, profile_t *profile2,
			  int nPos, int nConstraints,
			  distance_matrix_t *distance_matrix,
			  double weight1);

/* PosteriorProfile() is like AverageProfile() but it computes posterior probabilities
   rather than an average
*/
profile_t *PosteriorProfile(profile_t *profile1, profile_t *profile2,
			    double len1, double len2,
			    /*OPTIONAL*/transition_matrix_t *transmat,
			    rates_t *rates,
			    int nPos, int nConstraints);

/* Set a node's profile from its children.
   Deletes the previous profile if it exists
   Use -1.0 for a balanced join
   Fails unless the node has two children (e.g., no leaves or root)
*/
void SetProfile(/*IN/OUT*/NJ_t *NJ, int node, double weight1);

/* OutProfile does an unweighted combination of nodes to create the
   out-profile. It always sets code to NOCODE so that UpdateOutProfile
   can work.
*/
profile_t *OutProfile(profile_t **profiles, int nProfiles,
		      int nPos, int nConstraints,
		      distance_matrix_t *distance_matrix);

void UpdateOutProfile(/*UPDATE*/profile_t *out, profile_t *old1, profile_t *old2,
		      profile_t *new, int nActiveOld,
		      int nPos, int nConstraints,
		      distance_matrix_t *distance_matrix);

profile_t *NewProfile(int nPos, int nConstraints); /* returned has no vectors */
profile_t *FreeProfile(profile_t *profile, int nPos, int nConstraints); /* returns NULL */

void AllocRateCategories(/*IN/OUT*/rates_t *rates, int nRateCategories, int nPos);

/* f1 can be NULL if code1 != NOCODE, and similarly for f2
   Or, if (say) weight1 was 0, then can have code1==NOCODE *and* f1==NULL
   In that case, returns an arbitrary large number.
*/
double ProfileDistPiece(unsigned int code1, unsigned int code2,
			numeric_t *f1, numeric_t *f2, 
			/*OPTIONAL*/distance_matrix_t *dmat,
			/*OPTIONAL*/numeric_t *codeDist2);

/* Adds (or subtracts, if weight is negative) fIn/codeIn from fOut
   fOut is assumed to exist (as from an outprofile)
   do not call unless weight of input profile > 0
 */
void AddToFreq(/*IN/OUT*/numeric_t *fOut, double weight,
	       unsigned int codeIn, /*OPTIONAL*/numeric_t *fIn,
	       /*OPTIONAL*/distance_matrix_t *dmat);

/* Divide the vector (of length nCodes) by a constant
   so that the total (unrotated) frequency is 1.0 */
void NormalizeFreq(/*IN/OUT*/numeric_t *freq, distance_matrix_t *distance_matrix);

/* Allocate, if necessary, and recompute the codeDist*/
void SetCodeDist(/*IN/OUT*/profile_t *profile, int nPos, distance_matrix_t *dmat);

/* The allhits list contains the distances of the node to all other active nodes
   This is useful for the "reset" improvement to the visible set
   Note that the following routines do not handle the tophits heuristic
   and assume that out-distances are up to date.
*/
void SetBestHit(int node, NJ_t *NJ, int nActive,
		/*OUT*/besthit_t *bestjoin,
		/*OUT OPTIONAL*/besthit_t *allhits);
void ExhaustiveNJSearch(NJ_t *NJ, int nActive, /*OUT*/besthit_t *bestjoin);

/* Searches the visible set */
void FastNJSearch(NJ_t *NJ, int nActive, /*UPDATE*/besthit_t *visible, /*OUT*/besthit_t *bestjoin);

/* Subroutines for handling the tophits heuristic */

top_hits_t *InitTopHits(NJ_t *NJ, int m);
top_hits_t *FreeTopHits(top_hits_t *tophits); /* returns NULL */

/* Before we do any joins -- sets tophits and visible
   NJ may be modified by setting out-distances
 */
void SetAllLeafTopHits(/*IN/UPDATE*/NJ_t *NJ, /*IN/OUT*/top_hits_t *tophits);

/* 
   Find the best join to do. 
   Find best hit to do in O(N*log(N) + m*L*log(N)) time, by
   copying and sorting the visible list
   updating out-distances for the top (up to m) candidates
   selecting the best hit
   if !fastest then
      local hill-climbing for a better join,
      using best-hit lists only, and updating
      all out-distances in every best-hit list
*/
void TopHitNJSearch(/*IN/UPDATE*/NJ_t *NJ,
		    int nActive,
		    /*IN/OUT*/top_hits_t *tophits,
		    /*OUT*/besthit_t *bestjoin);

/* Returns the best hit within top hits
   NJ may be modified because it updates out-distances if they are too stale
   Does *not* update visible set
   Updates out-distances but does not reset or update visible set
   
*/
void GetBestFromTopHits(int iNode, /*IN/UPDATE*/NJ_t *NJ, int nActive,
			/*IN*/top_hits_t *tophits,
			/*OUT*/besthit_t *bestjoin);

/* visible set is modifiable so that we can reset it more globally when we do
   a "refresh", but we also set the visible set for newnode and do any
   "reset" updates too. And, we update many outdistances.

   Create a top hit list for the new node, either
   from children (if there are enough best hits left) or by a "refresh"
   Also set visible set for newnode
   Also update visible set for other nodes if we stumble across a "better" hit
 */
void TopHitJoin(int newnode,
		/*IN/UPDATE*/NJ_t *NJ, int nActive,
		/*IN/OUT*/top_hits_t *tophits);

/* Sort the input besthits by criterion
   and save the best nOut hits as a new array in top_hits_lists
   Does not update criterion or out-distances
   Ignores (silently removes) hit to self
   Saved list may be shorter than requested if there are insufficient entries
*/
void SortSaveBestHits(int iNode, /*IN/SORT*/besthit_t *besthits,
		      int nIn, int nOut,
		      /*IN/OUT*/top_hits_t *tophits);

/* Given candidate hits from one node, "transfer" them to another node:
   Stores them in a new place in the same order
   searches up to active nodes if hits involve non-active nodes
   If update flag is set, it also recomputes distance and criterion
   (and ensures that out-distances are updated); otherwise
   it sets dist to -1e20 and criterion to 1e20

 */
void TransferBestHits(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		      int iNode,
		      /*IN*/besthit_t *oldhits,
		      int nOldHits,
		      /*OUT*/besthit_t *newhits,
		      bool updateDistance);

/* Create best hit objects from 1 or more hits. Do not update out-distances or set criteria */
void HitsToBestHits(/*IN*/hit_t *hits, int nHits, int iNode, /*OUT*/besthit_t *newhits);
besthit_t HitToBestHit(int i, hit_t hit);

/* Given a set of besthit entries,
   look for improvements to the visible set of the j entries.
   Updates out-distances as it goes.
   Also replaces stale nodes with this node, because a join is usually
   how this happens (i.e. it does not need to walk up to ancestors).
   Note this calls UpdateTopVisible() on any change
*/
void UpdateVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		   /*IN*/besthit_t *tophitsNode,
		   int nTopHits,
		   /*IN/OUT*/top_hits_t *tophits);

/* Update the top-visible list to perhaps include this hit (O(sqrt(N)) time) */
void UpdateTopVisible(/*IN*/NJ_t * NJ, int nActive,
		      int iNode, /*IN*/hit_t *hit,
		      /*IN/OUT*/top_hits_t *tophits);

/* Recompute the top-visible subset of the visible set */
void ResetTopVisible(/*IN/UPDATE*/NJ_t *NJ,
		     int nActive,
		     /*IN/OUT*/top_hits_t *tophits);

/* Make a shorter list with only unique entries.
   Replaces any "dead" hits to nodes that have parents with their active ancestors
   and ignores any that become dead.
   Updates all criteria.
   Combined gets sorted by i & j
   The returned list is allocated to nCombined even though only *nUniqueOut entries are filled
*/
besthit_t *UniqueBestHits(/*IN/UPDATE*/NJ_t *NJ, int nActive,
			  /*IN/SORT*/besthit_t *combined, int nCombined,
			  /*OUT*/int *nUniqueOut);

nni_t ChooseNNI(profile_t *profiles[4],
		/*OPTIONAL*/distance_matrix_t *dmat,
		int nPos, int nConstraints,
		/*OUT*/double criteria[3]); /* The three internal branch lengths or log likelihoods*/

/* length[] is ordered as described by quartet_length_t, but after we do the swap
   of B with C (to give AC|BD) or B with D (to get AD|BC), if that is the returned choice
   bFast means do not consider NNIs if AB|CD is noticeably better than the star topology
   (as implemented by MLQuartetOptimize).
   If there are constraints, then the constraint penalty is included in criteria[]
*/
nni_t MLQuartetNNI(profile_t *profiles[4],
		   /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		   int nPos, int nConstraints,
		   /*OUT*/double criteria[3], /* The three potential quartet log-likelihoods */
		   /*IN/OUT*/numeric_t length[5],
		   bool bFast);

void OptimizeAllBranchLengths(/*IN/OUT*/NJ_t *NJ);
double TreeLogLk(/*IN*/NJ_t *NJ, /*OPTIONAL OUT*/double *site_loglk);
double MLQuartetLogLk(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
		      int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		      /*IN*/double branch_lengths[5],
		      /*OPTIONAL OUT*/double *site_likelihoods);

/* Given a topology and branch lengths, estimate rates & recompute profiles */
void SetMLRates(/*IN/OUT*/NJ_t *NJ, int nRateCategories);

/* Returns a set of nRateCategories potential rates; the caller must free it */
numeric_t *MLSiteRates(int nRateCategories);

/* returns site_loglk so that
   site_loglk[nPos*iRate + j] is the log likelihood of site j with rate iRate
   The caller must free it.
*/
double *MLSiteLikelihoodsByRate(/*IN*/NJ_t *NJ, /*IN*/numeric_t *rates, int nRateCategories);

typedef struct {
  double mult;			/* multiplier for the rates / divisor for the tree-length */
  double alpha;
  int nPos;
  int nRateCats;
  numeric_t *rates;
  double *site_loglk;
} siteratelk_t;

double GammaLogLk(/*IN*/siteratelk_t *s, /*OPTIONAL OUT*/double *gamma_loglk_sites);

/* Input site_loglk must be for each rate. Note that FastTree does not reoptimize
   the branch lengths under the Gamma model -- it optimizes the overall scale.
   Reports the gamma log likelihhod (and logs site likelihoods if fpLog is set),
   and reports the rescaling value.
*/
double RescaleGammaLogLk(int nPos, int nRateCats,
			/*IN*/numeric_t *rates, /*IN*/double *site_loglk,
			/*OPTIONAL*/FILE *fpLog);

/* P(value<=x) for the gamma distribution with shape parameter alpha and scale 1/alpha */
double PGamma(double x, double alpha);

/* Given a topology and branch lengths, optimize GTR rates and quickly reoptimize branch lengths
   If gtrfreq is NULL, then empirical frequencies are used
*/
void SetMLGtr(/*IN/OUT*/NJ_t *NJ, /*OPTIONAL IN*/double *gtrfreq, /*OPTIONAL WRITE*/FILE *fpLog);

/* P(A & B | len) = P(B | A, len) * P(A)
   If site_likelihoods is present, multiplies those values by the site likelihood at each point
   (Note it does not handle underflow)
 */
double PairLogLk(/*IN*/profile_t *p1, /*IN*/profile_t *p2, double length,
		 int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		 /*OPTIONAL IN/OUT*/double *site_likelihoods);

/* Branch lengths for 4-taxon tree ((A,B),C,D); I means internal */
typedef enum {LEN_A,LEN_B,LEN_C,LEN_D,LEN_I} quartet_length_t;

typedef struct {
  int nPos;
  transition_matrix_t *transmat;
  rates_t *rates;
  int nEval;			/* number of likelihood evaluations */
  /* The pair to optimize */
  profile_t *pair1;
  profile_t *pair2;
} quartet_opt_t;

double PairNegLogLk(double x, void *data); /* data must be a quartet_opt_t */

typedef struct {
  NJ_t *NJ;
  double freq[4];
  double rates[6];
  int iRate;			/* which rate to set x from */
  FILE *fpLog; /* OPTIONAL WRITE */
} gtr_opt_t;

/* Returns -log_likelihood for the tree with the given rates
   data must be a gtr_opt_t and x is used to set rate iRate
   Does not recompute profiles -- assumes that the caller will
*/
double GTRNegLogLk(double x, void *data);

/* Returns the resulting log likelihood. Optionally returns whether other
   topologies should be abandoned, based on the difference between AB|CD and
   the "star topology" (AB|CD with a branch length of MLMinBranchLength) exceeding
   closeLogLkLimit.
   If bStarTest is passed in, it only optimized the internal branch if
   the star test is true. Otherwise, it optimized all 5 branch lengths
   in turn.
 */
double MLQuartetOptimize(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
			 int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
			 /*IN/OUT*/double branch_lengths[5],
			 /*OPTIONAL OUT*/bool *pStarTest,
			 /*OPTIONAL OUT*/double *site_likelihoods);

/* Returns the resulting log likelihood */
double MLPairOptimize(profile_t *pA, profile_t *pB,
		      int nPos, /*OPTIONAL*/transition_matrix_t *transmat, rates_t *rates,
		      /*IN/OUT*/double *branch_length);

/* Returns the number of steps considered, with the actual steps in steps[]
   Modifies the tree by this chain of NNIs
*/
int FindSPRSteps(/*IN/OUT*/NJ_t *NJ, 
		 int node,
		 int parent,	/* sibling or parent of node to NNI to start the chain */
		 /*IN/OUT*/profile_t **upProfiles,
		 /*OUT*/spr_step_t *steps,
		 int maxSteps,
		 bool bFirstAC);

/* Undo a single NNI */
void UnwindSPRStep(/*IN/OUT*/NJ_t *NJ,
	       /*IN*/spr_step_t *step,
	       /*IN/OUT*/profile_t **upProfiles);


/* Update the profile of node and its ancestor, and delete nearby out-profiles */
void UpdateForNNI(/*IN/OUT*/NJ_t *NJ, int node, /*IN/OUT*/profile_t **upProfiles, bool useML);

/* Sets NJ->parent[newchild] and replaces oldchild with newchild
   in the list of children of parent
*/
void ReplaceChild(/*IN/OUT*/NJ_t *NJ, int parent, int oldchild, int newchild);

int CompareHitsByCriterion(const void *c1, const void *c2);
int CompareHitsByIJ(const void *c1, const void *c2);

int NGaps(NJ_t *NJ, int node);	/* only handles leaf sequences */

/* node is the parent of AB, sibling of C
   node cannot be root or a leaf
   If node is the child of root, then D is the other sibling of node,
   and the 4th profile is D's profile.
   Otherwise, D is the parent of node, and we use its upprofile
   Call this with profiles=NULL to get the nodes, without fetching or
   computing profiles
*/
void SetupABCD(NJ_t *NJ, int node,
	       /* the 4 profiles for ABCD; the last one is an upprofile */
	       /*OPTIONAL OUT*/profile_t *profiles[4], 
	       /*OPTIONAL IN/OUT*/profile_t **upProfiles,
	       /*OUT*/int nodeABCD[4],
	       bool useML);

int Sibling(NJ_t *NJ, int node); /* At root, no unique sibling so returns -1 */
void RootSiblings(NJ_t *NJ, int node, /*OUT*/int sibs[2]);

/* JC probability of nucleotide not changing, for each rate category */
double *PSameVector(double length, rates_t *rates);

/* JC probability of nucleotide not changing, for each rate category */
double *PDiffVector(double *pSame, rates_t *rates);

/* expeigen[iRate*nCodes + j] = exp(length * rate iRate * eigenvalue j) */
numeric_t *ExpEigenRates(double length, transition_matrix_t *transmat, rates_t *rates);

/* Print a progress report if more than 0.1 second has gone by since the progress report */
/* Format should include 0-4 %d references and no newlines */
void ProgressReport(char *format, int iArg1, int iArg2, int iArg3, int iArg4);
void LogTree(char *format, int round, /*OPTIONAL WRITE*/FILE *fp, NJ_t *NJ, char **names, uniquify_t *unique, bool bQuote);
void LogMLRates(/*OPTIONAL WRITE*/FILE *fpLog, NJ_t *NJ);

void *mymalloc(size_t sz);       /* Prints "Out of memory" and exits on failure */
void *myfree(void *, size_t sz); /* Always returns NULL */

/* One-dimensional minimization using brent's function, with
   a fractional and an absolute tolerance */
double onedimenmin(double xmin, double xguess, double xmax, double (*f)(double,void*), void *data,
		   double ftol, double atol,
		   /*OUT*/double *fx, /*OUT*/double *f2x);

double brent(double ax, double bx, double cx, double (*f)(double, void *), void *data,
	     double ftol, double atol,
	     double *foptx, double *f2optx, double fax, double fbx, double fcx);

/* Vector operations, either using SSE3 or not
   Code assumes that vectors are a multiple of 4 in size
*/
void vector_multiply(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, int n, /*OUT*/numeric_t *fOut);
numeric_t vector_multiply_sum(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, int n);
void vector_add_mult(/*IN/OUT*/numeric_t *f, /*IN*/numeric_t *add, numeric_t weight, int n);

/* multiply the transpose of a matrix by a vector */
void matrixt_by_vector4(/*IN*/numeric_t mat[4][MAXCODES], /*IN*/numeric_t vec[4], /*OUT*/numeric_t out[4]);

/* sum(f1*fBy)*sum(f2*fBy) */
numeric_t vector_dot_product_rot(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, /*IN*/numeric_t* fBy, int n);

/* sum(f1*f2*f3) */
numeric_t vector_multiply3_sum(/*IN*/numeric_t *f1, /*IN*/numeric_t *f2, /*IN*/numeric_t* f3, int n);

numeric_t vector_sum(/*IN*/numeric_t *f1, int n);
void vector_multiply_by(/*IN/OUT*/numeric_t *f, /*IN*/numeric_t fBy, int n);

double clockDiff(/*IN*/struct timeval *clock_start);
int timeval_subtract (/*OUT*/struct timeval *result, /*IN*/struct timeval *x, /*IN*/struct timeval *y);

char *OpenMPString(void);

void ran_start(long seed);
double knuth_rand();		/* Random number between 0 and 1 */
void tred2 (double *a, const int n, const int np, double *d, double *e);
double pythag(double a, double b);
void tqli(double *d, double *e, int n, int np, double *z);

/* Like mymalloc; duplicates the input (returns NULL if given NULL) */
void *mymemdup(void *data, size_t sz);
void *myrealloc(void *data, size_t szOld, size_t szNew, bool bCopy);

double pnorm(double z);		/* Probability(value <=z)  */

hashstrings_t *MakeHashtable(char **strings, int nStrings);
hashstrings_t *FreeHashtable(hashstrings_t* hash); /*returns NULL*/
hashiterator_t FindMatch(hashstrings_t *hash, char *string);

/* Return NULL if we have run out of values */
char *GetHashString(hashstrings_t *hash, hashiterator_t hi);
int HashCount(hashstrings_t *hash, hashiterator_t hi);
int HashFirst(hashstrings_t *hash, hashiterator_t hi);

void PrintNJ(/*WRITE*/FILE *, NJ_t *NJ, char **names, uniquify_t *unique, bool bShowSupport, bool bQuoteNames);

/* Print topology using node indices as node names */
void PrintNJInternal(/*WRITE*/FILE *, NJ_t *NJ, bool useLen);

uniquify_t *UniquifyAln(/*IN*/alignment_t *aln);
uniquify_t *FreeUniquify(uniquify_t *);	/* returns NULL */

/* Convert a constraint alignment to a list of sequences. The returned array is indexed
   by iUnique and points to values in the input alignment
*/
char **AlnToConstraints(alignment_t *constraints, uniquify_t *unique, hashstrings_t *hashnames);

/* ReadTree ignores non-unique leaves after the first instance.
   At the end, it prunes the tree to ignore empty children and it
   unroots the tree if necessary.
*/
void ReadTree(/*IN/OUT*/NJ_t *NJ,
	      /*IN*/uniquify_t *unique,
	      /*IN*/hashstrings_t *hashnames,
	      /*READ*/FILE *fpInTree);
char *ReadTreeToken(/*READ*/FILE *fp); /* returns a static array, or NULL on EOF */
void ReadTreeAddChild(int parent, int child, /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children);
/* Do not add the leaf if we already set this unique-set to another parent */
void ReadTreeMaybeAddLeaf(int parent, char *name,
			  hashstrings_t *hashnames, uniquify_t *unique,
			  /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children);
void ReadTreeRemove(/*IN/OUT*/int *parents, /*IN/OUT*/children_t *children, int node);

/* Routines to support tree traversal and prevent visiting a node >1 time
   (esp. if topology changes).
*/
typedef bool *traversal_t;
traversal_t InitTraversal(NJ_t*);
void SkipTraversalInto(int node, /*IN/OUT*/traversal_t traversal);
traversal_t FreeTraversal(traversal_t, NJ_t*); /*returns NULL*/

/* returns new node, or -1 if nothing left to do. Use root for the first call.
   Will return every node and then root.
   Uses postorder tree traversal (depth-first search going down to leaves first)
   Keeps track of which nodes are visited, so even after an NNI that swaps a
   visited child with an unvisited uncle, the next call will visit the
   was-uncle-now-child. (However, after SPR moves, there is no such guarantee.)

   If pUp is not NULL, then, if going "back up" through a previously visited node
   (presumably due to an NNI), then it will return the node another time,
   with *pUp = true.
*/
int TraversePostorder(int lastnode, NJ_t *NJ, /*IN/OUT*/traversal_t,
		      /*OUT OPTIONAL*/bool *pUp);

/* Routines to support storing up-profiles during tree traversal
   Eventually these should be smart enough to do weighted joins and
   to minimize memory usage
*/
profile_t **UpProfiles(NJ_t *NJ);
profile_t *GetUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node, bool useML);
profile_t *DeleteUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node); /* returns NULL */
profile_t **FreeUpProfiles(profile_t **upProfiles, NJ_t *NJ); /* returns NULL */

/* Recomputes the profile for a node, presumably to reflect topology changes
   If bionj is set, does a weighted join -- which requires using upProfiles
   If useML is set, computes the posterior probability instead of averaging
 */
void RecomputeProfile(/*IN/OUT*/NJ_t *NJ, /*IN/OUT*/profile_t **upProfiles, int node, bool useML);

/* Recompute profiles going up from the leaves, using the provided distance matrix
   and unweighted joins
*/
void RecomputeProfiles(/*IN/OUT*/NJ_t *NJ, /*OPTIONAL*/distance_matrix_t *dmat);

void RecomputeMLProfiles(/*IN/OUT*/NJ_t *NJ);

/* If bionj is set, computes the weight to be given to A when computing the
   profile for the ancestor of A and B. C and D are the other profiles in the quartet
   If bionj is not set, returns -1 (which means unweighted in AverageProfile).
   (A and B are the first two profiles in the array)
*/
double QuartetWeight(profile_t *profiles[4], distance_matrix_t *dmat, int nPos);

/* Returns a list of nodes, starting with node and ending with root */
int *PathToRoot(NJ_t *NJ, int node, /*OUT*/int *depth);
int *FreePath(int *path, NJ_t *NJ); /* returns NULL */

// Knuth code - random number generator
void ran_array(long aa[],int n);
long ran_arr_cycle();

void ReadTreeError(char *err, char *token);

#endif