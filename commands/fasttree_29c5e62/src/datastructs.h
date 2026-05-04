#ifndef DATASTRUCTS_H 
#define DATASTRUCTS_H

#ifdef USE_DOUBLE
#define SSE_STRING "Double precision"
typedef double numeric_t;
#define ScanNumericSpec "%lf"
#else
typedef float numeric_t;
#define ScanNumericSpec "%f"
#endif

#define MAXCODES 20
#define NOCODE 127

typedef struct {
  int nPos;
  int nSeq;
  char **names;
  char **seqs;
  int nSaved; /* actual allocated size of names and seqs */
} alignment_t;

/* For each position in a profile, we have a weight (% non-gapped) and a
   frequency vector. (If using a matrix, the frequency vector is in eigenspace).
   We also store codes for simple profile positions (all gaps or only 1 value)
   If weight[pos] > 0 && codes[pos] == NOCODE then we store the vector
   vectors itself is sets of nCodes long, so the vector for the ith nonconstant position
   starts at &vectors[nCodes*i]
   
   To speed up comparison of outprofile to a sequence or other simple profile, we also
   (for outprofiles) store codeDist[iPos*nCodes+k] = dist(k,profile[iPos])

   For constraints, we store a vector of nOn and nOff
   If not using constraints, those will be NULL
*/
typedef struct {
  /* alignment profile */
  numeric_t *weights;
  unsigned char *codes;
  numeric_t *vectors;		/* NULL if no non-constant positions, e.g. for leaves */
  int nVectors;
  numeric_t *codeDist;		/* Optional -- distance to each code at each position */

  /* constraint profile */
  int *nOn;
  int *nOff;
} profile_t;

/* A visible node is a pair of nodes i, j such that j is the best hit of i,
   using the neighbor-joining criterion, at the time the comparison was made,
   or approximately so since then.

   Note that variance = dist because in BIONJ, constant factors of variance do not matter,
   and because we weight ungapped sequences higher naturally when averaging profiles,
   so we do not take this into account in the computation of "lambda" for BIONJ.

   For the top-hit list heuristic, if the top hit list becomes "too short",
   we store invalid entries with i=j=-1 and dist/criterion very high.
*/
typedef struct {
  int i, j;
  numeric_t weight;			/* Total product of weights (maximum value is nPos)
				   This is needed for weighted joins and for pseudocounts,
				   but not in most other places.
				   For example, it is not maintained by the top hits code */
  numeric_t dist;			/* The uncorrected distance (includes diameter correction) */
  numeric_t criterion;		/* changes when we update the out-profile or change nActive */
} besthit_t;

typedef struct {
  int nChild;
  int child[3];
} children_t;

typedef struct {
  /* Distances between amino acids */
  numeric_t distances[MAXCODES][MAXCODES];

  /* Inverse of the eigenvalue matrix, for rotating a frequency vector
     into eigenspace so that profile similarity computations are
     O(alphabet) not O(alphabet*alphabet) time.
  */
  numeric_t eigeninv[MAXCODES][MAXCODES];
  numeric_t eigenval[MAXCODES];	/* eigenvalues */


  /* eigentot=eigeninv times the all-1s frequency vector
     useful for normalizing rotated frequency vectors
  */
  numeric_t eigentot[MAXCODES];	

  /* codeFreq is the transpose of the eigeninv matrix is
     the rotated frequency vector for each code */
  numeric_t codeFreq[MAXCODES][MAXCODES];
  numeric_t gapFreq[MAXCODES];
} distance_matrix_t;


/* A transition matrix gives the instantaneous rate of change of frequencies
   df/dt = M . f
   which is solved by
   f(t) = exp(M) . f(0)
   and which is not a symmetric matrix because of
   non-uniform stationary frequencies stat, so that
   M stat = 0
   M(i,j) is instantaneous rate of j -> i, not of i -> j

   S = diag(sqrt(stat)) is a correction so that
   M' = S**-1 M S is symmetric
   Let W L W**-1 = M' be an eigendecomposition of M'
   Because M' is symmetric, W can be a rotation, and W**-1 = t(W)
   Set V = S*W
   M = V L V**-1 is an eigendecomposition of M
   Note V**-1 = W**-1 S**-1 = t(W) S**-1
   
   Evolution by time t is given by

   exp(M*t) = V exp(L*t) V**-1
   P(A & B | t) = B . exp(M*t) . (A * stat)
   note this is *not* the same as P(A->B | t)

   and we can reduce some of the computations from O(a**2) to O(a) time,
   where a is the alphabet size, by storing frequency vectors as
   t(V) . f = t(W) . t(S) . f

   Then
   P(f0 & f1 | t) = f1 . exp(M*t) . f0 * (f0 . stat) = sum(r0j * r1j * exp(l_j*t))
   where r0 and r1 are the transformed vectors

   Posterior distribution of P given children f0 and f1 is given by
   P(i | f0, f1, t0, t1) = stat * P(i->f0 | t0) * P(i->f1 | t1)
   = P(i & f0 | t0) * P(i & f1 | t1) / stat
   ~ (V . exp(t0*L) . r0) * (V . exp(t1*L) . r1) / stat

   When normalize this posterior distribution (to sum to 1), divide by stat,
   and transform by t(V) -- this is the "profile" of internal nodes

   To eliminate the O(N**2) step of transforming by t(V), if the posterior
   distribution of an amino acid is near 1 then we can approximate it by
   P(i) ~= (i==A) * w + nearP(i) * (1-w), where
   w is fit so that P(i==A) is correct
   nearP = Posterior(i | i, i, 0.1, 0.1) [0.1 is an arbitrary choice]
   and we confirm that the approximation works well before we use it.

   Given this parameter w we can set
   rotated_posterior = rotation(w * (i==A)/stat + (1-w) * nearP/stat)
   = codeFreq(A) * w/stat(A) + nearFreq(A) * (1-w)
 */
typedef struct {
  numeric_t stat[MAXCODES]; /* The stationary distribution */
  numeric_t statinv[MAXCODES];	/* 1/stat */
  /* the eigenmatrix, with the eigenvectors as columns and rotations of individual
     characters as rows. Also includes a NOCODE entry for gaps */
  numeric_t codeFreq[NOCODE+1][MAXCODES];
  numeric_t eigeninv[MAXCODES][MAXCODES]; /* Inverse of eigenmatrix */
  numeric_t eigeninvT[MAXCODES][MAXCODES]; /* transpose of eigeninv */
  numeric_t eigenval[MAXCODES];	/* Eigenvalues  */
  /* These are for approximate posteriors (off by default) */
  numeric_t nearP[MAXCODES][MAXCODES]; /* nearP[i][j] = P(parent=j | both children are i, both lengths are 0.1 */
  numeric_t nearFreq[MAXCODES][MAXCODES]; /* rotation of nearP/stat */
} transition_matrix_t;

typedef struct {
  int nRateCategories;
  numeric_t *rates;			/* 1 per rate category */
  unsigned int *ratecat;	/* 1 category per position */
} rates_t;

typedef struct {
  /* The input */
  int nSeq;
  int nPos;
  char **seqs;			/* the aligment sequences array (not reallocated) */
  distance_matrix_t *distance_matrix; /* a pointer (not reallocated), or NULL if using %identity distance */
  transition_matrix_t *transmat; /* a pointer (is allocated), or NULL for Jukes-Cantor */
  /* Topological constraints are represented for each sequence as binary characters
     with values of '0', '1', or '-' (for missing data)
     Sequences that have no constraint may have a NULL string
  */
  int nConstraints;
  char **constraintSeqs;

  /* The profile data structures */
  int maxnode;			/* The next index to allocate */
  int maxnodes;			/* Space allocated in data structures below */
  profile_t **profiles;         /* Profiles of leaves and intermediate nodes */
  numeric_t *diameter;		/* To correct for distance "up" from children (if any) */
  numeric_t *varDiameter;		/* To correct variances for distance "up" */
  numeric_t *selfdist;		/* Saved for use in some formulas */
  numeric_t *selfweight;		/* Saved for use in some formulas */

  /* Average profile of all active nodes, the "outprofile"
   * If all inputs are ungapped, this has weight 1 (not nSequences) at each position
   * The frequencies all sum to one (or that is implied by the eigen-representation)
   */
  profile_t *outprofile;
  double totdiam;

  /* We sometimes use stale out-distances, so we remember what nActive was  */
  numeric_t *outDistances;		/* Sum of distances to other active (parent==-1) nodes */
  int *nOutDistActive;		/* What nActive was when this outDistance was computed */

  /* the inferred tree */
  int root;			/* index of the root. Unlike other internal nodes, it has 3 children */
  int *parent;			/* -1 or index of parent */
  children_t *child;
  numeric_t *branchlength;		/* Distance to parent */
  numeric_t *support;		/* 1 for high-confidence nodes */

  /* auxilliary data for maximum likelihood (defaults to 1 category of rate=1.0) */
  rates_t rates;
} NJ_t;

/* Uniquify sequences in an alignment -- map from indices
   in the alignment to unique indicies in a NJ_t
*/
typedef struct {
  int nSeq;
  int nUnique;
  int *uniqueFirst;		/* iUnique -> iAln */
  int *alnNext;			/* iAln -> next, or -1  */
  int *alnToUniq;		/* iAln -> iUnique, or -1 if another was the exemplar */
  char **uniqueSeq;		/* indexed by iUniq -- points to strings allocated elsewhere */
} uniquify_t;

/* Describes which switch to do */
typedef enum {ABvsCD,ACvsBD,ADvsBC} nni_t;

/* A list of these describes a chain of NNI moves in a rooted tree,
   making up, in total, an SPR move
*/
typedef struct {
  int nodes[2];
  double deltaLength;		/* change in tree length for this step (lower is better) */
} spr_step_t;

/* Keep track of hits for the top-hits heuristic without wasting memory
   j = -1 means empty
   If j is an inactive node, this may be replaced by that node's parent (and dist recomputed)
 */
typedef struct {
  int j;
  numeric_t dist;
} hit_t;

typedef struct {
  int nHits;			/* the allocated and desired size; some of them may be empty */
  hit_t *hits;
  int hitSource;		/* where to refresh hits from if a 2nd-level top-hit list, or -1 */
  int age;			/* number of joins since a refresh */
} top_hits_list_t;

typedef struct {
  int m;			 /* size of a full top hits list, usually sqrt(N) */
  int q;			 /* size of a 2nd-level top hits, usually sqrt(m) */
  int maxnodes;
  top_hits_list_t *top_hits_lists; /* one per node */
  hit_t *visible;		/* the "visible" (very best) hit for each node */

  /* The top-visible set is a subset, usually of size m, of the visible set --
     it is the set of joins to select from
     Each entry is either a node whose visible set entry has a good (low) criterion,
     or -1 for empty, or is an obsolete node (which is effectively the same).
     Whenever we update the visible set, should also call UpdateTopVisible()
     which ensures that none of the topvisible set are stale (that is, they
     all point to an active node).
  */
  int nTopVisible;		/* nTopVisible = m * topvisibleMult */
  int *topvisible;

  int topvisibleAge;		/* joins since the top-visible list was recomputed */

#ifdef OPENMP
  /* 1 lock to read or write any top hits list, no thread grabs more than one */
  omp_lock_t *locks;
#endif
} top_hits_t;
#endif DATASTRUCTS_H