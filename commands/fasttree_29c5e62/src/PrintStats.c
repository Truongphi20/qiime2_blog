#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

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
)
{
    for (int i = 0; i < nFPs; i++) {
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
}