#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void MakeMatrix(alignment_t *aln, distance_matrix_t *distance_matrix)
{
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
}