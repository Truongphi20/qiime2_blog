#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

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