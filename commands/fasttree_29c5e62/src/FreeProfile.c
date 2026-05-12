#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

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