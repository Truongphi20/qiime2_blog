#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void ReadConstraints(
    alignment_t **constraints, 
    char ***uniqConstraints, 
    hashstrings_t **hashnames, 
    char **constraintsFile,
    FILE **fpConstraints,
    bool bQuote,
    int iAln,
    uniquify_t **unique 
)
{
    if (*constraintsFile != NULL) {
        *constraints = ReadAlignment(*fpConstraints, bQuote);
        if ((*constraints)->nSeq < 4) {
            fprintf(stderr, "Warning: constraints file with less than 4 sequences ignored:\nalignment #%d in %s\n",
            iAln+1, *constraintsFile);
            *constraints = FreeAlignment(*constraints);
        } else {
            *uniqConstraints = AlnToConstraints(*constraints, *unique, *hashnames);
            ProgressReport("Read the constraints",0,0,0,0);
        }
    }	/* end load constraints */
}