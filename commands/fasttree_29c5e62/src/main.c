#include "include_stuff.h"
#include "datastructs.h"
#include "hyper_parameters.h"
#include "support_functions.h"


int main(int argc, char **argv) {
    FastTreeOptions_t opt;
    InitOptions(&opt);
    ParseCommandLine(argc, argv, &opt);
    
    FastTree(opt);
}