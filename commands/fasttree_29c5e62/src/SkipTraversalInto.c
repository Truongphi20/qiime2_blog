#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

void SkipTraversalInto(int node, /*IN/OUT*/traversal_t traversal) {
  traversal[node] = true;
}