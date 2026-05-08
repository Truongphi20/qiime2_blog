import biom
import pandas as pd
import numpy as np
import skbio
from skbio.stats.distance import DistanceMatrix
from scipy.spatial.distance import _METRIC_ALIAS
from scipy._lib._util import _asarray_validated

import sys
sys.path.insert(0, "/workspaces/qiime2_blog/commands/scipy_7dcd8c5_src")

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/scipy/spatial/distance.py:1864
def pdist(X, metric='euclidean', *, out=None, **kwargs):
    X = _asarray_validated(X, sparse_ok=False, objects_ok=True, mask_ok=True,
                           check_finite=False)
    s = X.shape
    m, n = s
    mstr = metric.lower()
    metric_info = _METRIC_ALIAS.get(mstr, None)
    pdist_fn = metric_info.pdist_func
    return pdist_fn(X, out=out, **kwargs)

# commands/scipy_7dcd8c5_src/distance_impl.h:706
def dist_to_squareform_from_vector_double(M_flat, X, d):
    # it1 is the start of the horizontal fill (diagonal + 1)
    # v is the index into our source vector X
    v = 0
    
    for i in range(1, d):
        # 1. The Horizontal Fill (memcpy equivalent)
        # C code: memcpy(it1, v, (n - i) * sizeof(double))
        it1 = (i - 1) * d + i 
        length = d - i
        M_flat[it1 : it1 + length] = X[v : v + length]
        
        # 2. The Vertical Fill (Symmetry)
        # C code: it2 = M + i * (n + 1) - 1
        it2 = i * (d + 1) - 1
        
        for j in range(i, d):
            M_flat[it2] = X[v] # Mirror the single value
            v += 1             # Move the source pointer
            it2 += d           # Jump down one row (same column)

    return M_flat.reshape(d, d)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/scipy/spatial/distance.py:2196
def squareform(X, force="no", checks=True):
    X = np.ascontiguousarray(X)
    s = X.shape

    # Grab the closest value to the square root of the number
    # of elements times 2 to see if the number of elements
    # is indeed a binomial coefficient.
    d = int(np.ceil(np.sqrt(s[0] * 2)))

    # Allocate memory for the distance matrix.
    M = np.zeros(d**2, dtype=X.dtype)

    # Fill in the values of the distance matrix.
    M = dist_to_squareform_from_vector_double(M, X, d)

    # Return the distance matrix.
    return M


# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/sklearn/metrics/pairwise.py:2168
def pairwise_distances(X, Y=None, metric="euclidean", *, n_jobs=None, force_all_finite=True, **kwds):
    # precompute data-derived metric params
    params = {}
    kwds.update(**params)
    return squareform(pdist(X, metric=metric, **kwds))

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(metric, counts, ids=None, pairwise_func=None):
    counts = (counts > 0.0)
    distances = pairwise_func(counts, metric=metric)
    return DistanceMatrix(distances, ids)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/q2_diversity_lib/beta.py:187
def jaccard(table: biom.Table, n_jobs: int = 1) -> skbio.DistanceMatrix:
    counts = table.matrix_data.toarray().T.copy()
    sample_ids = table.ids(axis='sample')
    jaccard_table = beta_diversity(
            metric='jaccard', 
            counts=counts, 
            ids=sample_ids,
            pairwise_func=pairwise_distances
    )
    return jaccard_table.to_data_frame()

if __name__ == "__main__":
    # Mini table (10 ASVs x 11 samples) extracted from the "table.qza" of DADA2 output 
    table_path = "/workspaces/qiime2_blog/support_data/mini-feature-table.tsv"

    # Read table
    with open(table_path, 'r') as f:
        table = biom.Table.from_tsv(f, None, None, lambda x: x)

    # Compute
    results = jaccard(table)
    print(results)

# Expected output:
#          L1S105    L1S140    L1S208    L1S257    L1S281     L1S57     L1S76      L1S8    L2S155    L2S175    L2S204
# L1S105  0.000000  0.600000  0.714286  0.666667  0.600000  0.250000  0.250000  0.500000  0.857143  0.875000  0.666667
# L1S140  0.600000  0.000000  0.400000  0.250000  0.000000  0.500000  0.500000  0.666667  1.000000  1.000000  1.000000
# L1S208  0.714286  0.400000  0.000000  0.200000  0.400000  0.666667  0.666667  0.571429  0.714286  0.750000  0.875000
# L1S257  0.666667  0.250000  0.200000  0.000000  0.250000  0.600000  0.600000  0.714286  0.857143  0.875000  0.857143
# L1S281  0.600000  0.000000  0.400000  0.250000  0.000000  0.500000  0.500000  0.666667  1.000000  1.000000  1.000000
# L1S57   0.250000  0.500000  0.666667  0.600000  0.500000  0.000000  0.000000  0.400000  0.833333  0.857143  0.833333
# L1S76   0.250000  0.500000  0.666667  0.600000  0.500000  0.000000  0.000000  0.400000  0.833333  0.857143  0.833333
# L1S8    0.500000  0.666667  0.571429  0.714286  0.666667  0.400000  0.400000  0.000000  0.714286  0.571429  0.875000
# L2S155  0.857143  1.000000  0.714286  0.857143  1.000000  0.833333  0.833333  0.714286  0.000000  0.200000  0.400000
# L2S175  0.875000  1.000000  0.750000  0.875000  1.000000  0.857143  0.857143  0.571429  0.200000  0.000000  0.500000
# L2S204  0.666667  1.000000  0.875000  0.857143  1.000000  0.833333  0.833333  0.875000  0.400000  0.500000  0.000000