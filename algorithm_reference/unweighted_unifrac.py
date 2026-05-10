import biom
import skbio
from skbio.diversity import beta_diversity
from skbio.diversity.beta import _unifrac
import pandas as pd
import numpy as np

import sys
sys.path.insert(0, "/workspaces/qiime2_blog/commands/scipy_7dcd8c5_src")
import _distance_wrap # type: ignore

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/stats/distance/_base.py:47
def squareform(X):
    s = X.shape

    # Grab the closest value to the square root of the number
    # of elements times 2 to see if the number of elements
    # is indeed a binomial coefficient.
    d = int(np.ceil(np.sqrt(s[0] * 2)))

    # Allocate memory for the distance matrix.
    M = np.zeros((d, d), dtype=X.dtype)

    # Fill in the values of the distance matrix.
    _distance_wrap.to_squareform_from_vector_wrap(M, X)

    return M

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/scipy/spatial/distance.py:2627
def _pdist_callable(X, metric):
    n = X.shape[0]
    out_size = (n * (n - 1)) // 2
    dm = np.empty((out_size,), dtype=np.float64)
    k = 0
    for i in range(X.shape[0] - 1):
        for j in range(i + 1, X.shape[0]):
            dm[k] = metric(X[i], X[j])
            k += 1
    return dm


# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/scipy/spatial/distance.py:1864
def pdist(X, metric='euclidean'):
    s = X.shape
    m, n = s

    return _pdist_callable(X, metric=metric)

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(
    metric, counts, ids, taxa, tree, validate=True
):
    metric, counts_by_node = _unifrac._setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
    counts = counts_by_node
    
    distances = pdist(counts, metric=metric)
    data = squareform(distances)

    return pd.DataFrame(data, columns=ids, index=ids)

def calculate_unweighted_unifrac(table, tree):
    """
    Calculates the Unweighted UniFrac distance matrix.
    """
    # Prepare Table Data
    counts = table.matrix_data.toarray().T.astype(int)
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # Calculate Distance Matrix
    dm = beta_diversity(
        metric='unweighted_unifrac',
        counts=counts,
        ids=sample_ids,
        taxa=otu_ids,
        tree=tree
    )

    # Convert to a readable DataFrame
    return dm

if __name__ == "__main__":
    BIOM_FILE = "/workspaces/qiime2_blog/support_data/feature-table-rarefied.biom"
    TREE_FILE = "/workspaces/qiime2_blog/support_data/tree.nwk"

    # Load data
    table = biom.load_table(BIOM_FILE)
    tree = skbio.TreeNode.read(TREE_FILE)
    
    # Calculate
    unweighted_unifrac_df = calculate_unweighted_unifrac(table, tree)
    
    print(unweighted_unifrac_df.iloc[:5, :5])

# Expected output:
#             q2..L1S105  q2..L1S140  q2..L1S208  q2..L1S257  q2..L1S281
# q2..L1S105    0.000000    0.441940    0.545205    0.540621    0.575789
# q2..L1S140    0.441940    0.000000    0.523709    0.534961    0.461661
# q2..L1S208    0.545205    0.523709    0.000000    0.346487    0.320912
# q2..L1S257    0.540621    0.534961    0.346487    0.000000    0.393405
# q2..L1S281    0.575789    0.461661    0.320912    0.393405    0.000000