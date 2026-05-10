import biom
import skbio
from skbio.diversity.beta import _unifrac
import pandas as pd
import numpy as np
import functools

import sys
sys.path.extend(["/workspaces/qiime2_blog/commands/scipy_7dcd8c5_src", "/workspaces/qiime2_blog/commands/scikit-bio-0.6.2"])
import _distance_wrap # type: ignore
import _phylogenetic # type: ignore

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

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_util.py:186
def _vectorize_counts_and_tree(counts, taxa, tree):
    tree_index = tree.to_array(nan_length_value=0.0)
    taxa = np.asarray(taxa)
    counts = np.atleast_2d(counts)
    counts_by_node = _phylogenetic._nodes_by_counts(counts, taxa, tree_index)
    branch_lengths = tree_index["length"]

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, tree_index, branch_lengths

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/beta/_unifrac.py:495
def _setup_multiple_unifrac(counts, taxa, tree):

    counts_by_node, tree_index, branch_lengths = _vectorize_counts_and_tree(
        counts, taxa, tree
    )

    return counts_by_node, tree_index, branch_lengths

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/beta/_unifrac.py:506
def _setup_multiple_unweighted_unifrac(counts, taxa, tree):
    counts_by_node, _, branch_lengths = _setup_multiple_unifrac(
        counts, taxa, tree
    )

    f = functools.partial(_unifrac._unweighted_unifrac, branch_lengths=branch_lengths)
    return f, counts_by_node 

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(metric, counts, ids, taxa, tree):
    
    metric, counts_by_node = _setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree
        )
    counts = counts_by_node
    
    distances = _pdist_callable(counts, metric=metric)
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