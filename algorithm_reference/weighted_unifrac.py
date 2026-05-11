import biom
import skbio
from  scipy.spatial import distance
from skbio.diversity.beta import _unifrac
import pandas as pd
import numpy as np

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/beta/_unifrac.py:542
def _setup_multiple_weighted_unifrac(counts, taxa, tree, normalized, validate):
    counts_by_node, tree_index, branch_lengths = _unifrac._setup_multiple_unifrac(
        counts, taxa, tree, validate
    )
    tip_indices = np.array(
        [n.id for n in tree_index["id_index"].values() if n.is_tip()], dtype=np.intp
    )

    def f(u_node_counts, v_node_counts):
            u_total_count = np.take(u_node_counts, tip_indices).sum()
            v_total_count = np.take(v_node_counts, tip_indices).sum()
            u, _, _ = _unifrac._weighted_unifrac(
                u_node_counts,
                v_node_counts,
                u_total_count,
                v_total_count,
                branch_lengths,
            )
            return u
    
    return f, counts_by_node

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/scipy/spatial/distance.py:2627
def _pdist_callable(X, metric, **kwargs):
    n = X.shape[0]
    out_size = (n * (n - 1)) // 2
    dm = np.empty((out_size,), dtype=np.float64)
    k = 0
    for i in range(X.shape[0] - 1):
        for j in range(i + 1, X.shape[0]):
            dm[k] = metric(X[i], X[j], **kwargs)
            k += 1
    return dm

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(
    metric, counts, taxa, tree, ids=None, validate=True, pairwise_func=None, **kwargs
):
    metric, counts_by_node = _setup_multiple_weighted_unifrac(
            counts, taxa=taxa, tree=tree, normalized=False, validate=validate
        )
    counts = counts_by_node
    distances = _pdist_callable(counts, metric=metric, **kwargs)

    data = distance.squareform(distances, force="tomatrix", checks=False)

    return pd.DataFrame(data, columns=ids, index=ids)

def calculate_weighted_unifrac(table, tree):
    """
    Calculates the weighted UniFrac distance matrix.
    """
    # Prepare Table Data
    counts = table.matrix_data.toarray().T.astype(int)
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # Calculate Distance Matrix
    dm = beta_diversity(
        metric='weighted_unifrac',
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
    weighted_unifrac_df = calculate_weighted_unifrac(table, tree)
    
    print(weighted_unifrac_df.iloc[:5, :5])

# Expected output:
#             q2..L1S105  q2..L1S140  q2..L1S208  q2..L1S257  q2..L1S281
# q2..L1S105    0.000000    0.105481    0.101449    0.129749    0.137230
# q2..L1S140    0.105481    0.000000    0.103910    0.151287    0.156065
# q2..L1S208    0.101449    0.103910    0.000000    0.066404    0.083570
# q2..L1S257    0.129749    0.151287    0.066404    0.000000    0.080584
# q2..L1S281    0.137230    0.156065    0.083570    0.080584    0.000000