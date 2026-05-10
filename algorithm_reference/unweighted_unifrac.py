import biom
import skbio
import pandas as pd
import numpy as np
import functools

import sys
sys.path.extend(["/workspaces/qiime2_blog/commands/scipy_7dcd8c5_src"])
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

def _traverse_reduce(child_index, count_array):
    """
    Propagate counts from children to parents.
    """

    for row in child_index:
        parent = row[0]

        for child in row[1:]:
            if child != -1:
                count_array[parent] += count_array[child]

# commands/scikit-bio-0.6.2/_phylogenetic.pyx:142
def _nodes_by_counts(counts, tip_ids, indexed):
    """Construct the count array and propagate counts up the tree.
    """

    nodes = indexed['name']

    # Allow counts to be a vector
    counts = np.atleast_2d(counts)

    # Determine observed IDs
    observed_indices = counts.sum(axis=0).nonzero()[0]
    observed_ids = tip_ids[observed_indices]
    observed_ids_set = set(observed_ids)

    # Map observed node names to their positions
    node_lookup = {}
    for i in range(nodes.shape[0]):
        n = nodes[i]
        if n in observed_ids_set:
            node_lookup[n] = i

    # Determine positions of observed IDs in nodes
    taxa_in_nodes = np.zeros(observed_ids.shape[0], dtype=np.intp)

    for i in range(observed_ids.shape[0]):
        n = observed_ids[i]
        taxa_in_nodes[i] = node_lookup[n]

    # count_array has:
    # rows = nodes
    # cols = environments
    n_count_vectors = counts.shape[0]

    count_array = np.zeros(
        (nodes.shape[0], n_count_vectors),
        dtype=np.intp
    )

    # Populate counts
    counts_t = counts.transpose()
    n_count_taxa = taxa_in_nodes.shape[0]

    for i in range(n_count_taxa):
        for j in range(n_count_vectors):
            count_array[taxa_in_nodes[i], j] = (
                counts_t[observed_indices[i], j]
            )

    # Propagate counts up the tree
    child_index = indexed['child_index'].astype(np.intp, copy=False)

    _traverse_reduce(child_index, count_array)

    return count_array

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_util.py:186
def _vectorize_counts_and_tree(counts, taxa, tree):
    tree_index = tree.to_array(nan_length_value=0.0)
    taxa = np.asarray(taxa)
    counts = np.atleast_2d(counts)
    counts_by_node = _nodes_by_counts(counts, taxa, tree_index)
    branch_lengths = tree_index["length"]

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, branch_lengths

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/beta/_unifrac.py:392
def _unweighted_unifrac(u_node_counts, v_node_counts, branch_lengths):

    unique_nodes = np.logical_xor(u_node_counts, v_node_counts)
    observed_nodes = np.logical_or(u_node_counts, v_node_counts)
    unique_branch_length = (branch_lengths * unique_nodes).sum()
    observed_branch_length = (branch_lengths * observed_nodes).sum()

    return unique_branch_length / observed_branch_length

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(counts, ids, taxa, tree):
    
    counts_by_node, branch_lengths = _vectorize_counts_and_tree(counts, taxa, tree)

    counts = counts_by_node
    metric = functools.partial(_unweighted_unifrac, branch_lengths=branch_lengths)
    
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