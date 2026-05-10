import biom
import skbio
import pandas as pd
import numpy as np

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

    # Equivalent to Cython DTYPE conversion
    counts = counts.astype(np.intp, copy=False)

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

    counts_by_node = _nodes_by_counts(counts, taxa, tree_index)
    branch_lengths = tree_index["length"]

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, tree_index, branch_lengths 

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/alpha/_pd.py:19
def _setup_pd(counts, taxa, tree, validate, rooted, single_sample):
    counts_by_node, _, branch_lengths = _vectorize_counts_and_tree(counts, taxa, tree)
    return counts_by_node, branch_lengths

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/alpha/_pd.py:53
def faith_pd(counts, taxa=None, tree=None, validate=True, otu_ids=None):
    
    counts_by_node, branch_lengths = _setup_pd(
        counts, taxa, tree, validate, rooted=True, single_sample=True
    )

    return (branch_lengths * (counts_by_node > 0)).sum()

def calculate_faith_pd(table: biom.Table, tree: skbio.TreeNode):
    
    # Convert to dense array: rows are ASVs, columns are Samples
    counts = table.matrix_data.toarray().T.astype(int) 
    sample_ids = table.ids(axis='sample')
    otu_ids = table.ids(axis='observation')

    # Calculate Faith's PD for each sample
    pd_results = []
    for i in range(len(sample_ids)):
        sample_counts = counts[i]
        
        # Calculate the metric
        result = faith_pd(sample_counts, otu_ids, tree)
        pd_results.append(result)

    return pd.Series(pd_results, index=sample_ids, name='faith_pd')

if __name__ == "__main__":
    BIOM_FILE = "/workspaces/qiime2_blog/support_data/feature-table-rarefied.biom"
    TREE_FILE = "/workspaces/qiime2_blog/support_data/tree.nwk"

    # Load the feature table
    table = biom.load_table(BIOM_FILE)

    # Load the phylogenetic tree
    tree = skbio.TreeNode.read(TREE_FILE)
    
    faith_pd_series = calculate_faith_pd(table, tree)
    print(faith_pd_series)

    # Expected output:
    # q2..L1S105     6.384621
    # q2..L1S140     5.900037
    # q2..L1S208     6.700060
    # q2..L1S257     7.362996
    # q2..L1S281     5.489152
    # q2..L1S57      5.572742
    # q2..L1S76      6.159797
    # q2..L1S8       4.973571
    # q2..L2S155    15.671622
    # q2..L2S175    13.210765
    # q2..L2S204    15.803438
    # q2..L2S222    18.029198
    # q2..L2S240     8.010033
    # q2..L2S309    13.595307
    # q2..L2S357     8.792001
    # q2..L2S382    12.776710
    # q2..L3S294     8.676022
    # q2..L3S313     7.818368
    # q2..L3S378     4.646372
    # q2..L4S112    15.468081
    # q2..L4S137    13.628922
    # q2..L4S63     21.770786
    # q2..L5S104     5.608632
    # q2..L5S155     4.361857
    # q2..L5S174     5.357451
    # q2..L5S203     5.903046
    # q2..L5S222     5.311506
    # q2..L5S240     5.060430
    # q2..L6S20      6.737120
    # q2..L6S68      6.374415
    # q2..L6S93      6.317089