import biom
import skbio
from skbio.diversity import beta_diversity
from skbio.diversity.beta import _unifrac
from skbio.stats.distance import DistanceMatrix
import scipy
import pandas as pd
import numpy as np

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_util.py:224
def _get_phylogenetic_kwargs(counts, **kwargs):
    taxa = kwargs.pop("taxa")
    tree = kwargs.pop("tree")
    return taxa, tree, kwargs

# /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/skbio/diversity/_driver.py:367
def beta_diversity(
    metric, counts, ids=None, validate=True, pairwise_func=None, **kwargs
):
    taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
    metric, counts_by_node = _unifrac._setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
    counts = counts_by_node
    pairwise_func = scipy.spatial.distance.pdist
    distances = pairwise_func(counts, metric=metric, **kwargs)
    return DistanceMatrix(distances, ids)

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
    return dm.to_data_frame()

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